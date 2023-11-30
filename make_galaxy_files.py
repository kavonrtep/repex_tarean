#!/usr/bin/env python
'''
Creates xml files for galaxy
'''

import subprocess
import urllib.request
import tempfile
import zipfile
import os
import glob

def parse_xml(template, output, keys):
    with open(template, 'r') as f:
        tmpl = f.read()
    for i in keys:
        tmpl = tmpl.replace(i, keys[i])
    with open(output, 'w') as f:
        f.write(tmpl)


long_hash = subprocess.check_output(
    "git log --pretty=format:'%H' -n 1  ", shell=True).decode('ascii').strip()
tag = subprocess.check_output(
    "git describe --tags --abbrev=0", shell=True).decode('ascii').strip()
galaxy_tool_version = "2" + tag[1:]

print("commit hash:", long_hash)

# write_version to file version_info.txt
subprocess.check_output("./get_version.sh")
print("\nUpdated version_info.txt :")
with open('version_info.txt') as f:
    for l in f:
        print(l)
version_info = l.strip()

print("\nGalaxy tool version:", galaxy_tool_version)

download_link = "https://bitbucket.org/petrnovak/repex_tarean/get/{}.zip".format(
    long_hash[0:12])
dir_name = "petrnovak-repex_tarean-{}".format(long_hash[0:12])
print("\nDownload link:\n", download_link, "\n")

print("Validation of download link")
filedata = urllib.request.urlopen(download_link)
print("  Downloading..")
data2write = filedata.read()
tmpzip = tempfile.NamedTemporaryFile(mode='wb')
with open(tmpzip.name, 'wb') as f:
    f.write(data2write)
print("  Download ok")

print("\nExtracting:")
tmpdir = tempfile.TemporaryDirectory()
with zipfile.ZipFile(tmpzip.name, 'r') as f:
    f.extractall(tmpdir.name)
print("  zip extracted.")
# check if dir exist:
if os.path.exists(tmpdir.name + "/" + dir_name):
    print("  directory {} OK".format(dir_name))


# remove old files:
files_to_remove = glob.glob("galaxy/tools/*") + glob.glob("galaxy/package/*") + glob.glob("galaxy/tools-cerit/*") 
# will delete all files except hiden (e.g. .shed.yml)
for i in files_to_remove:
    os.remove(i)


# create xml files
template = "./galaxy/templates/tool_dependencies_package_template.xml"
xml_out = "./galaxy/package/tool_dependencies.xml"
replacement = {
    '__DOWNLOAD_LINK__': download_link,
    '__DIR_NAME__': dir_name,
    '__VERSION_INFO__': version_info,
    '__GALAXY_TOOL_VERSION__': galaxy_tool_version,
    '__GALAXY_TOOL_VERSION_UNDERSCORE__': galaxy_tool_version.replace(".", "_"),
    '__CERIT_SPECIFIC__': "",
}
parse_xml("./galaxy/templates/tool_dependencies_package_template.xml",
          "./galaxy/package/tool_dependencies.xml", replacement)

parse_xml("./galaxy/templates/tool_dependencies_tool_template.xml",
          "./galaxy/tools/tool_dependencies.xml", replacement)

parse_xml("./galaxy/templates/repex_tarean_template.xml",
          "./galaxy/tools/repex_tarean.xml", replacement)

parse_xml("./galaxy/templates/repex_full_clustering_template.xml",
          "./galaxy/tools/repex_full_clustering.xml", replacement)

# build galaxy tarball - create package.tar.gz
subprocess.check_call("planemo shed_build galaxy/package", shell=True)
new_name = "galaxy/package_repex_tarean_{}_{}.tar.gz".format(
    galaxy_tool_version.replace(".", "_"),
    version_info.replace(" ", ""))
os.rename("galaxy/package.tar.gz", new_name)
print("package.tar gz renamed to ", new_name)


# build galaxy tarball - create tools.tar.gz
subprocess.check_call("planemo shed_build galaxy/tools", shell=True)
new_name = "galaxy/repeatexplorer2_{}.tar.gz".format(
              version_info.replace(" ", ""))
os.rename("galaxy/tools.tar.gz", new_name)
print("package.tar gz renamed to ", new_name)



# build packages for cerit galaxy
with open("galaxy/templates/cerit_specific_options.xml", 'r') as f:
        cerit_specific_options = f.read()
replacement["__CERIT_SPECIFIC__"] = cerit_specific_options
parse_xml("./galaxy/templates/repex_tarean_template.xml",
          "./galaxy/tools-cerit/repex_tarean.xml", replacement)

parse_xml("./galaxy/templates/repex_full_clustering_template.xml",
          "./galaxy/tools-cerit/repex_full_clustering.xml", replacement)

parse_xml("./galaxy/templates/tool_dependencies_tool_template.xml",
          "./galaxy/tools-cerit/tool_dependencies.xml", replacement)

# build galaxy tarball - create tools.tar.gz
subprocess.check_call("planemo shed_build galaxy/tools-cerit", shell=True)
new_name = "galaxy/repeatexplorer2_cerit_{}.tar.gz".format(
              version_info.replace(" ", ""))
os.rename("galaxy/tools-cerit.tar.gz", new_name)
print("package.tar gz renamed to ", new_name)

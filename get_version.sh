#!/bin/sh
# this file is copied into ./git/hooks.post-commit and post-checkout
branch=$(git rev-parse --abbrev-ref HEAD)
shorthash=$(git log --pretty=format:'%h' -n 1)
revcount=$(git log --oneline | wc -l)
tag=$(git describe --tags --abbrev=0)
echo  "version:" ${tag}"-"${revcount}"("$shorthash") branch:" $branch > version_info.txt


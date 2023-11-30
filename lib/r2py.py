#!/usr/bin/env python3
import os
import atexit
import socket
import time
import config
import pyRserve

def shutdown(port):
    try:
        conn = pyRserve.connect(port=port)
        print("Shutting down Rserv...", end="")
        conn.shutdown()
        print("Done")
    except pyRserve.rexceptions.RConnectionRefused:
        print("connection to Rserve refused, server is probably already down")

def get_open_port():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("", 0))
    s.listen(1)
    port = s.getsockname()[1]
    s.close()
    return port

    # find free port
def create_connection():
    '''Start R Rserver and test connection, port is
    stored in config.RSERVE_PORT
    '''
    config.RSERVE_PORT = get_open_port()
    print('Trying to start Rserve...',)
    os.system(
        "R CMD Rserve --RS-port {} -q --no-save ".format(config.RSERVE_PORT))
        # wait for server to start accepting connections
    time.sleep(1)
    try:
        conn = pyRserve.connect(port=config.RSERVE_PORT)
        print("connection OK")
        conn.close()
        atexit.register(shutdown, config.RSERVE_PORT)
        return config.RSERVE_PORT
    except:
        print("Connection with Rserve was not established!")
        raise


class setFunctionName():
    # decorator

    def __init__(self, f, name):
        self.f = f
        self.name = name

    def __call__(self, *args, **kwargs):
        return self.f(self.name, *args, **kwargs)


def convert_types(fn):
    ''' decorator function to convert type for r2py'''
    allowed_classes = [str, int, float, list, bool, type(None)]
    # everything else is converted to str

    def fn_wrapper(*args, **kwargs):
        new_args = list(args)
        new_kwargs = kwargs
        for i, value in enumerate(args):
            if any(type(value) is i for i in allowed_classes):
                new_args[i] = value
            else:
                new_args[i] = str(value)
        for i, value in kwargs.items():
            if any(type(value) is i for i in allowed_classes):
                new_kwargs[i] = value
            else:
                new_kwargs[i] = str(value)
        return fn(*new_args, **new_kwargs)

    return fn_wrapper


class R():

    def __init__(self, source, verbose=False):
        ''' Code in file should defined R functions which will be linked to python function
        purpose of this is to make it memory efficient - rserve connection will be closed
        after every exetion so memory is released.
        warning  - Source code is executed each time function is used so it should only
        contain function definition!!

        '''
        self.source = os.path.realpath(source)
        conn = pyRserve.connect(port=config.RSERVE_PORT)
        conn.voidEval("source('{}', chdir=TRUE)".format(self.source))
        # if single object is define then fn return str! conversion neccessary
        object_names = list(conn.r.ls())
        if verbose:
            print("R function loaded:", end=" ")
        for i in object_names:
            ## skip these objects - not compatible with older version of rserve
            ## and not needed to be accessed from python
            if i in ['DT_OPTIONS', 'HTMLHEADER', 'WD', 'options',
                     'htmlheader', 'options', 'xcolor_code', 'TANDEM_RANKS', 'RANKS_TANDEM',]:
                continue
            try:
                obj = getattr(conn.r, i)
                if isinstance(obj, pyRserve.rconn.RFuncProxy):
                    if verbose:
                        print(i, end=" ")
                    @convert_types
                    def rwrapper(fname, *args, **kwargs):
                        c = pyRserve.connect(port=config.RSERVE_PORT)
                        c.r.setwd(os.getcwd())
                        c.voidEval("source('{}',chdir=TRUE)".format(self.source))
                        fn = getattr(c.r, fname)
                        out = fn(*args, **kwargs)
                        c.close()
                        return out
                    rwrapper = setFunctionName(rwrapper, i)
                    setattr(self, i, rwrapper)
                    del(rwrapper)
            except:
                print("skipping :", i)
                pass
        if verbose:
            print("\r")

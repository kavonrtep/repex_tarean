#!/usr/bin/env python3
import multiprocessing
import os
import time
from itertools import cycle
'''
functions for parallel processing of data chunks using worker function
'''


def run_multiple_pbs_jobs(cmds, status_files, qsub_params=""):
    '''
    Example of pbs_params:
    -l walltime=1000:00:00,nodes=1:ppn=8,mem=15G
    -l walltime=150:00:00,nodes=1:ppn=1

    '''
    jobs = []
    status_function = []
    status_command = []
    for cmd, sf in zip(cmds, status_files):
        jobs.append(pbs_send_job(cmd, sf, qsub_params))
    for p in jobs:
        p.join()
        status_function.append(p.exitcode)
    # collect pbs run status
    for sf in status_files:
        with open(sf) as f:
            status_command.append(f.read().strip())
    status = {'function': status_function, 'command': status_command}
    return status


def pbs_send_job(cmd, status_file, qsub_params):
    ''' send job to pbs cluster, require status file'''
    p = multiprocessing.Process(target=pbs_run,
                                args=(cmd, status_file, qsub_params))
    p.start()
    return p


def pbs_run(cmd, status_file, qsub_params):
    '''
    run shell command cmd on pbs cluster, wait for job to finish
    and return status
    '''
    print(status_file)
    error_file = status_file + ".e"
    # test if writable
    try:
        f = open(status_file, 'w').close()
        f = open(error_file, 'w').close()
    except IOError:
        print("cannot write to status files, make sure path exists")
        raise IOError

    if os.path.exists(status_file):
        print("removing old status file")
        os.remove(status_file)
    cmd_full = ("echo '{cmd} && echo \"OK\" > {status_file} || echo \"ERROR\""
                " > {status_file}' | qsub -e {err}"
                " {qsub_params} ").format(cmd=cmd, status_file=status_file,
                                          err=error_file,
                                          qsub_params=qsub_params)
    os.system(cmd_full)

    while True:
        if os.path.exists(status_file):
            break
        else:
            time.sleep(3)
    with open(status_file) as f:
        status = f.read().strip()
    return status


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def get_max_proc():
    '''Number of cpu to ise in ether get from config.py is available or
    from global PROC or from environment variable PRCO or set to system max'''
    try:
        from config import PROC as max_proc
    except ImportError:
        if "PROC" in globals():
            max_proc = PROC
        elif "PROC" in os.environ:
            max_proc = int(os.environ["PROC"])

        else:
            max_proc = multiprocessing.cpu_count()
    return max_proc


def parmap2(f, X, groups, ppn):
    max_proc = get_max_proc()
    print("running in parallel using ", max_proc, "cpu(s)")
    process_pool = []
    output = [None] * len(X)
    # prepare processes
    for x, index in zip(X, list(range(len(X)))):
        # status:
        # 0: waiting, 1: running, 2:collected
        process_pool.append({
            'status': 0,
            'proc': None,
            'pipe': None,
            'index': index,
            'group': groups[index],
            'ppn': ppn[index]

        })

    # run processes
    running = 0
    finished = 0
    sleep_time = 0.001
    while True:
        # count alive processes
        if not sleep_time:
            sleep_time = 0.001
        for i in process_pool:
            if i['status'] == 1 and not (i['proc'].exitcode is None):
                sleep_time = 0.0
                # was running now finished --> collect
                i['status'] = 2
                running -= 1
                finished += 1
                output[i['index']] = collect(i['proc'], i['pipe'])
                del i['pipe']
                del i['proc']
            if i['status'] == 0 and running < max_proc:
                # waiting and free --> run
                # check if this group can be run
                running_groups = [pp['group']
                                  for pp in process_pool if pp['status'] == 1]
                # check max load  of concurent runs:
                current_load = sum([pp['ppn']
                                    for pp in process_pool if pp['status'] == 1])
                cond1 = (i['ppn'] + current_load) <= 1
                cond2 = not i['group'] in running_groups
                if cond1 and cond2:
                    sleep_time = 0.0
                    try:
                        i['pipe'] = multiprocessing.Pipe()
                    except OSError as e:
                        print('exception occured:',e)
                        continue
                    i['proc'] = multiprocessing.Process(
                        target=spawn(f),
                        args=(i['pipe'][1], X[i['index']]),
                        name=str(i['index']))
                    i['proc'].start()
                    i['status'] = 1
                    running += 1
        if finished == len(process_pool):
            break
        if sleep_time:
            # sleep only if nothing changed in the last cycle
            time.sleep(sleep_time)
            # sleep time gradually increase to 1 sec
            sleep_time = min(2 * sleep_time, 1)  
    return output


def print_status(pp):
    states = ['waiting', 'running', 'collected']
    print("___________________________________")
    print("jobid    status   group   ppn   exitcode")
    print("===================================")
    for i in pp:
        print(
            i['index'], "    ",
            states[i['status']], "    ",
            i['group'], "    ",
            i['ppn'], "    ",
            i['proc'].exitcode
        )


def collect(pf, pp):
    if pf.pid and not pf.exitcode and not pf.is_alive():
        returnvalue = pp[0].recv()
        pf.join()
        pp[0].close()
        pp[1].close()
        return returnvalue
    elif pf.exitcode:
        print("job finished with exit code {}".format(pf.exitcode))
        pf.join()
        pp[0].close()
        pp[1].close()
        return None
    # return None
    else:
        raise Exception('not collected')


def parmap(f, X):

    max_proc = get_max_proc()

    pipe = []
    proc = []
    returnvalue = {}

    for x, index in zip(X, list(range(len(X)))):
        pipe.append(multiprocessing.Pipe())
        proc.append(multiprocessing.Process(target=spawn(f),
                                            args=(pipe[-1][1], x), name=str(index)))
        p = proc[-1]
        # count alive processes
        while True:
            running = 0
            for i in proc:
                if i.is_alive():
                    running += 1
          #          print "running:"+str(running)
            if running < max_proc:
                break
            else:
                time.sleep(0.1)
        p.start()
        # print "process started:"+str(p.pid)
        # check for finished

        for pf, pp, index in zip(proc, pipe, range(len(pipe))):
            if pf.pid and not pf.exitcode and not pf.is_alive() and (pf.name not in returnvalue):
                pf.join()
                returnvalue[str(pf.name)] = pp[0].recv()
                pp[0].close()
                pp[1].close()
                # proc must be garbage collected - to free all file connection
                del proc[index]
                del pipe[index]

    # collect the rest:
    [pf.join() for pf in proc]
    for pf, pp in zip(proc, pipe):
        if pf.pid and not pf.exitcode and not pf.is_alive() and (pf.name not in returnvalue):
            returnvalue[str(pf.name)] = pp[0].recv()
            pp[0].close()
            pp[1].close()
    # convert to list in input correct order
    returnvalue = [returnvalue[str(i)] for i in range(len(X))]
    return returnvalue


def parallel2(command, *args, groups=None, ppn=None):
    ''' same as parallel but groups are used to identifie mutually
    exclusive jobs, jobs with the same goup id are never run together
    ppn params is 'load' of the job - sum of loads cannot exceed 1
    '''
    # check args, expand if necessary
    args = list(args)
    N = [len(i) for i in args]  # lengths of lists
    Mx = max(N)
    if len(set(N)) == 1:
        # all good
        pass
    elif set(N) == set([1, Mx]):
        # expand args of length 1
        for i in range(len(args)):
            if len(args[i]) == 1:
                args[i] = args[i] * Mx
    else:
        raise ValueError
    if not groups:
        groups = range(Mx)
    elif len(groups) != Mx:
        print("length of groups must be same as number of job or None")
        raise ValueError

    if not ppn:
        ppn = [0] * Mx
    elif len(ppn) != Mx:
        print("length of ppn must be same as number of job or None")
        raise ValueError
    elif max(ppn) > 1 and min(ppn):
        print("ppn values must be in 0 - 1 range")
        raise ValueError
    # convert argument to suitable format - 'transpose'
    argsTuples = list(zip(*args))
    args = [list(i) for i in argsTuples]

    # multiprocessing.Pool()

    def command_star(args):
        return(command(*args))

    x = parmap2(command_star,  argsTuples, groups, ppn)
    return x


def parallel(command, *args):
    ''' Execute command in parallel using multiprocessing
    command is the function to be executed
    args is list of list of arguments
    execution is :
        command(args[0][0],args[1][0],args[2][0],args[3][0],....)
        command(args[0][1],args[1][1],args[2][1],args[3][1],....)
        command(args[0][2],args[1][2],args[2][2],args[3][2],....)
        ...
    output of command is returned as list
    '''
    # check args, expand if necessary
    args = list(args)
    N = [len(i) for i in args]  # lengths of lists
    Mx = max(N)
    if len(set(N)) == 1:
        # all good
        pass
    elif set(N) == set([1, Mx]):
        # expand args of length 1
        for i in range(len(args)):
            if len(args[i]) == 1:
                args[i] = args[i] * Mx
    else:
        raise ValueError

    # convert argument to suitable format - 'transpose'
    argsTuples = list(zip(*args))
    args = [list(i) for i in argsTuples]

    multiprocessing.Pool()

    def command_star(args):
        return(command(*args))

    x = parmap(command_star, argsTuples)
    return x


def worker(*a):
    x = 0
    y = 0
    for i in a:
        if i == 1.1:
            print("raising exception")
            s = 1 / 0
        y += i
        for j in range(10):
            x += i
            for j in range(100000):
                x = 1.0 / (float(j) + 1.0)
    return(y)

# test
if __name__ == "__main__":
 #   x = parallel2(worker, [1], [2], [3], [4], [1], [1, 2, 3, 7, 10, 1.1, 20, 30, 40, 10, 30, 20, 40, 50, 50], [
 #       3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 5, 6, 4, 3, 2])

    x = parallel2(
        worker, [1], [2], [3], [4], [1],
        [1, 2, 3, 7, 10, 1.2, 20, 30, 40, 10, 30, 20, 40, 50, 50],
        [3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 5, 6, 4, 3, 2],
        groups=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        ppn=[0.6, 0.6, 0.2, 0.6, 0.2, 0.2, 0.4,
             0.1, 0.1, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1]
    )
    print(x)

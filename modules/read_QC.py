#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from multiprocessing import Pool
import time
import os
import sys
import subprocess


# default parameters
threads2use = '6'


def safeMkdir(dirName):
    work_dir = os.getcwd()
    if dirName in os.listdir():
        print("{} has existed.".format(dirName))
        return dirName
    os.mkdir(os.path.join(work_dir, dirName), 0o755)
    print("Create {}".format(dirName))
    return dirName


def getReadsFile():
    return sys.argv[1:]


def runFastqc(threads, fastqFile, outdir):
    cmd = " ".join(["fastqc", "-o", outdir, "-f", "fastq",
                    "-t", str(threads), fastqFile])
    subprocess.run(cmd, check=True, shell=True)


def main():
    lreads_file = getReadsFile()
    dpreQC_report = safeMkdir('preQC_report')
    processes = len(lreads_file)
    p = Pool(processes)
    for reads_file in lreads_file:
        p.apply_async(runFastqc, args=(threads2use, reads_file, dpreQC_report))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('All fastqc done.')


if __name__ == "__main__":
    start_time = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
    print("Read quality control start:\n{}".format(start_time))
    main()
    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
    print("Read quality control end:\n{}".format(end_time))

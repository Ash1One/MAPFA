#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
import sys
import subprocess

def runFastqc(threads, fastqFile, outdir):
    cmd = " ".join(["fastqc", "-o", outdir, "-f", "fastq",
                    "-t", str(threads), fastqFile])
    subprocess.run(cmd, check=True, shell=True) 

def runFastp(fq_R1, fq_R2, outdir):
    fq_qc_R1 = fq_R1.replace('.fq', '.qc.fq')
    fq_qc_R2 = fq_R2.replace('.fq', '.qc.fq')
    cmd = ' '.join(['fastp', '-i', fq_R1, '-o', fq_qc_R1, '-I', fq_R2, '-O', fq_qc_R2, '--cut_tail', '--length_required=50', '--correction'])
    subprocess.run(cmd, shell=True, check=True)


def rmHostGenome(ref, fq_R1, fq_R2, threads):
    fq = fq_R1.split('_')[:-1][0]
    fq_qc_R1 = fq_R1.replace('.fq', '.qc.fq')
    fq_qc_R2 = fq_R2.replace('.fq', '.qc.fq')
    #cmd_step0 = ' '.join(['bowtie2-build', ref, ref])
    cmd_step1 = ' '.join(['bowtie2', '-x', ref, '-1', fq_qc_R1,
                          '-2', fq_qc_R2, '-S', fq+'.sam', '-p', threads])
    cmd_step2 = ' '.join(['samtools', 'view', '-bS', fq +
                          '.sam', '>', fq+'.bam', '-@', threads])
    cmd_step3 = ' '.join(['samtools', 'view', '-b', '-f', '12', '-F',
                          '256', fq+'.bam', '>', fq+'.unmapped.bam', '-@', threads])
    cmd_step4 = ' '.join(['samtools', 'sort', '-n', '-o', fq +
                          '.unmapped.sorted.bam', fq+'.unmapped.bam', '-@', threads])
    cmd_step5 = ' '.join(['bamToFastq', '-i', fq+'.unmapped.sorted.bam',
                          '-fq', fq+'_R1.hostRemoved.fq', '-fq2', fq+'_R2.hostRemoved.fq'])
    for cmd in [cmd_step1, cmd_step2, cmd_step3, cmd_step4, cmd_step5]:
        subprocess.run(cmd, shell=True, check=True)

'''
if __name__ == "__main__":
    start_time = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
    print("Read quality control start:\n{}".format(start_time))
    main()
    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
    print("Read quality control end:\n{}".format(end_time))
'''
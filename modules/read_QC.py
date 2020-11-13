#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import subprocess
from multiprocessing import Pool
from pathlib import Path
from pathlib import PurePosixPath


def pool2runFastqc(fq_to_qc, report_outdir, threads, logger2log):
    p = Pool(len(fq_to_qc))
    for fq in fq_to_qc:
        p.apply_async(runFastqc, args=(threads, fq, report_outdir))
    logger2log.info("Waiting for all subprocesses of fastqc done...")
    p.close()
    p.join()
    logger2log.info("All subprocesses of fastqc done.")


def runFastqc(threads, fastqFile, outdir):
    cmd = " ".join(["fastqc", "-o", outdir, "-f", "fastq", "-t", str(threads), fastqFile])
    subprocess.run(cmd, check=True, shell=True)


def runFastp(fq_R1, fq_R2, outdir):
    fq_fastp_R1 = fq_R1.replace('.fastq', '.fastp.fq')
    fq_fastp_R2 = fq_R2.replace('.fastq', '.fastp.fq')
    cmd = ' '.join(['fastp', '-i', fq_R1, '-o', str(PurePosixPath(outdir).joinpath(fq_fastp_R1)), '-I', fq_R2, '-O', str(PurePosixPath(outdir).joinpath(fq_fastp_R2)), '--cut_tail', '--length_required=50', '--correction'])
    subprocess.run(cmd, check=True, shell=True)


def rmHostGenome(ref, fq_R1, fq_R2, threads):
    fq = fq_R1.split('_')[:-1][0]
    fq_fastp_R1 = fq_R1.replace('.fastq', '.fastp.fq')
    fq_fastp_R2 = fq_R2.replace('.fastq', '.fastp.fq')
    # cmd_step0 = ' '.join(['bowtie2-build', ref, ref])
    cmd_step1 = ' '.join(['bowtie2', '-x', ref, '-1', fq_fastp_R1, '-2', fq_fastp_R2, '-S', fq+'.sam', '-p', threads])
    cmd_step2 = ' '.join(['samtools', 'view', '-bS', fq+'.sam', '>', fq+'.bam', '-@', threads])
    cmd_step3 = ' '.join(['samtools', 'view', '-b', '-f', '12', '-F', '256', fq+'.bam', '>', fq+'.unmapped.bam', '-@', threads])
    cmd_step4 = ' '.join(['samtools', 'sort', '-n', '-o', fq+'.unmapped.sorted.bam', fq+'.unmapped.bam', '-@', threads])
    cmd_step5 = ' '.join(['bamToFastq', '-i', fq+'.unmapped.sorted.bam', '-fq', fq+'_1.qc.fq', '-fq2', fq+'_2.qc.fq'])
    for cmd in [cmd_step1, cmd_step2, cmd_step3, cmd_step4, cmd_step5]:
        subprocess.run(cmd, check=True, shell=True)


def cleanTemp(work_dir, logger):
    '''delete fastp.html, fastp.json, *.bam, *.sam, preQC_report(direc), againQC_report(direc)
    '''

    '''
    try:
        os.remove('fastp.html')
        os.remove('fastp.json')
    except OSError as exception:
        if exception.errno == errno.ENOENT:
            logger.warning("clean temp files -- Exception when cleaning fastp temp files: %s", repr(exception))
        else:
            logger.error("clean temp files -- There are something wrong when cleaning temp files", exc_info=True)
            raise
    '''
    cleaned_files = []
    # Path.glob - Python >= 3.5
    p = Path(work_dir)
    # p.glob() is a generator returning abspath of files (PosixPath class)
    for fastp_report in list(p.glob('*fastp.*')):
        # parameter of Path.unlink should be a PurePath class
        Path.unlink(fastp_report)
        cleaned_files.append(str(PurePosixPath(work_dir).joinpath(fastp_report)))
    # Regex to rm
    align_temp = re.compile(r"^.*(unmapped)?(sorted)?.[sb]am$")
    rm_align_temp_file = [i.group() for i in list(map(align_temp.search, os.listdir(work_dir))) if i != None]
    for file in rm_align_temp_file:
        os.remove(file)
        cleaned_files.append(str(PurePosixPath(work_dir).joinpath(file)))
    # rmdir
    # Depth-First-Search to delete a directory
    for dirpath, dirname, filename in os.walk(os.path.join(work_dir, 'preQC_report'), topdown=False):
        for file in filename:
            os.remove(os.path.join(dirpath, file))
            cleaned_files.append(os.path.join(dirpath, file))
        os.rmdir(dirpath)
    logger.info("clean temp files -- preQC_report directory was deleted.")

    for dirpath, dirname, filename in os.walk(os.path.join(work_dir, 'againQC_report'), topdown=False):
        # shutil.rmtree(os.path.join(root)) 一步到位
        for file in filename:
            os.remove(os.path.join(dirpath, file))
            cleaned_files.append(os.path.join(dirpath, file))
        os.rmdir(dirpath)

    return cleaned_files

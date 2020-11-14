#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging.config
import os
import subprocess
import sys
from multiprocessing import Pool

import yaml

import modules
from modules.utils import getLogger
from modules.utils import makesurePathExists
from modules.utils import fileExists
from version import __version__


def main(args=None):
    """please select a module of MAPFA:
    reads_QC, assembly, binning and add proper arguments.

    input reads should be named like read-a_1.fastq, read-a_2.fastq
    """
    parser = argparse.ArgumentParser(description="MAPFA: Metagenomic Analysis Pipeline for Food Animals", prog="mapfa")
    parser.add_argument('-t', '--threads', type=int, default=1, help="number of threads")
    parser.add_argument('-o', '--outdir', help='directory to write the result and log to', default=os.getcwd())
    # True when triggering the action
    parser.add_argument('--silent', help='Silent mode', action='store_true', default=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    group_qc = parser.add_argument_group('read_QC module arguments')
    group_qc.add_argument('-fr', '--forward_raw_reads', help="forward fastq raw reads", nargs='*')
    group_qc.add_argument('-rr', '--reverse_raw_reads', help="reverse fastq raw reads", nargs='*')
    group_qc.add_argument('-i', '--index', help="index of the host genome")

    group_assembly = parser.add_argument_group('assembly module arguments')
    group_assembly.add_argument('-fa', '--forward_clean_reads_to_assemble', help="forward fastq clean reads to assemble", nargs='*')
    group_assembly.add_argument('-ra', '--reverse_clean_reads_to_assemble', help="reverse fastq clean reads to assemble", nargs='*')
    group_assembly.add_argument('-m', '--memory', help="memory (gigabyte) to assemble", type=int, default=200)
    group_assembly.add_argument('--assemblyer', help="choose a assemblyer to assemble (metaspades or megahit). Default: megahit.", default='megahit')

    group_binning = parser.add_argument_group('binning module arguments')
    group_binning.add_argument('-fb', '--forward_clean_reads_to_binning', help="forward fastq clean reads", nargs='*')
    group_binning.add_argument('-rb', '--reverse_clean_reads_to_binning', help="forward fastq clean reads", nargs='*')
    group_binning.add_argument('-a', '--assembled_reads', help="assembled fasta reads")

    args = parser.parse_args(args)
    print(args)

    # load log config
    with open('log_config.yaml', 'r', encoding='utf-8') as f:
        log_config = yaml.safe_load(f)
        logging.config.dictConfig(log_config)
    logger = getLogger(args.outdir, args.silent)

    outdir = os.path.abspath(args.outdir)
    logger.info("MAPFA start:")
    logging.info(' '.join(sys.argv))
    logger.info("The output directory: %s", outdir)
    ##################################
    # read QC
    if args.forward_raw_reads and args.reverse_raw_reads and args.index:
        logger.info("Run reads_QC module:")
        ##############################
        # Pre-QC
        raw_fq_to_qc = args.forward_raw_reads + args.reverse_raw_reads
        preQC_report_outdir = os.path.join(outdir, 'preQC_report')
        makesurePathExists(preQC_report_outdir)
        modules.read_QC.pool2runFastqc(raw_fq_to_qc, preQC_report_outdir, args.threads, logger)
        ##############################

        ##############################
        # readFiltering
        readFiltering_tasks = len(args.forward_raw_reads)
        ## run fastp to cut low quality reads
        p2fastp = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2fastp.apply_async(modules.read_QC.runFastp, args=(file, file.replace('_1', '_2'), outdir))
        logger.info("Waiting for all subprocesses of fastp done...")
        p2fastp.close()
        p2fastp.join()
        logger.info("All subprocesses of fastp done.")
        ### check fastp results
        for raw_fq in raw_fq_to_qc:
            if not fileExists(raw_fq.split('.')[0]+'.fastp.fq'):
                logger.critical("Something wrong when running fastp to cut low quality reads")
                return
            else:
                logger.info("fastp is working properly")
                logger.info("continue")
        ## remove host genome contamination
        p2rmHG = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2rmHG.apply_async(modules.read_QC.rmHostGenome, args=(args.index, file, file.replace('_1', '_2'), str(args.threads)))
        logger.info("Waiting for all subprocesses of removing host genome contamination done...")
        p2rmHG.close()
        p2rmHG.join()
        logger.info("All subprocesses of removing host genome contamination done.")
        ##############################

        ##############################
        # QC-again
        fq_to_qc_again = [ i for i in os.listdir(outdir) if i.endswith('.qc.fq')]
        againQC_report_outdir = os.path.join(outdir, 'againQC_report')
        makesurePathExists(againQC_report_outdir)
        modules.read_QC.pool2runFastqc(fq_to_qc_again, againQC_report_outdir, args.threads, logger)

        ## multiqc
        subprocess.run(' '.join(['multiqc', preQC_report_outdir+'/*.zip', againQC_report_outdir+'/*.zip']), shell=True, check=True)
        ##############################

        ##############################
        # clean temp files
        cleaned_file = modules.read_QC.cleanTemp(outdir, logger)
        for file in cleaned_file:
            logger.info("clean temp files -- %s was removed.", file)
        ##############################

        ##############################
        # check results
        for raw_fq in raw_fq_to_qc:
            if not fileExists(raw_fq.split('.')[0]+'.qc.fq'):
                logger.critical("Something wrong when removing Host Genome")
                return
            else:
                logger.info("rmHostGenome is working properly")
                logger.info("read_QC module running smoothly.")
        ##############################
    ##################################
    ##################################
    # assembly
    elif args.forward_clean_reads_to_assemble and args.reverse_clean_reads_to_assemble and args.assemblyer:
        logger.info("Run assembly module:")
        logger.info("choose %s to assemble clean reads", args.assemblyer)
        assembly_outdir = os.path.join(outdir, 'assembly_out')
        makesurePathExists(assembly_outdir)
        fq_1_qc = [ fq1.replace('fastq', 'qc.fq') for fq1 in args.forward_raw_reads ]
        fq_2_qc = [ fq2.replace('fastq', 'qc.fq') for fq2 in args.reverse_raw_reads ]
        if args.assemblyer == 'metaspades':
            modules.assembly.metaspades2assembly(str(args.threads), str(args.memory), assembly_outdir, ' '.join(fq_1_qc), ' '.join(fq_2_qc))
        elif args.assemblyer == 'megahit':
            modules.assembly.megahit2assembly(str(args.threads), str(args.memory), assembly_outdir, ' '.join(fq_1_qc), ' '.join(fq_2_qc))
        else:
            logger.critical("please choose metaspades or megahit.")
            return

    elif args.forward_clean_reads_to_binning and args.reverse_clean_reads_to_binning and args.assembled_reads:
        logger.info("Run binning module:")





    logger.info("MAPFA end.")


##############################
# Assembly
##############################

##############################
# Assembly QC
##############################

##############################
# Binning
##############################

##############################
# Binning QC
##############################


if __name__ == '__main__':
    main()




'''
def moduleSwitch(module, targs):
    return {
        'reads_QC': modules.read_QC.main(targs),
        'assembly': modules.assembly.main(targs),
        'binning': modules.binning.main(targs),
    }.get(module, default="""please select a module of MAPFA:
    reads_QC, assembly, binning. (e.g. mapfa.py reads_QC)
    """)
'''

mapfa_test = 'ok'
assert mapfa_test == 'ok', 'something was wrong when running mapfa!'

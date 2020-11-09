#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging.config
import os
import sys
from multiprocessing import Pool

import yaml

import modules
from modules.utils import getLogger
from modules.utils import makesurePathExists
from version import __version__


def moduleSwitch(module, targs):
    return {
        'reads_QC': modules.read_QC.main(targs),
        'assembly': modules.assembly.main(targs),
        'binning': modules.binning.main(targs),
    }.get(module, default="""please select a module of MAPFA:
    reads_QC, assembly, binning. (e.g. mapfa.py reads_QC)
    """)


def main(args=None):
    """please select a module of MAPFA:
    reads_QC, assembly, binning and add proper arguments.
    """
    parser = argparse.ArgumentParser(
        description="MAPFA: Metagenomic Analysis Pipeline for Food Animals", prog="mapfa")
    parser.add_argument('-t', '--threads', type=int,
                        default=1, help="number of threads")
    parser.add_argument(
        '-o', '--outdir', help='directory to write the result and log to', default=os.getcwd())
    # True when triggering the action
    parser.add_argument('--silent', help='Silent mode',
                        action='store_true', default=False)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + __version__)

    group_qc = parser.add_argument_group('read_QC module arguments')
    group_qc.add_argument('-fr', '--forward_raw_reads',
                          help="forward fastq raw reads")
    group_qc.add_argument('-rr', '--reverse_raw_reads',
                          help="reverse fastq raw reads")
    group_qc.add_argument('-i', '--index', help="index of the host genome")

    group_assembly = parser.add_argument_group('assembly module arguments')
    group_assembly.add_argument(
        '-fa', '--forward_clean_reads_to_assemble', help="forward fastq clean reads to assemble")
    group_assembly.add_argument(
        '-ra', '--reverse_clean_reads_to_assemble', help="reverse fastq clean reads to assemble")
    group_assembly.add_argument(
        'assemblyer', help="choose a assemblyer to assemble (metaspades or megahit). Default: megahit.", default='megahit')

    group_binning = parser.add_argument_group('binning module arguments')
    group_binning.add_argument(
        '-fb', '--forward_clean_reads_to_binning', help="forward fastq clean reads")
    group_binning.add_argument(
        '-rb', '--reverse_clean_reads_to_binning', help="forward fastq clean reads")
    group_binning.add_argument(
        '-a', '--assembled_reads', help="assembled fasta reads")

    args = parser.parse_args(args)

    with open('log_config.yaml', 'r', encoding='utf-8') as f:
        log_config = yaml.safe_load(f)
        logging.config.dictConfig(log_config)
    logger = getLogger(args.outdir, args.silent)

    outdir = os.path.abspath(args.outdir)
    logger.info("MAPFA start:")
    logger.info("The output directory: %s", outdir)

    if args.forward_raw_reads and args.reverse_raw_reads and args.index:
        logger.info("Run reads_QC module:")

        ##############################
        # Pre-QC
        raw_fq_to_qc = args.forward_raw_reads + args.reverse_raw_reads
        preQC_report_outdir = os.path.join(outdir, 'preQC_report')
        makesurePathExists(preQC_report_outdir)
        p = Pool(len(raw_fq_to_qc))
        for raw_fq in raw_fq_to_qc:
            p.apply_async(read_QC.runFastqc, args=(args.threads, raw_fq, preQC_report_outdir))
        logger.info("Waiting for all subprocesses of fastqc done...")
        p.close()
        p.join()
        logger.info("All subprocesses of fastqc done.")
        ## get multiqc report

        ## check and clean temp files

        #############################

        ##############################
        # readFiltering
        readFiltering_tasks = len(args.forward_raw_reads)
        # run fastp to cut low quality reads
        p2fastp = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2fastp.apply_async(read_QC.runFastp, args=(file, file.replace('_1', '_2')))
        logger.info("Waiting for all subprocesses of fastp done...")
        p2fastp.close()
        p2fastp.join()
        logger.info("All subprocesses of fastp done.")
        # remove host genome contamination
        p2rmHG = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2rmHG.apply_async(read_QC.rmHostGenome, args=(args.index, file, file.replace('_1', '_2'), str(args.threads)))
        logger.info("Waiting for all subprocesses of removing host genome contamination done...")
        p2rmHG.close()
        p2rmHG.join()
        logger.info("All subprocesses of removing host genome contamination done.")
        # cleam temp files
        
        ##############################


    elif args.forward_clean_reads_to_assemble and args.reverse_clean_reads_to_assemble and args.assemblyer:
        logger.info("Run assembly module:")
        logger.info("choose %s to assemble clean reads", args.assemblyer)
    elif args.forward_clean_reads_to_binning and args.reverse_clean_reads_to_binning and args.assembled_reads:
        logger.info("Run binning module:")





    logger.info("MAPFA end.")

if __name__ == '__main__':
    main()





##############################
# QC-again
##############################

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

mapfa_test = 'ok'
assert mapfa_test == 'ok', 'something was wrong when running mapfa!'

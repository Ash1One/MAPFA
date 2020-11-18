#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging.config
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool

import yaml

from modules.read_QC import pool2runFastqc
from modules.read_QC import runFastp
from modules.read_QC import rmHostGenome
from modules.read_QC import cleanTemp
from modules.assembly import rmMetaspadesShortContigs
from modules.assembly import fixMegahitContigName
from modules.assembly import metaspades2assembly
from modules.assembly import megahit2assembly
from modules.assembly import quast2QC
from modules.utils import getLogger
from modules.utils import makesurePathExists
from modules.utils import fileExists
from version import __version__


def main(args=None):
    """please select a module of MAPFA:
    reads_QC, assembly, binning and add proper arguments.

    input reads should be named like read-a_1.fastq, read-a_2.fastq
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
    # QC group
    group_qc = parser.add_argument_group('read_QC module arguments')
    group_qc.add_argument('-fr', '--forward_raw_reads',
                          help="forward fastq raw reads", nargs='*')
    group_qc.add_argument('-rr', '--reverse_raw_reads',
                          help="reverse fastq raw reads", nargs='*')
    group_qc.add_argument('-i', '--index', help="index of the host genome")
    # assembly group
    group_assembly = parser.add_argument_group('assembly module arguments')
    group_assembly.add_argument('-fa', '--forward_clean_reads_to_assemble',
                                help="forward fastq clean reads to assemble", nargs='*')
    group_assembly.add_argument('-ra', '--reverse_clean_reads_to_assemble',
                                help="reverse fastq clean reads to assemble", nargs='*')
    group_assembly.add_argument(
        '-ma', '--memory4assemble', help="memory (gigabyte) for assembling. Default: 30 (GB)", type=int, default=30)
    group_assembly.add_argument('-la', '--minlength4a', help="keep long contigs over minlength. Default: 500 (bp).", default=500)
    group_assembly.add_argument(
        '--assemblyer', help="choose a assemblyer to assemble (metaspades or megahit). Default: megahit.", default='megahit')
    # binning group
    group_binning = parser.add_argument_group('binning module arguments')
    group_binning.add_argument(
        '-fb', '--forward_clean_reads_to_binning', help="forward fastq clean reads", nargs='*')
    group_binning.add_argument(
        '-rb', '--reverse_clean_reads_to_binning', help="forward fastq clean reads", nargs='*')
    group_binning.add_argument('-mb', '--memory4binning', help="memory (gigabyte) for binning. Default: 30 (GB).", default=30)
    group_binning.add_argument('-lb', '--minlength4b', help="contigs with minimum length to bin. Default: 1000 (bp).", default=1000)
    group_binning.add_argument(
        '-a', '--assembled_reads', help="assembled fasta reads")
    group_binning.add_argument('--metabat2', help="metabat2 to bin contigs", action="store_true", default=False)
    group_binning.add_argument('--maxbin2', help="maxbin2 to bin contigs", action="store_true", default=False)
    group_binning.add_argument('--groopm2', help="groopm2 to bin contigs", action="store_true", default=False)
    group_binning.add_argument('--concoct', help="concoct to bin contigs", action="store_true", default=False)
    #group_binning.add_argument('--vamb', help="vamb to bin contigs", action="store_true", default=False)

    if len(sys.argv) == 1:
        parser.print_help()
        return
    args = parser.parse_args(args)
    print(args)

    # load log config
    logger = getLogger(args.outdir, args.silent)

    outdir = os.path.abspath(args.outdir)
    logger.info("MAPFA start:")
    logging.info(' '.join(sys.argv))
    logger.info("The output directory: %s", outdir)
    ####################################################################
    # read QC
    if args.forward_raw_reads and args.reverse_raw_reads and args.index:
        logger.info("Run reads_QC module:")
        ##############################
        # Pre-QC
        raw_fq_to_qc = args.forward_raw_reads + args.reverse_raw_reads
        preQC_report_outdir = os.path.join(outdir, 'preQC_report')
        makesurePathExists(preQC_report_outdir)
        pool2runFastqc(
            raw_fq_to_qc, preQC_report_outdir, args.threads)
        ##############################

        ##############################
        # readFiltering
        readFiltering_tasks = len(args.forward_raw_reads)
        # run fastp to cut low quality reads
        p2fastp = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2fastp.apply_async(runFastp, args=(
                file, file.replace('_1', '_2'), outdir))
        logger.info("Waiting for all subprocesses of fastp done...")
        p2fastp.close()
        p2fastp.join()
        logger.info("All subprocesses of fastp done.")
        # check fastp results
        for raw_fq in raw_fq_to_qc:
            if not fileExists(raw_fq.split('.')[0]+'.fastp.fq'):
                logger.critical(
                    "Something wrong when running fastp to cut low quality reads")
                return
            else:
                logger.info("fastp is working properly")
                logger.info("continue")
        # remove host genome contamination
        p2rmHG = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2rmHG.apply_async(rmHostGenome, args=(
                args.index, file, file.replace('_1', '_2'), str(args.threads)))
        logger.info(
            "Waiting for all subprocesses of removing host genome contamination done...")
        p2rmHG.close()
        p2rmHG.join()
        logger.info(
            "All subprocesses of removing host genome contamination done.")
        ##############################

        ##############################
        # QC-again
        fq_to_qc_again = [i for i in os.listdir(
            outdir) if i.endswith('.qc.fq')]
        againQC_report_outdir = os.path.join(outdir, 'againQC_report')
        makesurePathExists(againQC_report_outdir)
        pool2runFastqc(
            fq_to_qc_again, againQC_report_outdir, args.threads)

        # multiqc
        subprocess.run(' '.join(['multiqc', preQC_report_outdir+'/*.zip',
                                 againQC_report_outdir+'/*.zip']), shell=True, check=True)
        ##############################

        ##############################
        # clean temp files
        cleaned_file = cleanTemp(outdir)
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
    ####################################################################
    ####################################################################
    # assembly
    elif args.forward_clean_reads_to_assemble and args.reverse_clean_reads_to_assemble and args.assemblyer:
        logger.info("Run assembly module:")
        logger.info("choose %s to assemble clean reads", args.assemblyer)
        # fq_1_qc = [ fq1.replace('fastq', 'qc.fq') for fq1 in args.forward_raw_reads ]
        # fq_2_qc = [ fq2.replace('fastq', 'qc.fq') for fq2 in args.reverse_raw_reads ]
        fq_1_qc = args.forward_clean_reads_to_assemble
        fq_2_qc = args.reverse_clean_reads_to_assemble
        assembled_path = ''

        # Assembling with metaspades
        if args.assemblyer == 'metaspades':
            metaspades_outdir = os.path.join(outdir, 'metaspades_out')
            metaspades_temp = os.path.join(outdir, 'metaspades_temp')
            makesurePathExists(metaspades_outdir)
            makesurePathExists(metaspades_temp)

            logger.info("Assembling with metaspades:")
            logger.info("%s and %s are used for assembling.",
                        ','.join(fq_1_qc), ','.join(fq_2_qc))

            if fileExists(os.path.join(metaspades_outdir, 'spades.log')):
                logger.info("metaspades restart from last running.")
                metaspades2assembly(1, str(args.threads), str(
                    args.memory4assemble), metaspades_outdir, metaspades_temp, ' '.join(fq_1_qc), ' '.join(fq_2_qc))
            else:
                metaspades2assembly(0, str(args.threads), str(
                    args.memory4assemble), metaspades_outdir, metaspades_temp, ' '.join(fq_1_qc), ' '.join(fq_2_qc))
            ## check the result        
            if not fileExists(os.path.join(metaspades_outdir, 'scaffolds.fasta')):
                logger.critical("Something wrong when assembling with metaspades!")
                return
            else:
                assembled_path = os.path.join(metaspades_outdir, 'scaffolds.fasta')
                try:
                    shutil.rmtree(metaspades_temp)
                except Exception:
                    logger.error("There are something wrong when remove metaspades tempdirectory", exc_info=True)
                    raise
                else:
                    logger.info("Clean reads are assembled with metaspades.")
                
        # Assembling with megahit
        elif args.assemblyer == 'megahit':
            megahit_outdir = os.path.join(outdir, 'megahit_out')
            megahit_temp = os.path.join(outdir, 'megahit_temp')
            makesurePathExists(megahit_outdir)
            makesurePathExists(megahit_temp)
            logger.info("Assembling with megahit:")
            megahit2assembly(str(args.threads), str(
                args.memory4assemble), megahit_outdir, megahit_temp, ' '.join(fq_1_qc), ' '.join(fq_2_qc))

            ## check the result
            if not fileExists(os.path.join(megahit_outdir, 'final.contigs.fa')):
                logger.critical("Something wrong when assembling with megahit!")
                return
            else:
                assembled_path = os.path.join(megahit_outdir, 'final.contigs.fa')
                try:
                    shutil.rmtree(megahit_temp)
                except Exception:
                    logger.error("There are something wrong when remove megahit tempdirectory", exc_info=True)
                    raise
                else:
                    logger.info("Clean reads are assembled with megahit.")

        else:
            logger.critical("please choose metaspades or megahit for assembling.")
            return

        # FORMAT the result
        if args.assemblyer == 'metaspades':
            rmMetaspadesShortContigs(args.minlength4a, assembled_path)
        elif args.assemblyer == 'megahit':
            fixMegahitContigName(args.minlength4a, assembled_path)
        
        # check the result
        if not fileExists(os.path.join(outdir, 'assembly.fasta')):
            logger.critical("Something wrong when removing short contigs")
            return
        else:
            logger.info("Assembly result is located in %s", os.path.join(outdir, 'assembly,fasta'))
        
        # assembly QC with QUAST
        logger.info("Running assembly QC with quast:")
        quast_outdir = os.path.join(outdir, 'quast_out')
        makesurePathExists(quast_outdir)
        quast2QC(args.threads, quast_outdir, os.path.join(outdir, 'assembly.fasta'))
        if not fileExists(os.path.join(quast_outdir, 'report.html')):
            logger.critical("Something wrong when assembly QC with quast!")
            return
        else:
            shutil.copy(os.path.join(quast_outdir, 'report.html'), os.path.join(outdir, 'assembly_QCreport.html'))
            logger.info("Assembly QC report is located in %s", os.path.join(outdir, 'assembly_QCreport.html'))

        logger.info("assembly module running smoothly.")

    elif args.forward_clean_reads_to_binning and args.reverse_clean_reads_to_binning and args.assembled_reads:
        if not (args.metabat2 or args.maxbin2 or args.groopm2 or args.concoct):
            logger.error("please select at least one binning tool. e.g. --metabat2")
            return
        

        logger.info("Run binning module:")

    else:
        logger.warning("Please choose a module for MAPFA.")
        return

    logger.info("MAPFA end.")

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
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging.config
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import yaml
import modules
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
    group_assembly.add_argument(
        '-la', '--minlength4a', help="keep long contigs over minlength. Default: 500 (bp).", type=int, default=500)
    group_assembly.add_argument(
        '--assemblyer', help="choose a assemblyer to assemble (metaspades or megahit). Default: megahit.", default='megahit')
    # binning group
    group_binning = parser.add_argument_group('binning module arguments')
    group_binning.add_argument(
        '-fb', '--forward_clean_reads_to_binning', help="forward fastq clean reads", nargs='*')
    group_binning.add_argument(
        '-rb', '--reverse_clean_reads_to_binning', help="forward fastq clean reads", nargs='*')
    group_binning.add_argument(
        '-mb', '--memory4binning', help="memory (gigabyte) for binning. Default: 30 (GB).", type=int, default=30)
    group_binning.add_argument(
        '-lb', '--minlength4b', help="contigs with minimum length to bin. Default: 1000 (bp).", type=int, default=1000)
    group_binning.add_argument(
        '-a', '--assembled_contigs', help="assembled fasta contigs")
    group_binning.add_argument(
        '--metabat2', help="metabat2 to bin contigs", action="store_true", default=False)
    group_binning.add_argument(
        '--maxbin2', help="maxbin2 to bin contigs", action="store_true", default=False)
    group_binning.add_argument(
        '--groopm2', help="groopm2 to bin contigs", action="store_true", default=False)
    group_binning.add_argument(
        '--concoct', help="concoct to bin contigs", action="store_true", default=False)
    # group_binning.add_argument('--vamb', help="vamb to bin contigs", action="store_true", default=False)

    if len(sys.argv) == 1:
        parser.print_help()
        return
    args = parser.parse_args(args)
    print(args)

    # load log config
    logger = modules.getLogger(args.outdir, args.silent)
    # check raw reads
    if not args.forward_raw_reads or not args.reverse_raw_reads:
        logger.critical("please ensure that input raw reads files exist!")
        parser.error("please input raw reads files!")
    raw_fq_to_qc = args.forward_raw_reads + args.reverse_raw_reads
    if False in list(map(modules.fileExists, raw_fq_to_qc)):
        logger.critical("please ensure that input raw reads files exist!")
        parser.error("please ensure that input raw reads files exist!")

    outdir = os.path.abspath(args.outdir)
    logger.info("MAPFA start:")
    logging.info(' '.join(sys.argv))
    logger.info("The output directory: %s", outdir)
    ####################################################################
    ####################################################################
    # read QC
    if args.forward_raw_reads and args.reverse_raw_reads and args.index:
        logger.info("Run reads_QC module:")
        ##############################
        # Pre-QC
        preQC_report_outdir = os.path.join(outdir, 'preQC_report')
        modules.makesurePathExists(preQC_report_outdir)
        modules.pool2runFastqc(
            raw_fq_to_qc, preQC_report_outdir, args.threads)
        ##############################

        ##############################
        # readFiltering
        readFiltering_tasks = len(args.forward_raw_reads)
        # run fastp to cut low quality reads
        p2fastp = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2fastp.apply_async(modules.runFastp, args=(
                file, file.replace('_1', '_2'), outdir))
        logger.info("Waiting for all subprocesses of fastp done...")
        p2fastp.close()
        p2fastp.join()
        logger.info("All subprocesses of fastp done.")
        # check fastp results
        for raw_fq in raw_fq_to_qc:
            if not modules.fileExists(raw_fq.split('.')[0]+'.fastp.fq'):
                logger.critical(
                    "Something wrong when running fastp to cut low quality reads")
                return
            else:
                logger.info("fastp is working properly")
                logger.info("continue")
        # remove host genome contamination
        p2rmHG = Pool(readFiltering_tasks)
        for file in args.forward_raw_reads:
            p2rmHG.apply_async(modules.rmHostGenome, args=(
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
        modules.makesurePathExists(againQC_report_outdir)
        modules.pool2runFastqc(
            fq_to_qc_again, againQC_report_outdir, args.threads)

        # multiqc
        subprocess.run(' '.join(['multiqc', preQC_report_outdir+'/*.zip',
                                 againQC_report_outdir+'/*.zip']), shell=True, check=True)
        ##############################

        ##############################
        # clean temp files
        cleaned_file = modules.cleanTemp(outdir)
        for file in cleaned_file:
            logger.info("clean temp files -- %s was removed.", file)
        ##############################

        ##############################
        # check results
        for raw_fq in raw_fq_to_qc:
            if not modules.fileExists(raw_fq.split('.')[0]+'.qc.fq'):
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
        fq_1_qc = args.forward_clean_reads_to_assemble
        fq_2_qc = args.reverse_clean_reads_to_assemble
        # check input clean reads files
        if False in list(map(modules.fileExists, fq_1_qc+fq_2_qc)):
            logger.critical("please ensure that input clean reads files exist!")
            parser.error("please ensure that input clean reads files exist!")

        logger.info("Run assembly module:")
        logger.info("choose %s to assemble clean reads", args.assemblyer)
        # fq_1_qc = [ fq1.replace('fastq', 'qc.fq') for fq1 in args.forward_raw_reads ]
        # fq_2_qc = [ fq2.replace('fastq', 'qc.fq') for fq2 in args.reverse_raw_reads ]

        assembled_path = ''

        # Assembling with metaspades
        if args.assemblyer == 'metaspades':
            metaspades_outdir = os.path.join(outdir, 'metaspades_out')
            metaspades_temp = os.path.join(outdir, 'metaspades_temp')
            modules.makesurePathExists(metaspades_outdir)
            modules.makesurePathExists(metaspades_temp)

            logger.info("Assembling with metaspades:")
            logger.info("%s and %s are used for assembling.",
                        ','.join(fq_1_qc), ','.join(fq_2_qc))

            if modules.fileExists(os.path.join(metaspades_outdir, 'spades.log')):
                logger.info("metaspades restart from last running.")
                modules.metaspades2assembly(1, str(args.threads), str(
                    args.memory4assemble), metaspades_outdir, metaspades_temp, ' '.join(fq_1_qc), ' '.join(fq_2_qc))
            else:
                modules.metaspades2assembly(0, str(args.threads), str(
                    args.memory4assemble), metaspades_outdir, metaspades_temp, ' '.join(fq_1_qc), ' '.join(fq_2_qc))
            # check the result
            if not modules.fileExists(os.path.join(metaspades_outdir, 'scaffolds.fasta')):
                logger.critical(
                    "Something wrong when assembling with metaspades!")
                return
            else:
                assembled_path = os.path.join(
                    metaspades_outdir, 'scaffolds.fasta')
                try:
                    shutil.rmtree(metaspades_temp)
                except Exception:
                    logger.error(
                        "There are something wrong when remove metaspades tempdirectory", exc_info=True)
                    raise
                else:
                    logger.info("Clean reads are assembled with metaspades.")

        # Assembling with megahit
        elif args.assemblyer == 'megahit':
            megahit_outdir = os.path.join(outdir, 'megahit_out')
            megahit_temp = os.path.join(outdir, 'megahit_temp')
            modules.makesurePathExists(megahit_outdir)
            modules.makesurePathExists(megahit_temp)
            logger.info("Assembling with megahit:")
            modules.megahit2assembly(str(args.threads), str(
                args.memory4assemble), megahit_outdir, megahit_temp, ' '.join(fq_1_qc), ' '.join(fq_2_qc))

            # check the result
            if not modules.fileExists(os.path.join(megahit_outdir, 'final.contigs.fa')):
                logger.critical(
                    "Something wrong when assembling with megahit!")
                return
            else:
                assembled_path = os.path.join(
                    megahit_outdir, 'final.contigs.fa')
                try:
                    shutil.rmtree(megahit_temp)
                except Exception:
                    logger.error(
                        "There are something wrong when remove megahit tempdirectory", exc_info=True)
                    raise
                else:
                    logger.info("Clean reads are assembled with megahit.")

        else:
            logger.critical(
                "please choose metaspades or megahit for assembling.")
            parser.error("please choose metaspades or megahit for assembling.")

        # FORMAT the result
        if args.assemblyer == 'metaspades':
            modules.rmMetaspadesShortContigs(args.minlength4a, assembled_path)
        elif args.assemblyer == 'megahit':
            modules.fixMegahitContigName(args.minlength4a, assembled_path)

        # check the result
        if not modules.fileExists(os.path.join(outdir, 'assembly.fasta')):
            logger.critical("Something wrong when removing short contigs")
            return
        else:
            logger.info("Assembly result is located in %s",
                        os.path.join(outdir, 'assembly,fasta'))

        # assembly QC with QUAST
        logger.info("Running assembly QC with quast:")
        quast_outdir = os.path.join(outdir, 'quast_out')
        modules.makesurePathExists(quast_outdir)
        modules.quast2QC(args.threads, quast_outdir,
                         os.path.join(outdir, 'assembly.fasta'))
        if not modules.fileExists(os.path.join(quast_outdir, 'report.html')):
            logger.critical("Something wrong when assembly QC with quast!")
            return
        else:
            shutil.copy(os.path.join(quast_outdir, 'report.html'),
                        os.path.join(outdir, 'assembly_QCreport.html'))
            logger.info("Assembly QC report is located in %s",
                        os.path.join(outdir, 'assembly_QCreport.html'))

        logger.info("assembly module running smoothly.")
    ####################################################################
    ####################################################################
    # binning
    elif args.forward_clean_reads_to_binning and args.reverse_clean_reads_to_binning and args.assembled_contigs:
        if not (args.metabat2 or args.maxbin2 or args.groopm2 or args.concoct):
            logger.error(
                "please select at least one binning tool. e.g. --metabat2")
            parser.error("please select at least one binning tool. e.g. --metabat2")

        fq_1_4binning = args.forward_clean_reads_to_binning
        fq_2_4binning = args.reverse_clean_reads_to_binning
        # check the assembled_reads and input clean reads files
        if not modules.fileExists(args.assembled_contigs):
            logger.critical("please ensure that assembled contigs exist!")
        if False in list(map(modules.fileExists, fq_1_4binning+fq_2_4binning)):
            logger.critical("please ensure that input clean reads files exist!")
            parser.error("please ensure that input clean reads files exist!")        
        logger.info("Run binning module:")
        # binning work directory
        binning_outdir = os.path.join(outdir, 'binning')
        modules.makesurePathExists(binning_outdir)
        logger.info("Binning outdir locate in %s", binning_outdir)

        ## remove checkm directory if it existed
        if not list(Path(binning_outdir).glob('*checkm')):
            logger.warning("checkm results has been existed.")
            logger.info("Removing former checkm results.")
            for checkm_dir in list(Path(binning_outdir).glob('*checkm')):
                shutil.rmtree(str(checkm_dir))

        # align clean reads to assembly for making coverage files
        ## make the alignment output directory
        align_outdir = os.path.join(binning_outdir, 'align_outdir')
        modules.makesurePathExists(align_outdir)
        logger.info("align_outdir locate in %s", align_outdir)
        ### copy assembly file
        assembly4bin = os.path.join(align_outdir, 'assembly.fa')
        shutil.copy(args.assembled_contigs,
                        assembly4bin)
        logger.info("Copy assembly file to %s", align_outdir)
        ## index the assembly
        if not modules.fileExists(os.path.join(align_outdir, 'assembly.fa.bwt')):
            logger.info("index the assembly...")
            subprocess.run(' '.join(['bwa', 'index', assembly4bin]), shell=True, check=True)
            if not modules.fileExists(os.path.join(align_outdir, 'assembly.fa.bwt')):
                logger.error("Something wrong when index assembly files. Exit.")
                return
        else:
            logger.info("assembly index already existed. So skipping this step.")

        ## alignment (paired end reads)
        for read in args.forward_clean_reads_to_binning:
            if read.endswith('_1.fastq') or read.endswith('_1.qc.fq') or read.endswith('_1.fq') or read.endswith('_1.qc.fastq'):
                if not modules.fileExists(os.path.join(align_outdir, sample_name+'.bam')):
                    read_1 = read
                    read_2 = read.replace('_1', '_2')
                    sample_name = read.split('_')[0]
                    logger.info("making %s alignment files", sample_name)
                    subprocess.run(' '.join(['bwa', 'mem', '-t', args.threads, '-v', '1', assembly4bin, read_1, read_2, '>', os.path.join(align_outdir, sample_name+'.sam')]), shell=True, check=True)
                    
                    logger.info("sorting %s alignment files", sample_name)
                    subprocess.run(' '.join(['samtools', 'sort', '-@', args.threads, '-O', 'BAM', '-o', os.path.join(align_outdir, sample_name+'.bam'), os.path.join(align_outdir, sample_name+'.sam')]))

                else:
                    logger.info("%s bam files already existed. So skipping alignment.")

            else:
                logger.error("Please input reads files with correct name!")
                parser.error("Please input reads files with correct name!")
        
        try:
            os.remove(os.path.join(align_outdir, sample_name+'.sam'))
        except Exception:
            logger.error("It failed to remove sam files. Something wrong.", exc_info=True)

        tbam_files = tuple([str(bam) for bam in tuple(Path(align_outdir).glob('*.bam'))])
        # index the bam files
        for bam in tbam_files:
            subprocess.run(' '.join(['samtools', 'index', '-@', str(args.threads), bam]))

        # run binning tools
        ## run metabat2
        if args.metabat2:
            logger.info("Running metabat2")
            if args.minlength4b <= 1500:
                minlength4metabat = 1500
            else:
                minlength4metabat = args.minlength4b
            
            metabat2_bin_outdir = modules.metabat2bin(align_outdir, assembly4bin, minlength4metabat, args.threads, tbam_files)
            if metabat2_bin_outdir:
                # checkm
                pass

        ## run maxbin2
        if args.maxbin2:
            maxbin2_bin_outdir = modules.maxbin2bin(align_outdir, assembly4bin, args.minlength4b, args.threads, tbam_files)
            if maxbin2_bin_outdir:
                # checkm
                pass
            else:
                print("please check the log file.")
                return
        ## run groopm2
        if args.groopm2:
            groopm2_bin_outdir = modules.groopm2bin(align_outdir, assembly4bin, args.threads, tbam_files)
            if groopm2_bin_outdir:
                # checkm
                pass
            else:
                print("please check the log file.")
                return
        ## concocct
        if args.concoct:
            concoct_bin_outdir = modules.concoct2bin(align_outdir, assembly4bin, args.minlength4b, args.threads, tbam_files)
            if concoct_bin_outdir:
                # checkm
                pass
            else:
                print("please check the log file")
                return
    else:
        logger.warning("Please choose a module for MAPFA.")
        parser.error("Please choose a module for MAPFA.")

    logger.info("MAPFA end.")

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

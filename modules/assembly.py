import os
import subprocess
import sys
import textwrap

from .utils import fileExists
from .utils import makesurePathExists


def fixMegahitContigName(min_length, assembled_fasta):
    dic = {}
    tmp_contig = ""
    good = True
    workdir = os.path.dirname(os.path.dirname(
        os.path.abspath(assembled_fasta)))
    output = os.path.join(workdir, 'assembly.fasta')
    if fileExists(output):
        os.remove(output)
    makesurePathExists(output)

    with open(assembled_fasta, 'r') as f:
        for line in f.readlines():
            if line[0] == ">":
                if tmp_contig != "":
                    if good == True:
                        dic[name] = tmp_contig
                    tmp_contig = ""
                cut = line.strip()[1:].split(" ")

                if int(cut[3].split("=")[1]) < min_length:
                    good = False
                else:
                    good = True

                name = ">"+cut[0]+"_length_" + \
                    cut[3].split("=")[1]+"_cov_"+cut[2].split("=")[1]
            else:
                tmp_contig += line.strip()
        dic[name] = tmp_contig

    for k in sorted(dic, key=lambda k: len(dic[k]), reverse=True):
        wLine2file(
            k+'\n'+textwrap.fill(dic[k], 100, break_on_hyphens=False)+'\n', output)


def wLine2file(line, file):
    with open(file, 'w+') as f:
        f.write(line)


def rmMetaspadesShortContigs(min_length, assembled_fasta):
    workdir = os.path.dirname(os.path.dirname(
        os.path.abspath(assembled_fasta)))
    output = os.path.join(workdir, 'assembly.fasta')
    if fileExists(output):
        os.remove(output)
    makesurePathExists(output)
    with open(assembled_fasta, 'r') as f:
        for line in f.readlines():
            if not line.startswith(">"):
                wLine2file(line.strip+'\n', output)
            else:
                if int(line.split("_")[3]) < int(min_length):
                    break
                else:
                    wLine2file(line.strip()+'\n', output)


def metaspades2assembly(restart, threads, memory_gb, outdir, tempdir, fq_R1, fq_R2):
    if restart:
        cmd = ' '.join(['metaspades.py', '-t', threads, '-m', memory_gb,
                        '-o', outdir, '--tem-dir', tempdir, '-1', fq_R1, '-2', fq_R2, '--restart-from last'])
    else:
        cmd = ' '.join(['metaspades.py', '-t', threads, '-m', memory_gb,
                        '-o', outdir, '--tem-dir', tempdir, '-1', fq_R1, '-2', fq_R2])
    subprocess.run(cmd, shell=True, check=True)


def megahit2assembly(threads, memory_gb, outdir, tempdir, fq_R1, fq_R2):
    cmd = ' '.join(['megahit', '-1', fq_R1, '-2', fq_R2, '-o', outdir,
                    't', threads, '-m', memory_gb+'000000000', '--continue'])
    subprocess.run(cmd, shell=True, check=True)


def quast2QC(threads, outdir, fa_assembly):
    default_minlength = '500'
    cmd = ' '.join(['quast', '-t', threads, '-o',
                    outdir, '-m', default_minlength, fa_assembly])
    subprocess.run(cmd, shell=True, check=True)

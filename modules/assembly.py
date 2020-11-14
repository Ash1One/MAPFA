import subprocess


def metaspades2assembly(threads, memory_gb, outdir, fq_R1, fq_R2):
    cmd = ' '.join(['metaspades.py', '-t', threads, '-m', memory_gb,
                    '-o', outdir, '-1', fq_R1, '-2', fq_R2])
    subprocess.run(cmd, shell=True, check=True)

def megahit2assembly(threads, memory_gb, outdir, fq_R1, fq_R2):
    cmd = ' '.join(['megahit', '-1', fq_R1, '-2', fq_R2, '-o', outdir, 't', threads, '-m', memory_gb+'000000000', '--continue'])
    subprocess.run(cmd, shell=True, check=True)

def quast2QC(threads, outdir, fa_assembly):
    cmd = ' '.join(['quast', '-t', threads, '-o', outdir, '-m', '500', fa_assembly])
    subprocess.run(cmd, shell=True, check=True)
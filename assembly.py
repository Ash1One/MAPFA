import subprocess


# default
threads2use = '10'
memory_gb = '200'


def metaspades2assembly(threads, memory, outdir, fq_R1, fq_R2):
    cmd = ' '.join(['metaspades.py', '-t', threads2use, '-m', memory_gb,
                    '-o', outdir, '-1', fq_R1, '-2', fq_R2])
    subprocess.run(cmd, shell=True, check=True)


def megahit2assembly(threads, fq_R1, fq_R2, outdir):
    cmd = ' '.join(['megahit', '-1', fq_R1, '-2', fq_R2, '-o', outdir, 't', threads, '-m', memory_gb+'000000000', '--continue'])
    subprocess.run(cmd, shell=True, check=True)

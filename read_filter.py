import sys
import subprocess
from multiprocessing import Pool


threads = '10'
ref = 'ref.fa'


def runFastp(fq_R1, fq_R2):
    fq_qc_R1 = fq_R1.split('.')[0] + '.qc.fq'
    fq_qc_R2 = fq_R2.split('.')[0] + '.qc.fq'
    cmd = ' '.join(['fastp', '-i', fq_R1, '-o', fq_qc_R1, '-I', fq_R2, '-O',
                    fq_qc_R2, '--cut_tail', '--length_required=50', '--correction'])
    subprocess.run(cmd, shell=True, check=True)


def rmHostGenome(ref, fq_R1, fq_R2, thread):
    fq = fq_R1.split('.')[0][:-3]
    fq_qc_R1 = fq_R1.split('.')[0] + '.qc.fq'
    fq_qc_R2 = fq_R2.split('.')[0] + '.qc.fq'
    cmd_step1 = ' '.join(['bowtie2-build', ref, ref])
    cmd_step2 = ' '.join(['bowtie2', '-x', ref, '-1', fq_qc_R1,
                          '-2', fq_qc_R2, '-S', fq+'.sam', '-p', thread])
    cmd_step3 = ' '.join(['samtools', 'view', '-bS', fq +
                          '.sam', '>', fq+'.bam', '-@', thread])
    cmd_step4 = ' '.join(['samtools', 'view', '-b', '-f', '12', '-F',
                          '256', fq+'.bam', '>', fq+'.unmapped.bam', '-@', thread])
    cmd_step5 = ' '.join(['samtools', 'sort', '-n', '-o', fq +
                          '.unmapped.sorted.bam', fq+'.unmapped.bam', '-@', thread])
    cmd_step6 = ' '.join(['bamToFastq', '-i', fq+'.unmapped.sorted.bam',
                          '-fq', fq+'_R1.hostRemoved.fq', '-fq2', fq+'_R2.hostRemoved.fq'])
    for cmd in [cmd_step1, cmd_step2, cmd_step3, cmd_step4, cmd_step5, cmd_step6]:
        subprocess.run(cmd, shell=True, check=True)


def main():
    fq_file = sys.argv[1:]
    if len(fq_file) % 2 != 0 or len(fq_file) == 0:
        raise Exception('Please input paired-end fastq files.')

    parallel_task = int(len(fq_file)/2)

    """
    p2fastp = Pool(parallel_task)
    for i in range(0, parallel_task, 2):
        p2fastp.apply_async(runFastp, args=(fq_file[i], fq_file[i+1]))
    p2fastp.close()
    p2fastp.join()
    """   

    ref_genome = ref
    p2rm = Pool(parallel_task)
    for i in range(0, parallel_task, 2):
        p2rm.apply_async(rmHostGenome, args=(
            ref_genome, fq_file[i], fq_file[i+1], threads))
    p2rm.close()
    p2rm.join()

if __name__ == "__main__":
    main()
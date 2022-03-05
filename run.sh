outdir=$1
shift 1;
mkdir -p ${outdir}

DIR="$( cd "$( dirname "$0"  )" && pwd  )"
cluster="bsub -q high -n {threads} -J meta.{wildcards.sample}.{rule} -R \"span[hosts=1]\" -o LOG/{wildcards.sample}/{rule}.out -e {wildcards.sample}/{rule}.err"
snakemake --nolock -s ${DIR}/mapfa.snakemake --latency-wait 45 -k -j 9999 --cluster "$cluster" -C outdir=${outidr} $@
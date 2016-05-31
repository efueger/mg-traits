#!/bin/bash

frag_gene_scan=/bioinf/software/fraggenescan/fraggenescan-1.19/run_FragGeneScan.pl

NAM=$1
NSLOTS=$2
NSEQ=$3

echo $NAM $NSLOTS $NSEQ

asplit '^>' ${NSEQ} spout_ < ../${NAM}.SR.qc.fasta
NFILES=$(ls spout_* | wc -l)

cat > fgs_runner << EOF
#!/bin/bash
#$ -j y
#$ -t 1-${NFILES}
#$ -pe threaded 4
#$ -cwd

############################
# define functions
############################

function time_start {
  date +%s.%N
}

function time_diff {
  END=$(date +%s.%N); DIFF=$(echo "$END - $1" | bc) 
  echo -e $DIFF 
}

############################
# run fgs
############################

#START=$( time_start )
$frag_gene_scan -genome=spout_\${SGE_TASK_ID} -out=\${SGE_TASK_ID}.genes10 -complete=0 -train=illumina_5 -thread=${NSLOTS}
#DIFF=$( time_diff $START ); echo -e "${NAM}\t${NSLOTS}\tfgs_time_\${SGE_TASK_ID}\t${DIFF}" >> ${RES}/time.logs


EOF

qsub ./fgs_runner
rm ./fgs_runner



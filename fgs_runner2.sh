#!/bin/bash
#$ -j y
#$ -cwd

source ./01-environment
source /bioinf/projects/megx/mg-traits/resources/config_files/config.functions.bash

############################
# run fgs
############################


IN_FASTA_FILE="05-part-${SGE_TASK_ID}.fasta"

${frag_gene_scan} -genome=${IN_FASTA_FILE} -out=${IN_FASTA_FILE}.genes10 -complete=0 -train=illumina_5 -thread="${NSLOTS}"

if [[ $? -ne "0" ]]; then
  email_comm "frag_gene_scan failed: ${frag_gene_scan} -genome=${IN_FASTA_FILE} -out=${IN_FASTA_FILE}.genes10 ..."
  db_error_comm "frag_gene_scan failed. File ${IN_FASTA_FILE}. Job $JOB_ID"
  qdel -u megxnet
  exit 2;
fi  






#!/bin/bash
#$ -cwd 
#$ -j y 
#$ -tc 3

source ./01-environment
source /bioinf/projects/megx/mg-traits/resources/config_files/config.functions.bash

############################
# sina files
############################

IN_FASTA_FILE="06-part-${SGE_TASK_ID}.fasta"
SINA_OUTFILE_SCREEN="06-part-${SGE_TASK_ID}.16S.screen.fasta"
SINA_OUTFILE_ALIGN="06-part-${SGE_TASK_ID}.16S.align.fasta"
SINA_OUTFILE_CLASSIFY="06-part-${SGE_TASK_ID}.16S.classify.fasta"
SINA_SCREEN_LOG="${SINA_LOG_DIR}/06-part-${SGE_TASK_ID}.16S.screen.log"
SINA_ALIGN_LOG="${SINA_LOG_DIR}/06-part-${SGE_TASK_ID}.16S.align.log"
SINA_CLASSIFY_LOG="${SINA_LOG_DIR}/06-part-${SGE_TASK_ID}.16S.classify.log"
SINA_SCREEN_RUN_LOG="${SINA_LOG_DIR}/06-part-${SGE_TASK_ID}.16S.screen.run.log"
SINA_ALIGN_RUN_LOG="${SINA_LOG_DIR}/06-part-${SGE_TASK_ID}.16S.align.run.log"
SINA_CLASSIFY_RUN_LOG="${SINA_LOG_DIR}/06-part-${SGE_TASK_ID}.16S.classify.run.log"

################
# aling 
################

${sina} -i "${IN_FASTA_FILE}" -o "${SINA_OUTFILE_ALIGN}" --intype fasta --ptdb ${sina_seed} --ptport ${SINA_SOCKET} \
        --fs-min 40 --fs-max 40 --fs-req=1 --fs-kmer-no-fast \
        --fs-min-len=50 --fs-req-full=0 --min-idty 60 \
        --meta-fmt comment \
        --show-conf \
        --log-file="${SINA_ALIGN_LOG}" \
         2> "${SINA_ALIGN_RUN_LOG}"        

if [[ $? -ne "0" ]]; then 
  email_comm "SINA alignment file ${IN_FASTA_FILE} failed. Please contact administrator"
  db_error_comm "SINA alignment failed. Please contact administrator"
 # qdel -u megxnet
  exit 2
fi  

               
NUM_RNA_ALIGN=$(grep -c '>' $SINA_OUTFILE_ALIGN)        
if [[ "${NUM_RNA_ALIGN}" -eq 0 ]]; then 
  email_comm "No aligned RNA sequences by sina in file ${IN_FASTA_FILE}"
  db_error_comm "No aligned RNA sequences by sina"
#  qdel -u megxnet
  exit 2;
fi  

##################
# classify         
#################

${sina} -i "${SINA_OUTFILE_ALIGN}" -o "${SINA_OUTFILE_CLASSIFY}" --ptdb ${sina_ref} --ptport ${SINA_SOCKET} \
    --prealigned \
    --meta-fmt comment \
    --search \
    --search-db ${sina_ref} \
    --lca-fields tax_slv \
    --show-conf \
    --log-file="${SINA_CLASSIFY_LOG}" \
    2> "${SINA_CLASSIFY_RUN_LOG}"

if [[ $? -ne "0" ]]; then 
  email_comm "SINA classify file ${IN_FASTA_FILE} failed. Please contact administrator"
  db_error_comm "SINA classify failed. Please contact administrator"
  #qdel -u megxnet
  exit 2
fi
    
NUM_RNA_CLASSIFY=$(grep -c '>' $SINA_OUTFILE_CLASSIFY)        
if [[ "${NUM_RNA_CLASSIFY}" -eq 0 ]]; then 
  email_comm "No classified RNA sequences by sina in file ${IN_FASTA_FILE}"
  db_error_comm "No classified RNA sequences by sina"
  #qdel -u megxnet
  exit 2;
fi  
    
   
###########################################################################################################
###########################################################################################################
########################### DATA BASE COMMUNICATION #######################################################
###########################################################################################################
###########################################################################################################

source ~/.bashrc
source config.bash
source tmp.vars


############################################################################################################
# 1 - fgs check
############################################################################################################

if [[ "${ERROR_FGS}" ]]; then

  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'FragGeneScan failed. Please contact adminitrator.\
  ' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  
  mail -s "mg_traits:${THIS_JOB_ID} subtask ${SGE_TASK_ID} failed" "${mt_admin_mail}" <<EOF
  "${frag_gene_scan}" -genome="${IN_FASTA_FILE}" -out="${IN_FASTA_FILE}".genes10 -complete=0 -train=sanger_10
  exited with RC "${ERROR_FGS}" in job "${JOB_ID}".
EOF
 
  qdel -u megxnet
  exit 2

fi


############################################################################################################
# 2 - sormerna check
############################################################################################################

if [[ "${ERROR_SORTMERNA}" ]]; then

  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'sormerna failed. Please contact adminitrator.\
  ' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  
  mail -s "mg_traits:${THIS_JOB_ID} subtask ${SGE_TASK_ID} failed" "${mt_admin_mail}" <<EOF
  "${sortmerna}" --reads "${SE}" -a "${NSLOTS}" --ref "${DB}"/rRNA_databases/silva-bac-16s-id90.fasta,\
  "${DB}"/index/silva-bac-16s-db:"${DB}"/rRNA_databases/silva-arc-16s-id95.fasta,\
  "${DB}"/index/silva-arc-16s-db:"${DB}"/rRNA_databases/silva-euk-18s-id95.fasta,\
  "${DB}"/index/silva-euk-18s-db \
  --blast 1 --fastx --aligned "${RES}/${NAM}.sortmerna.rDNA" -v --log -m "${MEM}" --best 1
  exited with RC "${ERROR_SORTMERNA}" in job "${JOB_ID}".

EOF
 
  qdel -u megxnet
  exit 2

fi


############################################################################################################
# 3 - sormerna check
############################################################################################################

if [[ "${ERROR_SINA}" ]]; then

  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'sina failed. Please contact adminitrator.\
  ' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  
  mail -s "mg_traits:${THIS_JOB_ID} subtask ${SGE_TASK_ID} failed" "${mt_admin_mail}" <<EOF
  "${sina}" -i "${RES}"/split_smr/spout_\${SGE_TASK_ID} -o "${RES}"/split_smr/\${SGE_TASK_ID}.16S.align.fasta --intype fasta --ptdb "${sina_seed}" --ptport "${SINA_SOCKET}" \
  --fs-min 40 --fs-max 40 --fs-req=1 --fs-kmer-no-fast \
  --fs-min-len=50 --fs-req-full=0 --min-idty 60 \
  --meta-fmt comment \
  --show-conf \
  --log-file="${RES}"/split_smr/\${SGE_TASK_ID}.16S.align.log \
  2> "${RES}"/split_smr/\${SGE_TASK_ID}.16S.align.run.log        
        
  "${sina}" -i "${RES}"/split_smr/\${SGE_TASK_ID}.16S.align.fasta -o "${RES}"/split_smr/\${SGE_TASK_ID}.16S.classify.fasta --ptdb "${sina_ref}" --ptport "${SINA_SOCKET}" \
  --prealigned \
  --meta-fmt comment \
  --search \
  --search-db "${sina_ref}" \
  --lca-fields tax_slv \
  --show-conf \
  --log-file="${RES}"/split_smr/\${SGE_TASK_ID}.16S.classify.log \
  2> "${RES}/split_smr/\${SGE_TASK_ID}.16S.classify.run.log
  
  exited with RC "${ERROR_SINA}" in job "${JOB_ID}".

EOF
 
  qdel -u megxnet
  exit 2

fi


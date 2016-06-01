###########################################################################################################
###########################################################################################################
########################### DATA BASE COMMUNICATION #######################################################
###########################################################################################################
###########################################################################################################

source ~/.bashrc
source config.bash
source tmp.vars



###########################################################################################################
# 1 - Check database connection
###########################################################################################################

if [[ -n "${DB_CONNECTION}" ]]; then
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  Cannot connect to database. Output:
  "${DB_RESULT}"
EOF
  exit 2
fi   

if [[ "${DB_RESULT}" != "UPDATE 1" ]]; then
  echo "sample name '${SAMPLE_LABEL}' is not in database"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  sample name '${SAMPLE_LABEL}' is not in database
  Result: "${DB_RESULT}"
EOF
  exit 2
fi
  

###########################################################################################################
#  2 - Download file from MG URL
###########################################################################################################

# validate MG URL 
if [[ -n "${ERROR_MG_URL}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Not a valid URL' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
Invalid URL.
EOF

 exit 1;

fi


# check if it already exist on our DB
if [[ -n "${ERROR_URLDB}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'The URL $MG_URL has been already succesfully crunched. If the file is different please change the file name.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
The URL "${MG_URL}" has been already succesfully crunched. If the file is different please change the file name.
EOF
  
  exit 1;

fi


# download MG_URL
if [[ -n "${ERROR_MG_DOWNLOAD}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve $MG_URL' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';"\
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  Couldn't retrive "${MG_URL}"
EOF

 exit 1;
  
fi


###########################################################################################################
# 3 -  Validate file
###########################################################################################################

# fasta file check 1
if [[ "${FASTA_ERROR_CODE}" -eq "1" ]]; then 
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = '${MG_URL} is not a valid FASTA file. FASTA validation failed at sequence ${FASTA_BAD_HEADER}.\
  ' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  "${MG_URL}" is not a valid FASTA file. Sequence validation failed.
  "${FASTA_BAD_SEQ}"
EOF

 exit 1;

fi

# fasta file check 2
if [[ "${FASTA_ERROR_CODE}" -eq "2" ]]; then  
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '${MG_URL} is not a valid FASTA file. Sequence ${FASTA_BAD_HEADER} too short.\
  ' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  "${MG_URL}" is not a valid FASTA file. Sequence too short found.
  "${FASTA_BAD_SEQ}"
EOF
 
 exit 1;
 
fi

# fasta file check 3
if [[ "${FASTA_ERROR_CODE}" -ne "0" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '${MG_URL} is not a valid FASTA file.' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  "${MG_URL}" is not a valid FASTA file.
EOF 

 exit 1;
 
fi 


###########################################################################################################
# 4 - Create job directory check
###########################################################################################################

if [[ "$(pwd)" != "${THIS_JOB_TMP_DIR}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'could not access job temp dir' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  Could not access job temp dir "${THIS_JOB_TMP_DIR}"
EOF
    exit 1;
fi

###########################################################################################################
# 5 - Check for utilities and directories
###########################################################################################################

if [[ -n "${ERROR_MESSAGE}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET return_code = 2, error_message = '${ERROR_MESSAGE}' WHERE sample_label = '${SAMPLE_LABEL}' AND id = ${MG_ID};" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:"${JOB_ID}" failed" "${mt_admin_mail}" << EOF
  "${ERROR_MESSAGE}"
EOF
  exit 1;
fi


###########################################################################################################
# 6 - Download data files from SVN
###########################################################################################################

# pfam downlaod check
if [[ -n "${ERROR_PFAM}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve ${PFAM_ACCESSIONS_URL}' WHERE sample_label = '${SAMPLE_LABEL}' AND id = ${MG_ID};" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
Could not retrive pfam_file
EOF

 exit 1;

fi

# tf_file download check
if [[ -n "${ERROR_TF}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve ${TFFILE_URL}' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
Could not retrive tf_file
EOF
  
  exit 1;

fi

# slv file download check
if [[ -n "${ERROR_SLV}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve ${SLV_TAX_URL}' WHERE sample_label = '${SAMPLE_LABEL}' AND id = ${MG_ID};" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
Could not retrive svl_file
EOF

  exit 1;

fi


###########################################################################################################
# 7 - Check for duplicates
###########################################################################################################

# cd_hit_dup check
if [[ -n "${CD_HIT_ERROR_CODE}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '${MG_URL} could not be processed by cd-hit-dup' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
"${cd_hit_dup}" -i "${RAW_FASTA}" -o /dev/null > "${UNIQUE_LOG}
exited with RC "${CD_HIT_ERROR_CODE}" in job "${JOB_ID}".
Files are at: "${FAILED_JOBS_DIR}"/job-"${JOB_ID}"
EOF

 exit 1;

fi

# uniq seqs check
if [[ -n "${ERROR_NUM_UNIQ}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '${MG_URL} contains duplicates. Please provide a pre-processed metagenome.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  "${MG_URL}" contains duplicates.
EOF
 
 exit 1;
 
fi


###########################################################################################################
# 8 - Calculate sequence statistics
###########################################################################################################

# infoseq check
if [[ -n "${ERROR_INFOSEQ}" ]]; then 
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Cannot calculate sequence statistics. Please contact adminitrator.' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  infoseq "${RAW_FASTA}" -only -pgc -length -noheading -auto > "${INFOSEQ_TMPFILE}"
  exited with RC "$?" in job "${JOB_ID}"
  Files are at: "${FAILED_JOBS_DIR}"/job-"${JOB_ID}"
EOF

 exit 1;

fi

# seq_stats.R check
if [[ -n "${SEQ_STATS_ERROR_CODE}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Cannot process sequence statistics. Please contact adminitrator.' WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
"${r_interpreter}" --vanilla --slave
exited with RC "${SEQ_STATS_ERROR_CODE}" in job "${JOB_ID}".
Infoseq script
Files are at: "${FAILED_JOBS_DIR}"/job-"${JOB_ID}"
EOF

  exit 1;
  
fi  


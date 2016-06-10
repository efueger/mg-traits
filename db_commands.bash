###########################################################################################################
###########################################################################################################
########################### DATA BASE COMMUNICATION #######################################################
###########################################################################################################
###########################################################################################################

###########################################################################################################
# write JobID and Hostname to database
###########################################################################################################

<<<<<<< HEAD

$test

=======
>>>>>>> 1394c56232ce53232d1d22957777437e5702cd85
source ~/.bashrc
source config.bash
source tmp.vars

echo "UPDATE mg_traits.mg_traits_jobs SET time_started = now(), job_id = $JOB_ID, cluster_node = '$HOSTNAME' WHERE sample_label = '$SAMPLE_LABEL' AND id = $MG_ID;"
DB_RESULT=$(echo "UPDATE mg_traits.mg_traits_jobs SET time_started = now(), job_id = $JOB_ID, cluster_node = '$HOSTNAME' WHERE sample_label = '$SAMPLE_LABEL' AND id = $MG_ID;" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name)

if [ "$?" -ne "0" ]; then
  echo "Cannot connect to database"
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
Cannot connect to database. Output:
$DB_RESULT
EOF
  exit 2
fi

if [ "$DB_RESULT" != "UPDATE 1" ]; then
	echo "sample name '$SAMPLE_LABEL' is not in database"
	mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
sample name '$SAMPLE_LABEL' is not in database
Result: $DB_RESULT
EOF
  exit 2
fi

###########################################################################################################
# ERROR_MESSAGE
###########################################################################################################

if [[ -n "$ERROR_MESSAGE" ]]; then
	echo "UPDATE mg_traits.mg_traits_jobs SET return_code = 2, error_message = '$ERROR_MESSAGE' WHERE sample_label = '$SAMPLE_LABEL' AND id = $MG_ID;" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
	mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" << EOF
        $ERROR_MESSAGE
EOF
  exit 2
fi

############### this needs to be updated 



###########################################################################################################
# FILE CHECK
###########################################################################################################

# link check
if [[ -n "${ERROR_MG_URL}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Not a valid URL' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
Invalid URL.
EOF
fi


# duplicate DB check
if [[ -n "${ERROR_URLDB}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'The URL $MG_URL has been already succesfully crunched. If the file is different please change the file name.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
The URL "${MG_URL}" has been already succesfully crunched. If the file is different please change the file name.
EOF
fi


# download check
if [[ -n "${ERROR_MG_DOWNLOAD}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve $MG_URL' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
Couldn't retrive "${MG_URL}"
EOF
fi

# fasta file check 1
if [[ "$FASTA_ERROR_CODE" -eq "1" ]]; then 
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = '$MG_URL is not a valid FASTA file. FASTA validation failed at sequence $FASTA_BAD_HEADER.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
$MG_URL is not a valid FASTA file. Sequence validation failed.
$FASTA_BAD_SEQ
EOF
fi

# fasta file check 2
if [[ "$FASTA_ERROR_CODE" -eq "2" ]]; then  
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '$MG_URL is not a valid FASTA file. Sequence $FASTA_BAD_HEADER too short.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
  $MG_URL is not a valid FASTA file. Sequence too short found.
  $FASTA_BAD_SEQ
EOF
fi

# fasta file check 3
if [[ "$FASTA_ERROR_CODE" -ne "0" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '$MG_URL is not a valid FASTA file.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
$MG_URL is not a valid FASTA file.
EOF 
fi 


###########################################################################################################
# THIS_JOB_TMP_DIR CHECK
###########################################################################################################

if [[ "$(pwd)" != "${THIS_JOB_TMP_DIR}" ]]; then
    echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'could not access job temp dir' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
    mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
Could not access job temp dir $THIS_JOB_TMP_DIR
EOF
    exit 2
fi


###########################################################################################################
# SVN DOWNLOADS CHECK
###########################################################################################################

# pfam downlaod check
if [[ -n "${ERROR_PFAM}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve $PFAM_ACCESSIONS_URL' WHERE sample_label = '$SAMPLE_LABEL' AND id = $MG_ID;" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
Could not retrive pfam_file
EOF
fi

# tf_file download check
if [[ -n "${ERROR_TF}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve $TFFILE_URL' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
Could not retrive tf_file
EOF
fi

# slv file download check
if [[ -n "${ERROR_SLV}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = 'Could not retrieve $SLV_TAX_URL' WHERE sample_label = '$SAMPLE_LABEL' AND id = $MG_ID;" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
Could not retrive svl_file
EOF
fi


###########################################################################################################
# PREPROCESSING CHECK
###########################################################################################################

# cd_hit_dup check
if [[ -n "${ERROR_CD_HIT}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '$MG_URL could not be processed by cd-hit-dup' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
"${cd_hit_dup}" -i "${RAW_FASTA}" -o /dev/null > "${UNIQUE_LOG}
exited with RC "${CD_HIT_ERROR_CODE}" in job "${JOB_ID}".
Files are at: "${FAILED_JOBS_DIR}"/job-"${JOB_ID}"
EOF
fi

# uniq seqs check
if [[ -n "${ERROR_NUM_UNIQ}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = '$MG_URL contains duplicates. Please provide a pre-processed metagenome.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
  "${MG_URL}" contains duplicates.
EOF
fi


###########################################################################################################
# SEQUENCE STATS CHECK
###########################################################################################################

# infoseq check
if [[ -n "${ERROR_INFOSEQ}}" ]]; then 
 echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Cannot calculate sequence statistics. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
 mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
infoseq "${RAW_FASTA}" -only -pgc -length -noheading -auto > "${INFOSEQ_TMPFILE}"
exited with RC "$?" in job "${JOB_ID}"
Files are at: "${FAILED_JOBS_DIR}"/job-"${JOB_ID}"
EOF
fi

# seq_stats.R check
if [[ -n "${ERROR_SEQ_STATS}" ]]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Cannot process sequence statistics. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
"${r_interpreter}" --vanilla --slave
exited with RC "${SEQ_STATS_ERROR_CODE}" in job "${JOB_ID}".
Infoseq script
Files are at: "${FAILED_JOBS_DIR}"/job-"${JOB_ID}"
EOF


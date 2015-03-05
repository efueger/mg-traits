#!/bin/bash
set -o pipefail
set -x
set -m

START_TIME=`date +%s.%N`

WORK_DIR=$1

echo "Work dir: $WORK_DIR"

cd $WORK_DIR

if [ ! -e ./00-environment ]; then
	echo "00-environment is missing from work dir: $WORK_DIR"
	exit 2
else
	cat 00-environment
	source ./00-environment
fi

IN_FASTA_FILE="$WORK_DIR/05-part-$SGE_TASK_ID.fasta"

# Gene calling via FGS
echo $frag_gene_scan -genome=$IN_FASTA_FILE -out=$IN_FASTA_FILE.genes10 -complete=0 -train=illumina_5 -threads=${NSLOTS}
$frag_gene_scan -genome=$IN_FASTA_FILE -out=$IN_FASTA_FILE.genes10 -complete=0 -train=illumina_5

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'FragGeneScan failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID subtask $SGE_TASK_ID failed" "$mt_admin_mail" <<EOF
$frag_gene_scan -genome=$IN_FASTA_FILE -out=$IN_FASTA_FILE.genes10 -complete=0 -train=sanger_10
exited with RC $? in job $JOB_ID.
EOF
  qdel -u megxnet
  exit 2
fi


END_TIME=`date +%s.%N`
RUN_TIME=`echo $END_TIME-$START_TIME | bc -l`

# update mg_traits_jobs
echo "UPDATE mg_traits.mg_traits_jobs SET total_run_time = total_run_time + $RUN_TIME, time_protocol = time_protocol || ('$JOB_ID', 'mg_traits_fgs:$SGE_TASK_ID', $RUN_TIME)::mg_traits.time_log_entry WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name

## permantly copying local cluster node files to NFS
#cp ${job_out_dir}/* $THIS_JOB_TMP_DIR

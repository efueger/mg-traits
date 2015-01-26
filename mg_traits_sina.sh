#!/bin/bash
set -x
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

export ARBHOME
export LD_LIBRARY_PATH

IN_FASTA_FILE="$WORK_DIR/05-part-$SGE_TASK_ID.fasta"
SINA_OUTFILE_SCREEN="$WORK_DIR/05-part-$SGE_TASK_ID.16S.screen.fasta"
SINA_OUTFILE_ALIGN="$WORK_DIR/05-part-$SGE_TASK_ID.16S.align.fasta"
SINA_OUTFILE_CLASSIFY="$WORK_DIR/05-part-$SGE_TASK_ID.16S.classify.fasta"
SINA_SCREEN_LOG="$SINA_LOG_DIR/05-part-$SGE_TASK_ID.16S.screen.log"
SINA_ALIGN_LOG="$SINA_LOG_DIR/05-part-$SGE_TASK_ID.16S.align.log"
SINA_CLASSIFY_LOG="$SINA_LOG_DIR/05-part-$SGE_TASK_ID.16S.classify.log"
SINA_SCREEN_RUN_LOG="$SINA_LOG_DIR/05-part-$SGE_TASK_ID.16S.screen.run.log"
SINA_ALIGN_RUN_LOG="$SINA_LOG_DIR/05-part-$SGE_TASK_ID.16S.align.run.log"
SINA_CLASSIFY_RUN_LOG="$SINA_LOG_DIR/05-part-$SGE_TASK_ID.16S.classify.run.log"
SINA_SOCKET=":/tmp/mg_traits_pt_"$(tr -cd '[:alnum:]' < /dev/urandom | fold -w 32 | head -n 1)

cleanup() {
    $sina_arb_pt_server -kill -D$sina_seed -T$SINA_SOCKET
    [[ -f "$SINA_SOCKET" ]] && rm -f "$SINA_SOCKET"
}

trap cleanup EXIT


# 16S identification using SINA
# Screening files for 16S
echo "SINA screening"
$sina -i $IN_FASTA_FILE -o $SINA_OUTFILE_SCREEN --ptdb $sina_seed --ptport $SINA_SOCKET \
    --fs-min 1 --fs-max 1 --fs-req=1 --fs-kmer-no-fast \
    --fs-min-len=50 --fs-req-full=0 --min-idty 60 \
    --log-file=$SINA_SCREEN_LOG \
    --show-conf \ 
    > $SINA_SCREEN_RUN_LOG

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'SINA screening failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID subtask $SGE_TASK_ID failed" "$mt_admin_mail" <<EOF
$sina -i $IN_FASTA_FILE -o $SINA_OUTFILE_SCREEN --ptdb $sina_seed --ptport $SINA_SOCKET \
    --fs-min 1 --fs-max 1 --fs-req=1 --fs-kmer-no-fast \
    --fs-min-len=50 --fs-req-full=0 --min-idty 60 \
    --show-conf \ 
    --log-file=$SINA_SCREEN_LOG $SINA_SCREEN_RUN_LOG
exited with RC $? in job $JOB_ID.
EOF
  qdel -u megxnet
  exit 2
fi

NUM_RNA_SCREEN=$(grep -c '>' $SINA_OUTFILE_SCREEN)
echo "16S RNA SCREENED: $NUM_RNA_SCREEN"
if [ "$NUM_RNA_SCREEN" -ne "0" ]; then
# Align the screened sequences
echo "SINA ALIGNMENT"
$sina -i $SINA_OUTFILE_SCREEN -o $SINA_OUTFILE_ALIGN --ptdb $sina_seed --ptport $SINA_SOCKET \
        --fs-min 40 --fs-max 40 --fs-req=1 --fs-kmer-no-fast \
        --fs-min-len=50 --fs-req-full=0 --min-idty 60 \
        --log-file=$SINA_ALIGN_LOG \
        --meta-fmt comment \
        --show-conf
        > $SINA_ALIGN_RUN_LOG

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'SINA alignment failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID subtask $SGE_TASK_ID failed" "$mt_admin_mail" <<EOF
$sina -i $SINA_OUTFILE_SCREEN -o $SINA_OUTFILE_ALIGN --ptdb $sina_seed --ptport $SINA_SOCKET \
        --fs-min 40 --fs-max 40 --fs-req=1 --fs-kmer-no-fast \
        --fs-min-len=50 --fs-req-full=0 --min-idty 60 \
        --log-file=$SINA_ALIGN_LOG \
        --meta-fmt comment \
        --show-conf \ 
        >$SINA_ALIGN_RUN_LOG
exited with RC $? in job $JOB_ID.
EOF
  qdel -u megxnet
  exit 2
fi
else
echo "NO 16S RNA aligned"
fi

NUM_RNA_ALIGN=$(grep -c '>' $SINA_OUTFILE_ALIGN)
echo "16S RNA ALIGNED: $NUM_RNA_ALIGN"
if [ "$NUM_RNA_SCREEN" -ne "0" ] && [ "$NUM_RNA_ALIGN" -ne "0" ]; then
# Classify the aligned sequences via SINA LCA

echo "SINA CLASSIFY"

$sina -i $SINA_OUTFILE_ALIGN -o $SINA_OUTFILE_CLASSIFY \
    --ptdb $sina_ref --ptport $SINA_SOCKET \
    --prealigned \
    --log-file=$SINA_CLASSIFY_LOG \
    --meta-fmt comment \
    --search \
    --search-db $sina_ref \
    --lca-fields tax_slv \
    --show-conf \ 
    >$SINA_CLASSIFY_RUN_LOG

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'SINA classification failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID subtask $SGE_TASK_ID failed" "$mt_admin_mail" <<EOF
$sina -i $SINA_OUTFILE_ALIGN -o $SINA_OUTFILE_CLASSIFY \
    --ptdb $sina_ref --ptport $SINA_SOCKET \
    --prealigned \
    --log-file=$SINA_CLASSIFY_LOG \
    --meta-fmt comment \
    --search \
    --search-db $sina_ref \
    --lca-fields tax_slv \
    --show-conf \ 
    >$SINA_CLASSIFY_RUN_LOG
exited with RC $? in job $JOB_ID.
EOF
  qdel -u megxnet
  exit 2
fi
else
echo "NO 16S RNA CLASSIFIED"
fi

END_TIME=`date +%s.%N`
RUN_TIME=`echo $END_TIME-$START_TIME | bc -l`

# update mg_traits_jobs
echo "UPDATE mg_traits.mg_traits_jobs SET total_run_time = total_run_time + $RUN_TIME, time_protocol = time_protocol || ('$JOB_ID', 'mg_traits_sina:$SGE_TASK_ID', $RUN_TIME)::mg_traits.time_log_entry WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name

## permantly copying local cluster node files to NFS
#cp ${job_out_dir}/* $THIS_JOB_TMP_DIR

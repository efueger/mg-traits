#!/bin/bash

frag_gene_scan=/bioinf/software/fraggenescan/fraggenescan-1.19/run_FragGeneScan.pl

NAM=$1
NSLOTS=$2
NSEQ=$3

echo $NAM $NSLOTS $NSEQ

#Split original
printf "Splitting file ("${NSEQ}" seqs file)..."
awk -vn="${NSEQ}" 'BEGIN {n_seq=0;partid=1;} /^>/ {if(n_seq%n==0){file=sprintf("05-part-%d.fasta",partid);partid++;} print >> file; n_seq++; next;} { print >> file; }' < "${RAW_FASTA}"

NFILES=$(ls -1 05-part*.fasta | wc -l)
echo "Split into ${NFILES} sub jobs..."


cat > fgs_runner << EOF
#!/bin/bash
#$ -j y
#$ -t 1-${NFILES}
#$ -pe threaded 4
#$ -cwd


############################
# run fgs
############################

START_TIME=\$( date +%s.%N )

${frag_gene_scan} -genome=spout_\${SGE_TASK_ID} -out=\${SGE_TASK_ID}.genes10 -complete=0 -train=illumina_5 -thread="${NSLOTS}"

END_TIME=\$( date +%s.%N )
RUN_TIME=\$( echo \${END_TIME} - \${START_TIME} | bc -l ) 


echo "UPDATE mg_traits.mg_traits_jobs SET total_run_time = total_run_time + \${RUN_TIME}, time_protocol = time_protocol || \
('\${JOB_ID}', 'mg_traits_fgs:\${SGE_TASK_ID}', \${RUN_TIME})::mg_traits.time_log_entry WHERE sample_label = '\${SAMPLE_LABEL}' AND id = '\${MG_ID}';" \
| psql -U \${target_db_user} -h \${target_db_host} -p \${target_db_port} -d \${target_db_name}

EOF





qsub ./fgs_runner
rm ./fgs_runner



#!/bin/bash

NAM=$1
NSLOTS=$2
nSEQ=$3
RES=$4

asplit '^>' ${nSEQ} spout_ < ../${NAM}.sortmerna.rDNA.fasta
NFILES=$(ls spout_* | wc -l)


sina="/bioinf/software/sina/sina-1.3.0rc/sina"
sina_arb_pt_server="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/lib/arb_pt_server"
sina_version="1.2.11"
sina_seed="/local/biodb/mg-traits/sina/ssu_seed_50_26_05_13_cut_t.arb"
sina_seed_version="ssu_seed_50_26_05_13_cut_t"
sina_ref="/local/biodb/mg-traits/sina/ssuref_silva_nr99_115_20_07_13.arb"
sina_ref_version="ssuref_silva_nr99_115_20_07_13.arb"
SINA_SOCKET=":/tmp/mg_traits_pt_"$(tr -cd '[:alnum:]' < /dev/urandom | fold -w 32 | head -n 1)

SGE_TASK_ID=1
cat > sina.runner << EOF
#!/bin/bash
#$ -t 1-${NFILES} 
#$ -cwd 
#$ -j y 
#$ -pe threaded ${NSLOTS}

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
# run sina
############################

#START=$( time_start )

$sina -i $RES/split_smr/spout_\${SGE_TASK_ID} -o $RES/split_smr/\${SGE_TASK_ID}.16S.align.fasta --intype fasta --ptdb $sina_seed --ptport $SINA_SOCKET \
        --fs-min 40 --fs-max 40 --fs-req=1 --fs-kmer-no-fast \
        --fs-min-len=50 --fs-req-full=0 --min-idty 60 \
        --meta-fmt comment \
        --show-conf \
        --log-file=$RES/split_smr/\${SGE_TASK_ID}.16S.align.log \
         2> $RES/split_smr/\${SGE_TASK_ID}.16S.align.run.log        
        
$sina -i $RES/split_smr/\${SGE_TASK_ID}.16S.align.fasta -o $RES/split_smr/\${SGE_TASK_ID}.16S.classify.fasta --ptdb $sina_ref --ptport $SINA_SOCKET \
    --prealigned \
    --meta-fmt comment \
    --search \
    --search-db $sina_ref \
    --lca-fields tax_slv \
    --show-conf \
    --log-file=$RES/split_smr/\${SGE_TASK_ID}.16S.classify.log \
    2> $RES/split_smr/\${SGE_TASK_ID}.16S.classify.run.log

#DIFF=$( time_diff $START ); echo -e "${NAM}\t${NSLOTS}\tsina_time_\${SGE_TASK_ID}\t${DIFF}" >> ${RES}/time.logs    
    
EOF

qsub ./sina.runner
rm ./sina.runner
   
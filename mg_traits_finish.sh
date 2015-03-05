#!/bin/bash
set -x
set -o pipefail

function cleanup {
if [ -d $THIS_JOB_TMP_DIR ];then
mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
fi
}

trap cleanup EXIT

START_TIME=`date +%s.%N`
WORK_DIR=$1
GENEAA="05-gene-aa-seqs"
GENENT="05-gene-nt-seqs"
PFAMDB="06-pfamdb"
PFAMFILERAW="06-pfam-raw"
PFAMFILE="06-pfam"
FUNCTIONALTABLE="06-pfam-functional-table"
CODONCUSP="06-codon.cusp"
TFPERC="06-tfperc"
CLPERC="06-perccl"
AA_TABLE="07-aa-table"
CODON_TABLE="07-codon-table"
ABRATIO_FILE="07-ab-ratio"
NUC_FREQS="07-nuc-freqs"
DINUC_FREQS="07-dinuc-freqs"
ODDS_TABLE="07-odds-table"
SLV_TAX_RAW="08-slv-tax-raw"
SLV_TAX_ORDER="08-slv-tax-order"
SLV_TMP_ORDER="08-slv-tmp-order"
SLV_TMP_UNIQUE="08-slv-tmp-unique"
PCA_CODON_FILE="09-pca-codon"
PCA_AA_FILE="09-pca-aa"
PCA_DINUC_FILE="09-pca-dinuc"
PCA_FUNCTIONAL_FILE="09-pca-functional"
PCA_TAXONOMY_FILE="09-pca-taxonomy"
PCA_CODON_DB="09-pca-codon-db"
PCA_AA_DB="09-pca-aa-db"
PCA_DINUC_DB="09-pca-dinuc-db"
PCA_FUNCTIONAL_DB="09-pca-functional-db"
PCA_TAXONOMY_DB="09-pca-taxonomy-db"

echo "Work dir: $WORK_DIR"

cd $WORK_DIR

if [ ! -e ./00-environment ]; then
	echo "00-environment is missing from work dir: $WORK_DIR"
	exit 2
else
	cat 00-environment
	source ./00-environment
fi

LD_LIBRARY_PATH=/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/lib:/bioinf/software/gcc/gcc-4.9/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
###########################################################################################################
# Check job array results
###########################################################################################################

FAA_RESULTS=$(ls -1 05-part*.faa | wc -l)
FFN_RESULTS=$(ls -1 05-part*.ffn | wc -l)
SLV_CLASSIFY_RESULTS=$(ls -1 05-part-*.screen.fasta | wc -l)
echo "subjobs: $SUBJOBS"
echo "FAA results found: $FAA_RESULTS"
echo "FFN results found: $FFN_RESULTS"

if [ "$FAA_RESULTS" -ne "$SUBJOBS" ] || [ "$FFN_RESULTS" -ne "$SUBJOBS" ] || [ "$SLV_CLASSIFY_RESULTS" -ne "$SUBJOBS" ]; then
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'At least one of the subjobs did not yield results. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$JOB_ID failed" "$mt_admin_mail" <<EOF
At least one of the subjobs did not yield results.
Files are at: $FAILED_JOBS_DIR/job-$JOB_ID
EOF
  exit 2
fi

cat 05-part*.faa > $GENEAA
cat 05-part*.ffn > $GENENT

NUM_GENES=$(grep -c '>' $GENENT)
printf "Number of genes: %d\n" $NUM_GENES

# Check if we have genes if not exit
if [ "$NUM_GENES" -eq "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'No genes predicted in the sample.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID subtask $SGE_TASK_ID failed" "$mt_admin_mail" <<EOF
No genes predicted in the sample
exited with RC $? in job $JOB_ID.
EOF
  exit 2
fi

NUM_RNA=$(cat 05-part-*.classify.fasta | grep -c '>')
printf "Number of genes: %d\n" $NUM_RNA

# Check if we have 16S if not exit
if [ "$NUM_RNA" -eq "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'No RNA found in the sample.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID subtask $SGE_TASK_ID failed" "$mt_admin_mail" <<EOF
No RNA found in the sample
exited with RC $? in job $JOB_ID.
EOF
  exit 2
fi

###########################################################################################################
# Functional annotation
###########################################################################################################

echo "Getting functional profile..."
echo "Running UProC on $NSLOTS cores..."
$uproc -t $NSLOTS -p -l -O 2 -P 3 -o $PFAMFILERAW $uproc_pfam $uproc_model $GENENT 
#exit 0
echo $?
if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'UProC error. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID subtask $SGE_TASK_ID failed" "$mt_admin_mail" <<EOF
$uproc -t $NSLOTS -p -s -P 3 -o $PFAMFILERAW $uproc_pfam $uproc_model $GENENT
exited with RC $? in job $JOB_ID.
EOF
  exit 2
fi

echo "Extracting PFAMS..."
cut -f2,7 -d ',' $PFAMFILERAW > $PFAMFILE 

echo "Calculating functional table..."
$r_interpreter --vanilla --slave <<RSCRIPT
n.genes<-$NUM_GENES
t<-read.table(file = '$PFAMFILE', header = F, stringsAsFactors=F, sep=",")
colnames(t)<-c("seq_id", "pfam_acc")
perc.cl<-(length(unique(t[,1]))/n.genes)*100
t<-subset(t, select = "pfam_acc")
p<-read.table(file = '$PFAM_ACCESSIONS', header = F, stringsAsFactors=F)
colnames(p)<-'pfam_acc'
tf<-read.table(file = '$TFFILE', header = F, stringsAsFactors=F)
colnames(tf)<-'pfam_acc'
t.t<-as.data.frame(table(t$pfam_acc))
colnames(t.t)<-c("pfam_acc", "counts")
t.m<-merge(p, t.t, all = T, by= "pfam_acc")
t.m[is.na(t.m)]<-0
colnames(t.m)<-c("pfam_acc", "counts")
tf.m<-merge(t.t, tf, all = F, by= "pfam_acc")
colnames(tf.m)<-c("pfam_acc", "counts")
perc.tf<-(sum(tf.m[,2])/sum(t.m[,2]))*100
write.table(t.m, file = '$FUNCTIONALTABLE', sep = "\t", row.names = F, quote = F, col.names = F)
write.table(perc.tf, file = '$TFPERC', sep = "\t", row.names = F, quote = F, col.names = F)
write.table(perc.cl, file = '$CLPERC', sep = "\t", row.names = F, quote = F, col.names = F)
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Could not process UPro output. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
$r_interpreter --vanilla --slave
exited with RC $? in job $THIS_JOB_ID.
Functional table script!
Files are at: $FAILED_JOBS_DIR/job-$THIS_JOB_ID
EOF
  exit 2
fi

sort -k1 $FUNCTIONALTABLE | sed -e 's/\t/=>/g' | tr '\n' ',' | sed -e 's/^/\"/' -e 's/,$/\"/' > $PFAMDB

cusp --auto -stdout $GENENT |awk '{if ($0 !~ "*" && $0 !~ /[:alphanum:]/ && $0 !~ /^$/){ print $1,$2,$5}}' > $CODONCUSP
if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'cusp failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
cusp --auto -stdout $GENENT |awk '{if ($0 !~ "*" && $0 !~ /[:alphanum:]/ && $0 !~ /^$/){ print $1,$2,$5}}' > $CODONCUSP
exited with RC $? in job $THIS_JOB_ID.
Files are at: $FAILED_JOBS_DIR/job-$THIS_JOB_ID
EOF
  exit 2
fi

echo "OK"

$r_interpreter --vanilla --slave <<RSCRIPT
codon<-read.table(file = "$CODONCUSP", header = F, stringsAsFactors = F, sep = ' ')
codon<-cbind(codon, codon\$V3/sum(codon\$V3))
colnames(codon)<-c("codon", "aa", "raw", "prop")
codon2<-as.data.frame(t(codon\$prop), stringsAsFactors = F)
colnames(codon2)<-codon\$codon
aa<-aggregate(raw ~ aa, data = codon, sum)
aa<-cbind(aa, (aa\$raw/sum(aa\$raw)))
colnames(aa)<-c("aa", "raw", "prop")
aa2<-as.data.frame(t(aa\$prop))
colnames(aa2)<-aa\$aa
ab<-(aa2\$D + aa2\$E)/(aa2\$H + aa2\$R + aa2\$K)
write.table(aa2, file = "$AA_TABLE", sep = "\t", row.names = F, quote = F, col.names  = T)
write.table(codon2, file = "$CODON_TABLE", sep = "\t", row.names = F, quote = F, col.names  = T)
write(ab, file = "$ABRATIO_FILE")
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'AB-Ratio script failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
$r_interpreter --vanilla --slave
exited with RC $? in job $THIS_JOB_ID.
Files are at: $FAILED_JOBS_DIR/job-$THIS_JOB_ID
EOF
  exit 2
fi

ABRATIO=$(cat $ABRATIO_FILE)
PERCTF=$(cat $TFPERC)
PERCCL=$(cat $CLPERC)
compseq --auto -stdout -word 1 $RAW_FASTA |awk '{if (NF == 5 && $0 ~ /^A|T|C|G/ && $0 !~ /[:alphanum:]/ ){print $1,$2,$3}}' > $NUC_FREQS
if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Compseq for nucleotide freqs failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
compseq --auto -stdout -word 1 $RAW_FASTA |awk '{if (NF == 5 && $0 ~ /^A|T|C|G/ && $0 !~ /[:alphanum:]/ ){print $1,$2,$3}}' > $NUC_FREQS
exited with RC $? in job $THIS_JOB_ID.
Files are at: $FAILED_JOBS_DIR/job-$THIS_JOB_ID
EOF
  exit 2
fi

compseq --auto -stdout -word 2 $RAW_FASTA |awk '{if (NF == 5 && $0 ~ /^A|T|C|G/ && $0 !~ /[:alphanum:]/ ){print $1,$2,$3}}' > $DINUC_FREQS
if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Compseq for dinucleotide freqs failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
compseq --auto -stdout -word 2 $RAW_FASTA |awk '{if (NF == 5 && $0 ~ /^A|T|C|G/ && $0 !~ /[:alphanum:]/ ){print $1,$2,$3}}' > $DINUC_FREQS
exited with RC $? in job $THIS_JOB_ID.
Files are at: $FAILED_JOBS_DIR/job-$THIS_JOB_ID
EOF
  exit 2
fi

$r_interpreter --vanilla --slave <<RSCRIPT
nuc<-read.table(file = "$NUC_FREQS", header = F, stringsAsFactors = F, sep = ' ')
rownames(nuc)<-nuc\$V1
nuc\$V1<-NULL
nuc<-as.data.frame(t(nuc))
dinuc<-read.table(file = "$DINUC_FREQS", header = F, stringsAsFactors = F, sep = ' ')
rownames(dinuc)<-dinuc\$V1
dinuc\$V1<-NULL
dinuc<-as.data.frame(t(dinuc))
#Forward strand f(X) when X={A,T,C,G} in S
fa<-nuc\$A[[2]]
ft<-nuc\$T[[2]]
fc<-nuc\$C[[2]]
fg<-nuc\$G[[2]]
#Frequencies when S + SI = S*; f*(X) when X= {A,T,C,G}
faR<-(fa+ft)/2
fcR<-(fc+fg)/2
fAA <- (dinuc\$AA[[2]] + dinuc\$TT[[2]])/2
fAC <- (dinuc\$AC[[2]] + dinuc\$GT[[2]])/2
fCC <- (dinuc\$CC[[2]] + dinuc\$GG[[2]])/2
fCA <- (dinuc\$CA[[2]] + dinuc\$TG[[2]])/2
fGA <- (dinuc\$GA[[2]] + dinuc\$TC[[2]])/2
fAG <- (dinuc\$AG[[2]] + dinuc\$CT[[2]])/2
pAA <- fAA/(faR * faR)
pAC <- fAC/(faR * fcR)
pCC <- fCC/(fcR * fcR)
pCA <- fCA/(faR * fcR)
pGA <- fGA/(faR * fcR)
pAG <- fAG/(faR * fcR)
pAT <- dinuc\$AT[[2]]/(faR * faR)
pCG <- dinuc\$CG[[2]]/(fcR * fcR)
pGC <- dinuc\$GC[[2]]/(fcR * fcR)
pTA <- dinuc\$TA[[2]]/(faR * faR)
odds<-cbind(pAA, pAC, pCC, pCA, pGA, pAG, pAT, pCG, pGC, pTA)
colnames(odds)<-c("pAA/pTT", "pAC/pGT", "pCC/pGG", "pCA/pTG", "pGA/pTC", "pAG/pCT", "pAT", "pCG", "pGC", "pTA")
write.table(odds, file = "$ODDS_TABLE", sep = "\t", row.names = F, quote = F, col.names  = F)
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Odds table script failed. Please contact adminitrator.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
$r_interpreter --vanilla --slave
exited with RC $? in job $THIS_JOB_ID.
Odds table script!
Files are at: $FAILED_JOBS_DIR/job-$THIS_JOB_ID
EOF
  exit 2
fi

###########################################################################################################
# Getting taxonomic classification from SINA output
###########################################################################################################
echo "Getting raw taxonomic classification from SINA output"
cat 05-part-*.classify.fasta | grep lca_tax_slv | cut -d '=' -f 2 | sort | uniq -c | awk '{print substr($0, index($0, $2))"=>"$1}' | tr '\n' ',' | sed -e 's/ /_/g' -e 's/^/\"/' -e 's/,$/\"/' > $SLV_TAX_RAW

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error getting raw taxonomic classification.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error getting raw taxonomic classification.
EOF
  exit 2
fi
echo "done"

echo "Getting taxonomic classification from SINA output for rank=order"
cat 05-part-*.classify.fasta | grep lca_tax_slv | cut -d '=' -f 2 | awk 'BEGIN{FS=";"}{if (NF > 4) {print $4}}' | sort | uniq -c | awk '{print substr($0, index($0, $2))"=>"$1}' > $SLV_TMP_ORDER

# We need at least 5 orders to proceed
NUM_ORDER=$(wc -l $SLV_TMP_ORDER | cut -f 1 -d ' ')

if [ "$NUM_ORDER" -lt "5" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Taxonomic analysis error: Not enough ORDERS found in the metagenome.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Taxonomic analysis error: Not enough ORDERS found in the metagenome.
EOF
  exit 2
fi
echo "done"


cut -f1 -d "=" $SLV_TMP_ORDER | cat - $SLV_FILE |sort |uniq -u | awk '{print $0"=>0"}' > $SLV_TMP_UNIQUE

cat $SLV_TMP_ORDER $SLV_TMP_UNIQUE |  tr '\n' ',' | sed -e 's/ /_/g' -e 's/^/\"/' -e 's/,$/\"/' > $SLV_TAX_ORDER

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error getting taxonomic classification for rank=order.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error getting taxonomic classification for rank=order.
EOF
  exit 2
fi
echo "done"

# load mg_traits_codon
tail -n1 $CODON_TABLE | awk -vI=$MG_ID -vO=$SAMPLE_LABEL '{print I"\t"O"\t"$0}' | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_codon FROM STDIN"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting CODON results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting CODON results.
EOF
  exit 2
fi

# load mg_traits_aa
tail -n1 $AA_TABLE | awk -vI=$MG_ID -vO=$SAMPLE_LABEL '{print O"\t"$0"\t"I}' | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_aa FROM STDIN"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting AA results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting AA results.
EOF
  exit 2
fi

# load mg_traits_pfam
printf "$MG_ID\t$SAMPLE_LABEL\t" | cat - $PFAMDB | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_functional FROM STDIN CSV delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting functional results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting functional results.
EOF
  exit 2
fi

# load mg_traits_taxonomy
paste $SLV_TAX_ORDER $SLV_TAX_RAW | awk -vI=$MG_ID -vO=$SAMPLE_LABEL '{print I"\t"O"\t"$0}' | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_taxonomy FROM STDIN CSV delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting taxonomy results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting taxonomy results.
EOF
  exit 2
fi

# load mg_traits_dinuc
printf "$SAMPLE_LABEL\t" | cat - $ODDS_TABLE | awk -vI=$MG_ID '{print $0"\t"I}' | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_dinuc FROM STDIN"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting DINUC results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting DINUC results.
EOF
  exit 2
fi

# insert into mg_traits_results
echo "INSERT INTO mg_traits.mg_traits_results (sample_label, gc_content, gc_variance, num_genes, total_mb, num_reads, ab_ratio, perc_tf, perc_classified, id) VALUES ('$SAMPLE_LABEL',$GC,$VARGC, $NUM_GENES, $NUM_BASES, $NUM_READS, $ABRATIO, $PERCTF, $PERCCL, $MG_ID);" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting SIMPLE TRAITS results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting SIMPLE TRAITS results.
EOF
  exit 2
fi


# Get existing data
# For CODON
psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY (select C.* from mg_traits.mg_traits_codon C inner join mg_traits.mg_traits_jobs_public P on C.id = P.id) TO $PCA_CODON_FILE CSV HEADER delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error exporting codon data.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error exporting codon data.
EOF
  exit 2
fi


# We calculate the PCAs if we have more than 30 metagenomes in the database
# For functional and taxonomy we apply hellinger transformation to the data before PCA
if [ $(wc -l $PCA_CODON_FILE | cut -f1 -d ' ') -ge "30" ]; then

# For AA
psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY (select A.* from mg_traits.mg_traits_aa A inner join mg_traits.mg_traits_jobs_public P on A.id = P.id) TO $PCA_AA_FILE CSV HEADER delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error exporting amino acid data.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error exporting amino acid data.
EOF
  exit 2
fi

# For DINUC
psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY (select D.* from mg_traits.mg_traits_aa D inner join mg_traits.mg_traits_jobs_public P on D.id = P.id) TO $PCA_DINUC_FILE CSV HEADER delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error exporting dinucleotide odds-ratio data.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error exporting dinucleotide odds-ratio data.
EOF
  exit 2
fi

# For FUNCTIONAL
psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY (SELECT F.id,(each(functional)).key as key, (each(functional)).value FROM mg_traits.mg_traits_functional F inner join mg_traits.mg_traits_jobs_public P on F.id = P.id order by id, key) TO $PCA_FUNCTIONAL_FILE CSV HEADER delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error exporting functional data.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error exporting functional data.
EOF
  exit 2
fi

# For TAXONOMY
psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY (SELECT T.id,(each(taxonomy_order)).key as key, (each(taxonomy_order)).value FROM mg_traits.mg_traits_taxonomy T inner join mg_traits.mg_traits_jobs_public P on T.id = P.id order by id, key) TO $PCA_TAXONOMY_FILE CSV HEADER delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error exporting taxonomic data.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error exporting taxonomic data.
EOF
  exit 2
fi


# Calculate PCA for CODON
$r_interpreter --vanilla --slave <<RSCRIPT
library(vegan)
t<-read.table(file="$PCA_CODON_FILE", sep = "\t", header = T, stringsAsFactors = F)
t\$sample_label <- NULL
row.names(t)<- t\$id
t\$id <- NULL
pca<-rda(t)
pca.scores<-as.data.frame(scores(pca)\$sites)
pca.scores<-cbind(rownames(pca.scores), "codon-usage", pca.scores, "$MG_ID")
write.table(pca.scores, file = "$PCA_CODON_DB", col.names = F, row.names = F, sep = "\t", quote = F)
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error calculating codon PCA.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error calculating codon PCA.
EOF
  exit 2
fi

# Calculate PCA for AA
$r_interpreter --vanilla --slave <<RSCRIPT
library(vegan)
t<-read.table(file="$PCA_AA_FILE", sep = "\t", header = T, stringsAsFactors = F)
t\$sample_label <- NULL
row.names(t)<- t\$id
t\$id <- NULL
pca<-rda(t)
pca.scores<-as.data.frame(scores(pca)\$sites)
pca.scores<-cbind(rownames(pca.scores), "amino-acid-content", pca.scores, "$MG_ID")
write.table(pca.scores, file = "$PCA_AA_DB", col.names = F, row.names = F, sep = "\t", quote = F)
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error calculating amino acid PCA.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error calculating amino acid PCA.
EOF
  exit 2
fi

# Calculate PCA for DINUC
$r_interpreter --vanilla --slave <<RSCRIPT
library(vegan)
t<-read.table(file="$PCA_DINUC_FILE", sep = "\t", header = T, stringsAsFactors = F)
t\$sample_label <- NULL
row.names(t)<- t\$id
t\$id <- NULL
pca<-rda(t)
pca.scores<-as.data.frame(scores(pca)\$sites)
pca.scores<-cbind(rownames(pca.scores), "di-nucleotide-odds-ratio", pca.scores, "$MG_ID")
write.table(pca.scores, file = "$PCA_DINUC_DB", col.names = F, row.names = F, sep = "\t", quote = F)
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error calculating dinucleotide odds-ratio PCA.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error calculating dinucleotide odds-ratio PCA.
EOF
  exit 2
fi

# Calculate PCA for FUNCTION
$r_interpreter --vanilla --slave <<RSCRIPT
library(vegan)
t<-read.table(file="$PCA_FUNCTIONAL_FILE", sep = "\t", header = T, stringsAsFactors = F)
t\$sample_label <- NULL
t.ids<-as.vector(unique(t\$id))
createTableAbun <- function(X){
t1<-subset(t, id == X)
s <- cbind(X,t(t1\$value))
colnames(s) <- c("id",as.vector(t(t1\$key)))
return(s)
}
ab.list <- lapply(t.ids, createTableAbun)
ab.table<-data.frame(do.call("rbind", ab.list), stringsAsFactors = F)
row.names(ab.table)<- ab.table\$id
ab.table\$id <- NULL
ab.table <- decostand(ab.table, method = "hellinger")
pca<-rda(ab.table)
pca.scores<-as.data.frame(scores(pca)\$sites)
pca.scores<-cbind(rownames(pca.scores), "functional-table", pca.scores, "$MG_ID")
write.table(pca.scores, file = "$PCA_FUNCTIONAL_DB", col.names = F, row.names = F, sep = "\t", quote = F)
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error calculating fucntional table PCA.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error calculating functional table PCA.
EOF
  exit 2
fi

# Calculate PCA for TAXONOMY
$r_interpreter --vanilla --slave <<RSCRIPT
library(vegan)
t<-read.table(file="$PCA_TAXONOMY_FILE", sep = "\t", header = T, stringsAsFactors = F)
t\$sample_label <- NULL
t.ids<-as.vector(unique(t\$id))
createTableAbun <- function(X){
t1<-subset(t, id == X)
s <- cbind(X,t(t1\$value))
colnames(s) <- c("id",as.vector(t(t1\$key)))
return(s)
}
ab.list <- lapply(t.ids, createTableAbun)
ab.table<-data.frame(do.call("rbind", ab.list), stringsAsFactors = F)
row.names(ab.table)<- ab.table\$id
ab.table\$id <- NULL
ab.table <- decostand(ab.table, method = "hellinger")
pca<-rda(ab.table)
pca.scores<-as.data.frame(scores(pca)\$sites)
pca.scores<-cbind(rownames(pca.scores), "taxonomic-table", pca.scores, "$MG_ID")
write.table(pca.scores, file = "$PCA_TAXONOMY_DB", col.names = F, row.names = F, sep = "\t", quote = F)
RSCRIPT

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error calculating taxonomic table PCA.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error calculating taxonomic table PCA.
EOF
  exit 2
fi


# Load codon PCA data in the DB
cat $PCA_CODON_DB | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_pca FROM STDIN CSV delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting codon PCA results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting codon PCA results.
EOF
  exit 2
fi

# Load amino acid PCA data in the DB
cat $PCA_AA_DB | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_pca FROM STDIN CSV delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting amino acid PCA results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting amino acid PCA results.
EOF
  exit 2
fi

# Load dinucleotide odd-ratio PCA data in the DB
cat $PCA_DINUC_DB | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_pca FROM STDIN CSV delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting dinucleotide odds-ratio PCA results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting dinucleotide odds-ratio PCA results.
EOF
  exit 2
fi

# Load functional table PCA data in the DB
cat $PCA_FUNCTIONAL_DB | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_pca FROM STDIN CSV delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting functional table PCA results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting functional table PCA results.
EOF
  exit 2
fi

# Load taxonomic data PCA data in the DB
cat $PCA_TAXONOMY_DB | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "\COPY mg_traits.mg_traits_pca FROM STDIN CSV delimiter E'\t'"

if [ "$?" -ne "0" ]; then
  echo "failed"
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 2, error_message = 'Error inserting taxonomic table PCA results.' WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name
  cd ..
  mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
  mail -s "mg_traits:$THIS_JOB_ID failed" "$mt_admin_mail" <<EOF
Error inserting taxonomic table PCA results.
EOF
  exit 2
fi
fi

END_TIME=`date +%s.%N`
RUN_TIME=`echo $END_TIME-$START_TIME | bc -l`

# update mg_traits_jobs
echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 0, total_run_time = total_run_time + $RUN_TIME, time_protocol = time_protocol || ('$JOB_ID', 'mg_traits_finish', $RUN_TIME)::mg_traits.time_log_entry  WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | psql -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name

cd ..
#rm -rf "job-$THIS_JOB_ID"
mv $THIS_JOB_TMP_DIR $FINISHED_JOBS_DIR

TOTAL_TIME=$( psql -t -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "select (time_finished - time_started) from mg_traits.mg_traits_jobs WHERE sample_label = '$SAMPLE_LABEL' AND id = '$MG_ID';" | tr -d ' ')
if [ "$?" -eq "0" ]; then
  mail -s "mg_traits:$Analysis of $SAMPLE_LABEL with job id $JOB_ID done." "$mt_admin_mail" <<EOF
Analysis of $SAMPLE_LABEL done in $TOTAL_TIME.
EOF
  exit 0
fi

## permantly copying local cluster node files to NFS
#cp ${job_out_dir}/* $THIS_JOB_TMP_DIR

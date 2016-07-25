function cleanup {
if [ -d $THIS_JOB_TMP_DIR ];then
mv $THIS_JOB_TMP_DIR $FAILED_JOBS_DIR
fi
}

#trap cleanup EXIT

START_TIME=`date +%s.%N`
WORK_DIR=$1
GENEAA="05-gene-aa-seqs"
GENENT="05-gene-nt-seqs"
PFAMDB="07-pfamdb"
PFAMFILERAW="07-pfam-raw"
PFAMFILE="07-pfam"
FUNCTIONALTABLE="07-pfam-functional-table"
CODONCUSP="07-codon.cusp"
TFPERC="07-tfperc"
CLPERC="07-perccl"
AA_TABLE="08-aa-table"
CODON_TABLE="08-codon-table"
ABRATIO_FILE="08-ab-ratio"
NUC_FREQS="08-nuc-freqs"
DINUC_FREQS="08-dinuc-freqs"
ODDS_TABLE="08-odds-table"
SLV_TAX_RAW="09-slv-tax-raw"
SLV_TAX_ORDER="09-slv-tax-order"
SLV_TMP_ORDER="09-slv-tmp-order"
SLV_TMP_UNIQUE="09-slv-tmp-unique"
PCA_CODON_FILE="10-pca-codon"
PCA_AA_FILE="10-pca-aa"
PCA_DINUC_FILE="10-pca-dinuc"
PCA_FUNCTIONAL_FILE="10-pca-functional"
PCA_TAXONOMY_FILE="10-pca-taxonomy"
PCA_CODON_DB="10-pca-codon-db"
PCA_AA_DB="10-pca-aa-db"
PCA_DINUC_DB="10-pca-dinuc-db"
PCA_FUNCTIONAL_DB="10-pca-functional-db"
PCA_TAXONOMY_DB="10-pca-taxonomy-db"

LD_LIBRARY_PATH="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/lib:/bioinf/software/gcc/gcc-4.9/lib64:/usr/lib/libgomp.so.1:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH


function email_comm() {
mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
"${1}"
EOF
}

phen=$1
seed=$2
suffix=".n_remove_10000.0.seed_${seed}"
echo "... Running MTAG on ${phen} ..."
gsutil -m cp  gs://nbaya/sexdiff/prs/${phen}.gwas.{fe,}male${suffix}.tsv.bgz ~/data/

python ~/mtag/mtag.py --sumstats ~/data/${phen}.gwas.female${suffix}.tsv.bgz,~/data/${phen}.gwas.male${suffix}.tsv.bgz \
	--p_name pval --out ~/${phen}${suffix}.mtag_def
python ~/mtag/mtag.py --sumstats ~/data/${phen}.gwas.female${suffix}.tsv.bgz,~/data/${phen}.gwas.male${suffix}.tsv.bgz \
	--p_name pval --perfect_gencov --out ~/${phen}${suffix}.mtag_rg1

gzip -f ~/${phen}*.txt

mv ~/${phen}${suffix}.mtag_def_trait_1.txt.gz ~/${phen}.gwas.female${suffix}.mtag_def.tsv.bgz
mv ~/${phen}${suffix}.mtag_def_trait_2.txt.gz ~/${phen}.gwas.male${suffix}.mtag_def.tsv.bgz
mv ~/${phen}${suffix}.mtag_rg1_trait_1.txt.gz ~/${phen}.gwas.both_sexes${suffix}.mtag_rg1.tsv.bgz
rm ~/${phen}${suffix}.mtag_rg1_trait_2.txt.gz
gsutil -m cp ~/${phen}*gz gs://nbaya/sexdiff/prs/



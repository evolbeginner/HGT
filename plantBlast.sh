#! /bin/bash


for i in Arabidopsis_thaliana Populus_trichocarpa Glycine_max Vitis_vinifera Oryza_sativa Sorghum_bicolor \
	Selaginella_moellendorffii Physcomitrella_patens \
	Micromonas_pusilla Volvox_carteri Chlamydomonas_reinhardtii; do
	#blastp -db db/blastp/ref.pep -query prot/$i.fas -out res/blastp/$i.blast8 -evalue 1e-10 -num_threads 24 -outfmt 6
	diamond blastp -d db/diamond/ref.pep.dmnd -q prot/$i.fas -o res/diamond/$i.dmnd8 -p 12 -e 1e-10 -k 1000
done



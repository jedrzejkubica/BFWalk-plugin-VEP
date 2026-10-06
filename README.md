# BFWalk-plugin-VEP

plugin template downloaded from:

https://github.com/Ensembl/VEP_plugins/blob/release/116/pLI.pm



example VCF

`wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`

`wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi`

`bcftools view -G -r 22:17000000-18000000 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o test.vcf.gz`



example scores_MMAF.txt and scores_NOA.txt

/home/kubicaj/BFWalk-validation/BFWalk-output/Bioinformatics-v4-260917/MMAF/scores.txt

/home/kubicaj/BFWalk-validation/BFWalk-output/Bioinformatics-v4-260917/NOA/scores.txt


run VEP on example data:


`vep -i test.vcf.gz --cache /home/nthierry/.vep/ --offline -o test.out`



run vep in the same dir as BFWalk.pm, scores.tsv and uniprot_parsed.tsv

`vep -i test.vcf.gz -o test4.out --cache /home/nthierry/.vep/ --offline --dir_plugins /home/kubicaj/BFWalk-plugin-VEP/ --plugin BFWalk,file=scores_MMAF.tsv,map=uniprot_parsed.tsv,label=MMAF`
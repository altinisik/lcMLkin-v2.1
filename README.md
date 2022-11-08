lcMLkin v2.1 
================
#### implemented for Python3

lcMLkin is a powerful tool to estimate kinship coefficient on
low-coverage genome data. [Current release](https://github.com/kveeramah/lcMLkin_v2) has been implemented on
Python2. Here, I slightly polished the script and implemented it for
Python3. Now, one can use it with both PL or GL according to your called
vcf. The most important difference between v1 (Lipatov *et al.* 2015) and v2 (Žegarac *et al.* 2021) of lcMLkin is the genotyping pipeline. In [v1](https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation), the authors provided a genotyping script, however this is extremely slow. In [v2](https://github.com/kveeramah/lcMLkin_v2), one can easily calculate their own genotype likelihoods. Below, a standard pipeline using `bcftools` has been provided.

## Running lcMLkin v2.1

It is a simple python script. You need a vcf file (zipped or not) and plink files to be calculate AFs (or pre-calculated allele frequencies). In this version, lcMLkin does LD-pruning for each pair. That's why binary plink files are needed to run the script. Other flags are optional. See help output below. 

#### On frequencies

See the original papers for a detailed explanation. To accurately calculate kinship coefficients on low-coverage data, background allele frequencies are important. For ancient DNA data or for a couple of individuals, generally it is not possible to calculate highly accurate allele frequency of the population. Providing background allele frequencies is a work-around to solve this problem. You can use either present-day high-coverage populations or pooled ancient populations (see [Altınışık *et al.* 2022](https://www.science.org/doi/10.1126/sciadv.abo3609) for a possible approach). These populations should not deviate a lot from the population that you analyse (Žegarac *et al.* 2021). In plink files that was provided to `lcmlkinv2.py`, you should have only individuals that you intend to use for allele frequency calculation.

```bash
python lcmlkinv2.py -h
```
```
usage: lcmlkinv2.py [-h] -v VCF -p PLINK [-f [FST]] [-t [NBTHREADS]] [-o OUT]

lcmlkin v2 adapted for python3

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     vcf file must have PL or GL fields.
  -p PLINK, --plink PLINK
                        plink file prefix to calculate allele frequencies.
                        Even if pre-calculated frequencies are the folder
                        (prefix.frq), one should have plink files for LD-
                        pruning.
  -f [FST], --fst [FST]
                        a priori defined fst (see Zegarac 2021)
  -t [NBTHREADS], --nbthreads [NBTHREADS]
                        number of threads
  -o OUT, --out OUT     output file name. if you do not provide a output name,
                        script will generate a name.
```

A sample output:

```
Ind1	Ind2	Z0ag	Z1ag	Z2ag	PI_HATag	nbSNP
Inda	Indb	1.0	0.0	0.0	0.0	1195
Inda	Indc	0.61	0.37	0.03	0.21	236
Inda	Indd	1.0	0.0	0.0	0.0	110
Indb	Indc	1.0	0.0	0.0	0.0	114
Indb	Indd	0.77	0.22	0.01	0.116	471
```


## Genotyping Pipeline

1. First, download 1000G sites from the website. This is needed for correct ref/alt info of positions.

```bash
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz{,.tbi}
```
2. Prepare inputs and arrange variables 
  + `snplist.txt`: tab-separated file with two columns, `CHR` and `Position`. See `vcftools` documentation for more info about the format. You can use pre-selected SNPs, such as 1240K, or you can generate your own positions. To have some idea about the SNP selection, see Altınışık *et al.* 2022 (Material-Methods). 
  + `bam.list`: path to bam files one in each line.

```bash
vcftools=/PATH/TO/vcftools
bcftools=/PATH/TO/bcftools
reffa=/PATH/TO/reference.fa
lcmlkin=/PATH/TO/lcmlkinv2.py
out=outputname
```
3. Subset bi-allelic SNPs which are referred in `snplist.txt` from whole 1000G dataset. Then, index the file.

```bash
$vcftools --gzvcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz --positions snplist.txt --remove-indels --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | bgzip -c > ${out}.vcf.gz &
$bcftools index -f ${out}.vcf.gz &
```

4. Create a file including correct order of ref/alt alleles. 

```bash
$bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${out}.vcf.gz | bgzip -c > ${out}.tsv.gz &
tabix -s1 -b2 -e2 ${out}.tsv.gz
```

5. Calculate genotype likelihoods of selected SNPs in your own data. Since lcMLkin uses info from low-quality SNPs (see [this](https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation/wiki/Best-Practices)), it is better to use low quality thresholds, but you can use your own parameters. Be careful not use these less confident genotype likelihoods in other analyses. 

```bash
$bcftools mpileup -f ${reffa} -B -q20 -Q5 -I -a 'FORMAT/DP' -T {out}.vcf.gz -b bam.list -Ou | \
$bcftools call -Am -C alleles -f GP,GQ -T {out}.tsv.gz -Oz -o mydata.vcf.gz &
```

6. Run lcMLkin v2.1

```bash
python lcmlkinv2.py -v mydata.vcf.gz -p PlinkPrefix 
```

## Other softwares for kinship analysis on low-coverage data

1. [READ](https://bitbucket.org/tguenther/read/src/master/)
2. [TKGWV2](https://github.com/danimfernandes/tkgwv2)
3. [ngsrelate](https://github.com/ANGSD/NgsRelate)


## Citation

If you use this script, please cite original papers. 
- The first version was published in Lipatov *et al.* 2015: https://doi.org/10.1101/023374
- The second version was published in Žegarac *et al.* 2021: https://doi.org/10.1038/s41598-021-89090-x

The SNP list ascertained to Yoruba population used in [Altınışık *et al.* 2022](https://www.science.org/doi/10.1126/sciadv.abo3609) is openly available [here](https://zenodo.org/record/7305608). Please see the paper for a detailed description of SNP list preparation procedure.

# BFWalk-plugin-VEP

This repository contains the [VEP](https://github.com/ensembl/ensembl-vep) plugin for [BFWalk](https://github.com/jedrzejkubica/BFWalk).


## How it works

- **One file per phenotype:** Run BFWalk once per phenotype and pass the scores file with `file=`. Each phenotype gets an output column.
- **Matching:** Running BFWalk on an interactome scores proteins. The plugin matches them to transcripts through the UniProt/Swiss-Prot accessions from the VEP cache.
- **Per-transcript scores:** Only transcripts whose protein has a Swiss-Prot accession get a score. Non-coding transcripts, and coding transcripts mapped only to TrEMBL, get no value. If a transcript has several Swiss-Prot accessions, the highest score is reported.


## Requirements

- [VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) with cache (or database) that includes UniProt annotation
- [BFWalk](https://github.com/jedrzejkubica/BFWalk)


## How to install

```bash
git clone https://github.com/jedrzejkubica/BFWalk-plugin-VEP.git
```

Then copy `BFWalk.pm` into the VEP plugins directory (e.g. `~/.vep/Plugins/`) or point VEP to this repository with `--dir_plugins`.


## How to use

**1. Generate scores with BFWalk** (once per phenotype)

```bash
python BFWalk.py \
    --network interactome.sif \
    --seeds seeds.txt \
    1> scores.tsv
```

The output is a tab-separated file with a header `NODE<tab>SCORE`, one row per node in the network.

**2. Run VEP with the BFWalk plugin**

```bash
# single phenotype
vep -i variants.vcf \
    --cache \
    --offline \
    --plugin BFWalk,file=scores.tsv

# many phenotypes (here two)
vep -i variants.vcf \
    --cache \
    --offline \
    --plugin BFWalk,file=MMAF:scores_MMAF.tsv,file=NOA:scores_NOA.tsv
```


### Parameters

Parameters are passed to the plugin as comma-separated key=value pairs after `--plugin BFWalk`.

| Parameter | Description |
|---|---|
| `file=[PHENOTYPE:]FILE` | BFWalk scores (header `NODE<tab>SCORE`) for `PHENOTYPE`. The name can be omitted when only one file is given. |
| `format=FMT` | `sprintf` format for scores (default `%.3g`) |

### Output

The plugin adds one column per phenotype with scores to each transcript’s annotation in the VEP output.

| Column | Description |
|---|---|
| `BFWalk_<PHENOTYPE>_score` | BFWalk score for the transcript's protein, one column per phenotype |
| `BFWalk_score` | The same, when a single file is given without a name |
 
## How to test

This example annotates a 1 Mb region of chromosome 22 from the 1000 Genomes Project with BFWalk scores for two male infertility phenotypes (MMAF and NOA).

**1. Download test data** (1000 Genomes Project, chromosome 22, GRCh37):
 
```bash
# VCF file
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# index file (.tbi)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
```
 
**2. Extract a small region without genotypes:**
 
```bash
bcftools view -G -r 22:17000000-18000000 \
    ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
    -Oz -o test.vcf.gz
```
 
**3. Get example score files:** place `scores_MMAF.tsv` and `scores_NOA.tsv` in this directory.
 
**4. Run VEP with the plugin:**
 
```bash
vep -i test.vcf.gz \
    --offline \
    --cache --dir_cache <path/to/vep_cache> \
    --dir_plugins <path/to/BFWalk-plugin-VEP> \
    --plugin BFWalk,file=MMAF:scores_MMAF.tsv,file=NOA:scores_NOA.tsv \
    -o test.out \
```
 
> The example data is on GRCh37, so use a GRCh37 VEP cache (add `--assembly GRCh37` if you have more than one installed).
 
This test run produces two files: `test.out` (annotations with `BFWalk_MMAF_score` and `BFWalk_NOA_score`) and `test.out_summary.html` (VEP run summary).
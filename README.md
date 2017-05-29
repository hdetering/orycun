This repository describes the analysis pipeline used in an association study of an Australian rabbit population.

This repository is structured as follows:
```
/
|- meta         - meta information about the population studied
|- results      - output files generated by the analysis steps
|- scripts      - analysis scripts used
```

# Outline

In this study, DNA from 81 individual rabbits was extracted and processed using the Genotyping-by-Sequencing ([Elshire2011][1]) protocol. Subsequently, individuals were pooled and sequenced on an Illumina HiSeq2000 (5 Lanes, 100bp paired-end).

A genome-wide association study was carried out to determine potential genetic factors associated with the survival of infection with rabbit heamorrhagic disease (RHD) virus.

# Quality Control and De-multiplexing

Sequencing reads were trimmed to 90 base pairs (bp) length by removing the last 10 bp using [fastx-trimmer](http://hannonlab.cshl.edu/fastx_toolkit/). This QC approach was chosen because the downstream tool for processing of the loci (Stacks) did accept only reads of uniform length at the time of the analyis.

```bash
cat ./raw/120808_I162_FCC10KDACXX_L3_1.fq | \
    fastx_trimmer -l 90 > ./raw/120808_I162_FCC10KDACXX_L${PBS_ARRAYID}_1.trimmed.fq
cat ./raw/120808_I162_FCC10KDACXX_L$3_2.fq | \
    fastx_trimmer -l 90 > ./raw/120808_I162_FCC10KDACXX_L${PBS_ARRAYID}_2.trimmed.fq
```

Since individuals were pooled for sequencing after tagging them with custom [barcodes](meta/samples_lib_barcode.csv), sequencing reads were de-multiplexed using the [process_radtags](http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) script of the [Stacks](http://catchenlab.life.illinois.edu/stacks/) suite.

```bash
# Quality filter (checking barcodes and restriction sites as well)
# Since we use barcodes of variable length which STACKS does not support
# we need to run a filtering step for each barcode length.
#
# options:
#  -r : rescue barcodes and RAD-tags.
#  -c : clean data, remove any read with an uncalled base.
#  -q : discard reads with low quality scores.
#  -t : truncate final read length to this value.
#  -w : size of the sliding window (fraction of read length).
#  -s : quality score limit. If the average score within the sliding window drops below this value, the read is discarded.
for l in 4 5 6 7 8 9
do
  $HOME/tools/stacks/bin/process_radtags \
    -1 ./raw/120808_I162_FCC10KDACXX_L3_1.trimmed.fq \
    -2 ./raw/120808_I162_FCC10KDACXX_L3_2.trimmed.fq \
    -o ./samples/ -b ./barcodes.$l.txt -e pstI -r -c -q -t 81 -w 0.15 -s 20
done
```

# Read Mapping

Reads were aligned to the reference genome [OryCun2.0](http://www.ensembl.org/Oryctolagus_cuniculus/Info/Index) using [BWA MEM](https://github.com/lh3/bwa). Subsequently, read alignments with quality score below 40 were filtered out.

```bash
# 1. Align paired-end reads for each sample.
# 2. Sort alignments
# 3. Filter by mapping quality (MAQ>=40)
while read sample group; do
  bwa mem ref/ensembl/oryCun2 fastq/${sample}.1.fq.gz fastq/${sample}.2.fq.gz \
  | samtools view -bS - > bam/${sample}.bam
  samtools sort -T ${sample} -O bam bam/${sample}.bam > bam/${sample}.s.bam
  samtools view -b -q 40T bam/${sample}.s.bam > bam/${sample}.s.f.bam
  samtools index bam/${sample}.s.f.bam
done < popmap.txt
```

Inspection of mapping results revealed that three samples (Y2510, pt1125, TF_pt448) had significantly lower mapping rates than the rest of the samples. Those three samples were excluded from further analysis, leaving **78 samples** (**36**: died from RHD, **42**: survivors of RHD).

The resulting BAM files have been uploaded to ENA under accession [PRJEB20958](http://www.ebi.ac.uk/ena/data/view/PRJEB20958).

# SNP Calling

Single nucleotide polymorphisms were identified using the Stacks ref_map.pl pipeline.

```bash
# calculate SNPs from mapped reads
# -T: number of threads
# -m: minimum coverage for each stack
# -B: database name
# -b: batch ID
# -a: batch run date
# -D: batch description
# -S: disable recording SQL data in database
# -O: population map
$HOME/bin/ref_map.pl -T 24 -m 10 \
  -B begendiv_rabbit_radtags -b 1 -a 2016-01-11 -D "78 samples reference-aligned using BWA MEM" \
  -S \
  -O ./popmap.txt \
  -X "populations:--vcf" \
  -X "populations:--plink" \
  -X "populations:--structure" \
  -o ./stacks \
  -s ./stacks/bam/pt1954.s.f.bam \
# [77 more samples...]
```

## Result:

**283602 SNPs** were be identified across 78 individuals.

# Association Test

Tests for association with survival of RHD were carried out using [PLINK](http://zzz.bwh.harvard.edu/plink/). SNPs were required to be located on autosomes and to be genotyped in at least 50% of samples.

```bash
# create phenotype file
awk -v "OFS=\t" '!/^#/{print $1,$2,$1}' batch_1.plink.ped > pheno.txt
# create a list of loci on X chromosome (will be blacklisted in later steps)
awk '$1=="X"{print $2}' batch_1.plink.map > snps_chrX.txt

# calculate association scores for SNPs
# --geno 0.5: filter out SNPs with >50% of individuals having a missing genotype
# --exclude snps_chrX.txt: filter out all SNPs on X chromosome
# --fisher:   perform association test using Fisher's exact test
# --model:    tests for association using all available models
plink --noweb --file batch_1.plink --missing-genotype 0 \
      --pheno pheno.txt \
      --geno 0.5 --exclude snps_chrX.txt \
      --fisher --model
```

## Result

**154546 SNPs** with association scores remain, of which p-values from *Fisher exact tests* under the *genotypic model* ("GENO") were further analyzed.

# Gene Assignment

SNPs were assigned to genes annotated on the OryCun2.0 assembly using the R package [LDsnpR](http://services.cbu.uib.no/software/ldsnpr). Please refer to script [snps2genes.R](scripts/snps2genes.R) for details.

## Result

**9219 SNPs** were assigned to at least one gene, **14962 genes** had at least one SNP assigned. Output can be found in [results/genes_scores_snps.40k.csv](results/genes_scores_snps.40k.csv) (raw) and [results/genes_scores_snps.40k.xls](results/genes_scores_snps.40k.xls) (with additional meta info).

[1]: https://doi.org/10.1371/journal.pone.0019379

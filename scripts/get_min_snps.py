#!/usr/bin/env python3
import math, sys

if len(sys.argv) < 4:
  print("usage: %s snp_scores.csv gene_scores.csv gene_snps.csv" % sys.argv[0], file=sys.stderr)
  sys.exit(0)

fn_snp_scores  = sys.argv[1]
fn_gene_scores = sys.argv[2]
fn_gene_snps   = sys.argv[3]

# read snp scores into map
snp2score = {}
with open(fn_snp_scores) as fp:
  line = fp.readline()
  while line.startswith('#'):
    line = fp.readline()
  while line:
    # [0]SNP_ID, [1]CHR, [2]POS, [3]STRAND, [4]A, [5]B, [6]TEST, [7]AFF, [8]UNAFF, [9]SCORE
    parts = line.strip().split('\t')
    snp_id, score = (parts[0], parts[9])
    snp2score[snp_id] = float(score)
    line = fp.readline()

print("read %d snp scores." % len(snp2score), file=sys.stderr)

# read gene scores into map
gene2score = {}
with open(fn_gene_scores) as fp:
  line = fp.readline()
  while line.startswith('#'):
    line = fp.readline()
  # ok to go to next line, skipping non-commented header
  for line in fp:
    # [0]ROW_ID, [1]CHR, [2]START, [3]END, [4]WIDTH, [5]ENSEMBL_GENE_ID, [6]SCORE
    parts = line.strip().split('\t')
    gene_id, score = (parts[5], parts[6])
    gene2score[gene_id] = float(score) if score != 'NA' else -1

print("read %d gene scores." % len(gene2score), file=sys.stderr)

# identify minimal SNPs
with open(fn_gene_snps) as fp:
  line = fp.readline()
  while line.startswith('#'):
    line = fp.readline()
  # ok to go to next line, skipping non-commented header
  for line in fp:
    # [0]ROW_ID, [1]CHR, [2]START, [3]END, [4]WIDTH, [5]ENSEMBL_GENE_ID, [6]SNP_LIST
    parts = line.strip().split('\t')
    gene_id, snp_list = (parts[5], parts[6])
    min_snp = '-'
    if (snp_list != 'NA'):
      for snp_chr in snp_list.split(';'):
        snp_num = math.floor(float(snp_chr))
        snp_id = math.floor(snp_num/100)
        offset = snp_num % 100
        snp_id = "%s_%d" % (snp_id, offset)
        if gene2score[gene_id] == snp2score[snp_id]:
          min_snp = snp_id
          break
    print("%s\t%s" % (gene_id, min_snp))
    gene2score[gene_id] = float(score) if score != 'NA' else -1

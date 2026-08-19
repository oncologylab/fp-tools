#!/usr/bin/env bash
set -euo pipefail

mkdir -p encode_data/bams encode_data/peaks

curl -fL --retry 3 -o encode_data/bams/ENCFF624SON.bam https://www.encodeproject.org/files/ENCFF624SON/@@download/ENCFF624SON.bam
curl -fL --retry 3 -o encode_data/bams/ENCFF926KFU.bam https://www.encodeproject.org/files/ENCFF926KFU/@@download/ENCFF926KFU.bam
curl -fL --retry 3 -o encode_data/bams/ENCFF990VCP.bam https://www.encodeproject.org/files/ENCFF990VCP/@@download/ENCFF990VCP.bam
curl -fL --retry 3 -o encode_data/bams/ENCFF077FBI.bam https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam
curl -fL --retry 3 -o encode_data/bams/ENCFF128WZG.bam https://www.encodeproject.org/files/ENCFF128WZG/@@download/ENCFF128WZG.bam
curl -fL --retry 3 -o encode_data/bams/ENCFF534DCE.bam https://www.encodeproject.org/files/ENCFF534DCE/@@download/ENCFF534DCE.bam
curl -fL --retry 3 -o encode_data/peaks/ENCFF536RJV.bed.gz https://www.encodeproject.org/files/ENCFF536RJV/@@download/ENCFF536RJV.bed.gz
curl -fL --retry 3 -o encode_data/peaks/ENCFF855PCP.bed.gz https://www.encodeproject.org/files/ENCFF855PCP/@@download/ENCFF855PCP.bed.gz

md5sum -c <<'CHECKSUMS'
4a18e40ee643905cf62ae1fbbd33fb87  encode_data/bams/ENCFF624SON.bam
debb5a616ce1ca26a0d956ce54eae02b  encode_data/bams/ENCFF926KFU.bam
7cc7c6739736b169f56e966784cefcde  encode_data/bams/ENCFF990VCP.bam
135e5a65e5eafac97a8a2242b4c13963  encode_data/bams/ENCFF077FBI.bam
ed3d553646995be008241f6e75f57042  encode_data/bams/ENCFF128WZG.bam
4bd9c1e48617afea2a2cf4834eed961b  encode_data/bams/ENCFF534DCE.bam
4020785aaa90efbebc5a3fc92b5a9cff  encode_data/peaks/ENCFF536RJV.bed.gz
bee81b0efe50b144c24d2060a87aa843  encode_data/peaks/ENCFF855PCP.bed.gz
CHECKSUMS

gzip -dc encode_data/peaks/ENCFF536RJV.bed.gz > encode_data/peaks/ENCFF536RJV.bed
gzip -dc encode_data/peaks/ENCFF855PCP.bed.gz > encode_data/peaks/ENCFF855PCP.bed

for bam in encode_data/bams/*.bam; do
  samtools index "$bam"
done

## 1. Map the reads to reference genome.

```sh
#!/bin/bash
echo "开始时间：`date '+%Y%m%d %H-%M-%S'`"

## index
bowtie2_index=/data/wangchao/reference/hg38/hg38

#需要比对样本
cat sample| while  read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}

echo " start map ...."
bowtie2  -p 15  --very-sensitive -X 2000 -x  $bowtie2_index -1 $fq1 -2 $fq2 |samtools sort  -O bam  -@ 5 -o - > ${sample}.raw.bam
samtools index ${sample}.raw.bam
bedtools bamtobed -i ${sample}.raw.bam  > ${sample}.raw.bed
samtools flagstat ${sample}.raw.bam  > ${sample}.raw.stat

sambamba markdup --overflow-list-size 600000  --tmpdir='./'  -r ${sample}.raw.bam  ${sample}.rmdup.bam
samtools index   ${sample}.rmdup.bam

## Calculate %mtDNA:
mtReads=$(samtools idxstats  ${sample}.rmdup.bam | grep 'chrM' | cut -f 3)
totalReads=$(samtools idxstats  ${sample}.rmdup.bam | awk '{SUM += $3} END {print SUM}')
echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'

samtools flagstat  ${sample}.rmdup.bam > ${sample}.rmdup.stat
samtools view  -h  -f 2 -q 30    ${sample}.rmdup.bam   |grep -v chrM |samtools sort  -O bam  -@ 5 -o - > ${sample}.last.bam
samtools index   ${sample}.last.bam
samtools flagstat  ${sample}.last.bam > ${sample}.last.stat
done

echo "结束时间：`date '+%Y%m%d %H-%M-%S'`"

## 合并不同index的去重复bam文件成一个文件、
nohup samtools merge 293T-output8.bam *last.bam &
```

## 2.Split the chromsome

```shell
# 拆分前建立索引
samtools index 293T-output8.last.bam
# 按染色体拆分并转bw，无需标准化。
for chrom in `seq 1 22` X Y ;
do
samtools view -bh 293T-output5.last.bam chr${chrom} | samtools sort -@ 8 -O bam -o 293T-output5-chr${chrom}.bam;
samtools index 293T-output5-chr${chrom}.bam;
done

for chrom in `seq 1 22` X Y ;
do
bamCoverage --extendReads --ignoreForNormalization chr$chrom -b 293T-output5-chr${chrom}.bam  -o 293T-output5-chr${chrom}.bw
done

#对未标准化的bw进行处理，评估shear pcr map gquad等的影响，下面以output5示例
for chrom in `seq 1 22` X Y;
do
cradle correctBias_stored -ctrlbw /data/wangchao/all_starr-seq/input5/input5/chr/293T-input5-chr${chrom}.bw  -expbw /data/wangchao/all_starr-seq/output5/output5/chr/293T-output5-chr${chrom}.bw  -r /data/wangchao/starr-seq/target/chr${chrom}.target.bed -p 15  -biasType shear pcr map gquad -covariDir /data/wangchao/starr-seq/covariate_files/hg38_fragLen300_kmer100 -genome /data/wangchao/starr-seq/hg38.2bit -o output5_correctBias_bwresult -bl /data/wangchao/starr-seq/hg38.blacklist.bed
done
```

## 3. Callpeack

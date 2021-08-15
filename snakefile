from __future__ import print_function
import os
import pandas as pd

LIB_FD = "data/libraries/"
GENOME = "data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
SAMPLES = os.listdir(LIB_FD)
CONTIGS = pd.read_table(GENOME + ".fai", header=None, usecols=[0], squeeze=True, dtype=str)

rule all:
    input:
        expand(LIB_FD + "{sample}/{sample}_sorted.bam", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.bam", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.log", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.grp", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD_recal.flagstat", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_Z.cov", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_Z.sex", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_SmTRP-PZQ.cov", sample=SAMPLES)

rule alignment:
    input:
        read1=LIB_FD + "{sample}/{sample}_R1.fastq.gz",
        read2=LIB_FD + "{sample}/{sample}_R2.fastq.gz",
        genome=GENOME
    output:
        temp(LIB_FD + "{sample}/{sample}_sorted.bam")
    params:
        rg=r"@RG\tID:{sample}\tPL:illumina\tLB:{sample}\tSM:{sample}"
    shell:
        'sleep $[ ( $RANDOM % 100 )  + 1 ]s ; '
        'bwa mem -t $(nproc) -M -R "{params.rg}" "{input.genome}" "{input.read1}" "{input.read2}" | samtools sort -@8 -o "{output}" -'

rule indexing1:
    input:
        LIB_FD + "{sample}/{sample}_sorted.bam"
    output:
        temp(LIB_FD + "{sample}/{sample}_sorted.bam.bai")
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule mark_duplicates:
    input:
        bam=LIB_FD + "{sample}/{sample}_sorted.bam",
        bai=LIB_FD + "{sample}/{sample}_sorted.bam.bai"
    output:
        bam=temp(LIB_FD + "{sample}/{sample}_sorted_MD.bam"),
        log=protected(LIB_FD + "{sample}/{sample}_sorted_MD.log")
    shell:
        'gatk --java-options "-Xmx2g" MarkDuplicates -I "{input.bam}" -O "{output.bam}" -M "{output.log}" --VALIDATION_STRINGENCY LENIENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP $(ulimit -n)'

rule indexing2:
    input:
        LIB_FD + "{sample}/{sample}_sorted_MD.bam"
    output:
        temp(LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai")
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule bsqr:
    input:
        bam=LIB_FD + "{sample}/{sample}_sorted_MD.bam",
        bai=LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai",
        genome=GENOME,
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        table=temp(LIB_FD + "{sample}/{sample}_sorted_MD.grp"),
        bam=protected(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam")
    params:
        table=LIB_FD + r"{sample}/{sample}_sorted_MD.grp"
    shell:
        'gatk --java-options "-Xmx2g" BaseRecalibrator -R "{input.genome}" -I "{input.bam}" --known-sites "{input.sites}" -O "{output.table}" &&\
         gatk --java-options "-Xmx2g" ApplyBQSR -R "{input.genome}" -I "{input.bam}" --bqsr-recal-file "{params.table}" -O "{output.bam}"'

rule indexing3:
    input:
        LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam"
    output:
        protected(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai")
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule stats:
    input:
        bam=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam",
        bai=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai"
    output:
        protected(LIB_FD + "{sample}/{sample}_sorted_MD_recal.flagstat")
    params:
        spl=LIB_FD + r"{sample}/{sample}_sorted_MD_recal.bam"
    shell:
        'echo -e "File: {params.spl}" > "{output}" ; \
         samtools flagstat "{input.bam}" >> "{output}" ; \
         echo  "\n" >> "{output}"'

rule z_cov:
    input:
        bam=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam",
        bai=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai"
    output:
        protected(LIB_FD + "{sample}/{sample}_Z.cov")
    shell:
        'samtools view -u "{input.bam}" SM_V7_ZW | bedtools genomecov -ibam - -d | grep -w "SM_V7_ZW" | awk \'$3 > 0 {{print}}\' > "{output}"'

rule z_sex:
    input:
        LIB_FD + "{sample}/{sample}_Z.cov"
    output:
        protected(LIB_FD + "{sample}/{sample}_Z.sex")
    shell:
        'scripts/schisto_sex_RD.R --file "{input}"'

rule SmTRP_PZQ_cov:
    input:
        bam=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam",
        bai=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai",
        gff="data/genome/Smp_246790.gff"
    output:
        protected(LIB_FD + "{sample}/{sample}_SmTRP-PZQ.cov")
    shell:
        'samtools depth -a -b <(awk \'$3 == "gene" {{print $1"\t"$4-1"\t"$5}}\' "{input.gff}") "{input.bam}" > "{output}"'

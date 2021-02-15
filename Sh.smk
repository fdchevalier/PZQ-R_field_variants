from __future__ import print_function
import os
import pandas as pd

LIB_FD = "data/libraries/2-Sh/"
GENOME = "data/genome/schistosoma_haematobium.PRJNA78265.WBPS15.genomic.fa"
SAMPLES = os.listdir(LIB_FD)
CONTIGS = pd.read_table(GENOME + ".fai", header=None, usecols=[0], squeeze=True, dtype=str)

rule all:
    input:
        expand(LIB_FD + "{sample}/{sample}_sorted.bam", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.bam", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.log", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai", sample=SAMPLES),
        # expand(LIB_FD + "{sample}/{sample}_sorted_MD.grp", sample=SAMPLES),
        # expand(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam", sample=SAMPLES),
        # expand(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai", sample=SAMPLES),
        # expand(LIB_FD + "{sample}/{sample}_sorted_MD_recal.flagstat", sample=SAMPLES),
        expand(LIB_FD + "{sample}/{sample}_sorted_MD.flagstat", sample=SAMPLES),
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
        bam=protected(LIB_FD + "{sample}/{sample}_sorted_MD.bam"),
        log=protected(LIB_FD + "{sample}/{sample}_sorted_MD.log")
    shell:
        'gatk --java-options "-Xmx2g" MarkDuplicates -I "{input.bam}" -O "{output.bam}" -M "{output.log}" --VALIDATION_STRINGENCY LENIENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP $(ulimit -n)'

rule indexing2:
    input:
        LIB_FD + "{sample}/{sample}_sorted_MD.bam"
    output:
        protected(LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai")
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

# rule bsqr:
#     input:
#         bam=LIB_FD + "{sample}/{sample}_sorted_MD.bam",
#         bai=LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai",
#         genome=GENOME,
#         sites="data/genome/sm_dbSNP_v7.vcf"
#     output:
#         table=temp(LIB_FD + "{sample}/{sample}_sorted_MD.grp"),
#         bam=protected(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam")
#     params:
#         table=LIB_FD + r"{sample}/{sample}_sorted_MD.grp"
#     shell:
#         'gatk --java-options "-Xmx2g" BaseRecalibrator -R "{input.genome}" -I "{input.bam}" --known-sites "{input.sites}" -O "{output.table}" &&\
#          gatk --java-options "-Xmx2g" ApplyBQSR -R "{input.genome}" -I "{input.bam}" --bqsr-recal-file "{params.table}" -O "{output.bam}"'

# rule indexing3:
#     input:
#         LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam"
#     output:
#         protected(LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai")
#     shell:
#         'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
#         'samtools index "{input}"'

# rule stats:
#     input:
#         bam=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam",
#         bai=LIB_FD + "{sample}/{sample}_sorted_MD_recal.bam.bai"
#     output:
#         protected(LIB_FD + "{sample}/{sample}_sorted_MD_recal.flagstat")
#     params:
#         spl=LIB_FD + r"{sample}/{sample}_sorted_MD_recal.bam"
#     shell:
#         'echo -e "File: {params.spl}" > "{output}" ; \
#          samtools flagstat "{input.bam}" >> "{output}" ; \
#          echo  "\n" >> "{output}"'

rule stats:
    input:
        bam=LIB_FD + "{sample}/{sample}_sorted_MD.bam",
        bai=LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai"
    output:
        protected(LIB_FD + "{sample}/{sample}_sorted_MD.flagstat")
    params:
        spl=LIB_FD + r"{sample}/{sample}_sorted_MD.bam"
    shell:
        'echo -e "File: {params.spl}" > "{output}" ; \
         samtools flagstat "{input.bam}" >> "{output}" ; \
         echo  "\n" >> "{output}"'

rule SmTRP_PZQ_cov:
    input:
        bam=LIB_FD + "{sample}/{sample}_sorted_MD.bam",
        bai=LIB_FD + "{sample}/{sample}_sorted_MD.bam.bai",
        gff=r"data/genome/MS3_0012599.gff"
    output:
        protected(LIB_FD + "{sample}/{sample}_SmTRP-PZQ.cov")
    shell:
        'samtools depth -a -b <(awk \'$3 == "gene" {{print $1"\t"$4-1"\t"$5}}\' "{input.gff}") "{input.bam}" > "{output}"'


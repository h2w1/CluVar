
configfile: "config/config.yaml"

OUTDIR = f"results/{config['run_id']}"   


rule all:
    input:
        f"{OUTDIR}/autoencoder/real_pred_labels.csv",
        f"{OUTDIR}/autoencoder/cleanup_done.flag"

 
from datetime import datetime

onstart:
    print(f"Workflow started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

onsuccess:
    print(f"Workflow completed successfully at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

onerror:
    print(f"Error occurred during workflow at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


rule make_cell_specific_bam:
    input:
        bam=f"{OUTDIR}/raw/{{sample}}.bam",
        barcodes=f"{OUTDIR}/raw/{{sample}}.barcodes.tsv"
    output:
        bam=temp(f"{OUTDIR}/filtered/{{sample}}.filtered.bam"),
        log=f"{OUTDIR}/logs/filtered_bam/{{sample}}.log"
    conda:
        "envs/python.yaml"
    threads: 2
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.log})
        python scripts/filter_bam_by_barcode.py \
            -i {input.bam} \
            -o {output.bam} \
            -b {input.barcodes} \
            &> {output.log}
        """

rule modify_bam_cb_tag:
    input:
        bam=f"{OUTDIR}/filtered/{{sample}}.filtered.bam"
    output:
        bam=temp(f"{OUTDIR}/filtered/{{sample}}.tagged.bam"),
        log=f"{OUTDIR}/logs/modify_cb_tag/{{sample}}.log"
    threads: 2
    singularity:
        "bio_pipeline.sif"
    shell:
        """
        mkdir -p $(dirname {output.log})
        samtools view -h {input.bam} | \
        awk -v sample="{wildcards.sample}" 'BEGIN {{OFS="\\t"}}
        {{
            if ($1 ~ /^@/) {{
                print $0
            }} else {{
                for (i=12; i<=NF; i++) {{
                    if ($i ~ /^CB:Z:/) {{
                        split($i, a, ":")
                        split(a[3], b, "-")
                        $i = "CB:Z:" b[1] "-" sample
                    }}
                }}
                print $0
            }}
        }}' | samtools view -b > {output.bam} 2> {output.log}
        """


rule modify_barcode_tag:
    input:
        barcodes=f"{OUTDIR}/raw/{{sample}}.barcodes.tsv"
    output:
        barcodes=f"{OUTDIR}/filtered/{{sample}}.barcodes.tsv"
    log:
        f"{OUTDIR}/logs/modify_barcode_tag/{{sample}}.log"
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.barcodes})
        mkdir -p $(dirname {log})
        sed 's/-1$/-{wildcards.sample}/' {input.barcodes} > {output.barcodes} 2> {log}
        """



rule bam_to_fastq:
    input:
        bam=f"{OUTDIR}/filtered/{{sample}}.tagged.bam",
        barcodes=f"{OUTDIR}/filtered/{{sample}}.barcodes.tsv"
    output:
        fastq=temp(f"{OUTDIR}/fastq/{{sample}}.fastq")
    log:
        f"{OUTDIR}/logs/bam_to_fastq/{{sample}}.log"
    conda:
        "envs/python.yaml"
    threads: 2
    shell:
        """
        mkdir -p $(dirname {output.fastq})
        mkdir -p $(dirname {log})
        python scripts/renamer.py \
            -f {input.bam} \
            -b {input.barcodes} \
            -o {output.fastq} \
            &> {log}
        """

rule minimap2_remap:
    input:
        fastq=f"{OUTDIR}/fastq/{{sample}}.fastq"
    output:
        bam=temp(f"{OUTDIR}/mapped/{{sample}}.bam")
    log:
        f"{OUTDIR}/logs/minimap2/{{sample}}.log"
    conda:
        "envs/minimap2.yaml"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        minimap2 -ax splice \
            -G50k -k 21 -w 11 --sr \
            -A2 -B8 -O12,32 -E2,1 \
            -r200 -p.5 -N20 -f1000,5000 -n2 \
            -m20 -s40 -g2000 -2K50m --secondary=no \
            -t {threads} \
            ./reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
            {input.fastq} | samtools view -b -o {output.bam} \
            &> {log}
        """


rule retag_bam:
    input:
        bam=f"{OUTDIR}/mapped/{{sample}}.bam"
    output:
        bam=temp(f"{OUTDIR}/mapped/{{sample}}.retagged.bam")
    log:
        f"{OUTDIR}/logs/retag_bam/{{sample}}.log"
    conda:
        "envs/python.yaml"
    threads: 2
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        python scripts/retag.py \
            -i {input.bam} \
            -o {output.bam} \
            &> {log}
        """


rule sort_bam:
    input:
        bam=f"{OUTDIR}/mapped/{{sample}}.retagged.bam"
    output:
        bam=temp(f"{OUTDIR}/processed/{{sample}}.sorted.bam")
    log:
        f"{OUTDIR}/logs/sort_bam/{{sample}}.log"
    singularity:
        "bio_pipeline.sif"
    threads: 2
    resources:
        mem_mb=8000
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        java -jar /opt/tools/picard.jar SortSam \
            -I {input.bam} \
            -O {output.bam} \
            -SORT_ORDER coordinate \
            &> {log}
        """


rule add_read_groups:
    input:
        bam=f"{OUTDIR}/processed/{{sample}}.sorted.bam"
    output:
        bam=temp(f"{OUTDIR}/processed/{{sample}}.sorted.rg.bam")
    log:
        f"{OUTDIR}/logs/add_rg/{{sample}}.log"
    singularity:
        "bio_pipeline.sif"
    threads: 2
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        java -jar /opt/tools/picard.jar AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -RGID {wildcards.sample} \
            -RGLB transcriptome \
            -RGPL ILLUMINA \
            -RGPU machine \
            -RGSM {wildcards.sample} \
            &> {log}
        """


rule mark_duplicates:
    input:
        bam=f"{OUTDIR}/processed/{{sample}}.sorted.rg.bam"
    output:
        deduped=temp(f"{OUTDIR}/processed/{{sample}}.sorted.rg.dedup.bam"),
        metrics=temp(f"{OUTDIR}/processed/{{sample}}.dedup.metrics.txt")
    log:
        f"{OUTDIR}/logs/markdup/{{sample}}.log"
    singularity:
        "bio_pipeline.sif"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.deduped})
        mkdir -p $(dirname {log})
        java -jar /opt/tools/picard.jar MarkDuplicates \
            I={input.bam} \
            O={output.deduped} \
            M={output.metrics} \
            &> {log}
        """


rule split_ncigar:
    input:
        bam=f"{OUTDIR}/processed/{{sample}}.sorted.rg.dedup.bam"
    output:
        bam=temp(f"{OUTDIR}/processed/{{sample}}.split.bam")
    log:
        f"{OUTDIR}/logs/split_ncigar/{{sample}}.log"
    singularity:
        "bio_pipeline.sif"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        /opt/tools/gatk/gatk SplitNCigarReads \
            -R ./reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
            -I {input.bam} \
            -O {output.bam} \
            &> {log}
        """


rule base_recal:
    input:
        bam=f"{OUTDIR}/processed/{{sample}}.split.bam",
        known="reference/ref_vcf_human/dbSNP_GRCh38_human.vcf.gz"
    output:
        recal_table=temp(f"{OUTDIR}/processed/{{sample}}.recal_data.table")
    log:
        f"{OUTDIR}/logs/bqsr/{{sample}}.log"
    singularity:
        "bio_pipeline.sif"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.recal_table})
        mkdir -p $(dirname {log})
        /opt/tools/gatk/gatk BaseRecalibrator \
            -I {input.bam} \
            -R ./reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
            --known-sites {input.known} \
            -O {output.recal_table} \
            &> {log}
        """

rule apply_bqsr:
    input:
        bam=f"{OUTDIR}/processed/{{sample}}.split.bam",
        bqsr=f"{OUTDIR}/processed/{{sample}}.recal_data.table"
    output:
        bam=f"{OUTDIR}/processed/{{sample}}.recal.bam"
    log:
        f"{OUTDIR}/logs/apply_bqsr/{{sample}}.log"
    singularity:
        "bio_pipeline.sif"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        /opt/tools/gatk/gatk ApplyBQSR \
            -I {input.bam} \
            --bqsr-recal-file {input.bqsr} \
            -R ./reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
            -O {output.bam} \
            &> {log}
        """

rule haplotype_caller:
    input:
        bam=f"{OUTDIR}/processed/{{sample}}.recal.bam"
    output:
        vcf=f"{OUTDIR}/haplotypecaller/{{sample}}.vcf"
    log:
        f"{OUTDIR}/logs/haplotypecaller/{{sample}}.log"
    singularity:
        "bio_pipeline.sif"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        /opt/tools/gatk/gatk HaplotypeCaller \
            -R ./reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
            -I {input.bam} \
            -O {output.vcf} \
            -ERC GVCF \
            --dont-use-soft-clipped-bases true \
            -stand-call-conf 20 \
            &> {log}
        """
rule generate_map_file:
    input:
        vcf=expand(f"{OUTDIR}/haplotypecaller/{{sample}}.vcf", sample=config["samples"])
    output:
        mapfile=f"{OUTDIR}/genomics_DB/map.txt"
    run:
        import os

        samples = config["samples"]
        outdir = OUTDIR
        map_path = output.mapfile
        map_dir = os.path.dirname(map_path)
        tmp_dir = f"{outdir}/genomics_DB/tmp"

        # tmp 디렉토리와 map 디렉토리 생성
        os.makedirs(map_dir, exist_ok=True)
        os.makedirs(tmp_dir, exist_ok=True)

        # map.txt 생성
        with open(map_path, "w") as f:
            for sample in samples:
                vcf_path = f"{outdir}/haplotypecaller/{sample}.vcf"
                f.write(f"{sample}\t{vcf_path}\n")

rule makegenomics_DB:
    input:
        map_file=f"{OUTDIR}/genomics_DB/map.txt"
    output:
        db_dir=directory(f"{OUTDIR}/genomics_DB/sample_genomics_DB")
    log:
        f"{OUTDIR}/logs/makegenomics_DB/log"
    singularity:
        "bio_pipeline.sif"
    threads: 3
    shell:
        """
        mkdir -p $(dirname {log})
        /opt/tools/gatk/gatk GenomicsDBImport \
            --genomicsdb-workspace-path {output.db_dir} \
            --batch-size 0 \
            --sample-name-map {input.map_file} \
            --tmp-dir {OUTDIR}/genomics_DB/tmp \
            --intervals ./intervals.list \
            &> {log}
        """

rule GenotypeGVCFs:
    input:
        gendb=f"{OUTDIR}/genomics_DB/sample_genomics_DB"
    output:
        vcf=f"{OUTDIR}/genomics_DB/genotypeGVCFs.vcf"
    log:
        f"{OUTDIR}/logs/GenotypeGVCFs/log"
    singularity:
        "bio_pipeline.sif"
    threads: 3
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.vcf})
        mkdir -p {OUTDIR}/genotypeGVCFs_tmp
        /opt/tools/gatk/gatk GenotypeGVCFs \
            -R ./reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
            -V gendb://{input.gendb} \
            -O {output.vcf} \
            --tmp-dir {OUTDIR}/genotypeGVCFs_tmp \
            &> {log}
        """

rule normalize_and_filter:
    input:
        vcf     = f"{OUTDIR}/genomics_DB/genotypeGVCFs.vcf",
        pon     = "variants_DB/PoN.scRNAseq.hg38.tsv",
        editing = "variants_DB/AllEditingSites.hg38.txt",
        gpon    = "variants_DB/1000g_pon.hg38.vcf"
    output:
        vcf = f"{OUTDIR}/vartrix/filtered_genotypeGVCFs.vcf"
    log:
        f"{OUTDIR}/logs/vartrix/filtered_vcf_log"
    singularity:
        "bio_pipeline.sif"
    threads: 1
    shell: r"""
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        mkdir -p {OUTDIR}/tmp_norm

        set -e

        # Normalize multi-allelics
        bcftools norm -m - -O v -o {OUTDIR}/tmp_norm/norm.vcf {input.vcf}

        # Create BED files from exclusion sources
        awk '$0!~/^#/ {{print $1"\t"$2-1"\t"$2}}' {input.pon}     > {OUTDIR}/tmp_norm/ex_pon.bed
        awk '$0!~/^#/ {{print $1"\t"$2-1"\t"$2}}' {input.editing} > {OUTDIR}/tmp_norm/ex_editing.bed
        awk '$0!~/^#/ && NF >= 5 {{print $1"\t"$2-1"\t"$2}}' {input.gpon} > {OUTDIR}/tmp_norm/ex_gpon.bed

        # Combine exclusion regions
        cat {OUTDIR}/tmp_norm/ex_pon.bed \
            {OUTDIR}/tmp_norm/ex_editing.bed \
            {OUTDIR}/tmp_norm/ex_gpon.bed > {OUTDIR}/tmp_norm/exclude_all.bed

        # Apply filtering
        bcftools view -O v -o {output.vcf} {OUTDIR}/tmp_norm/norm.vcf \
            -T ^{OUTDIR}/tmp_norm/exclude_all.bed \
            &> {log}

        # Optional cleanup
        # rm -rf {OUTDIR}/tmp_norm
    """

rule merge_bam_barcode:
    input:
        bam=expand(f"{OUTDIR}/processed/{{sample}}.recal.bam", sample=config["samples"]),
        barcode=expand(f"{OUTDIR}/filtered/{{sample}}.barcodes.tsv", sample=config["samples"])
    output:
        bam=f"{OUTDIR}/vartrix/merged.bam",
        bai=f"{OUTDIR}/vartrix/merged.bam.bai",
        barcode=f"{OUTDIR}/vartrix/merged.barcodes.tsv"
    log:
        f"{OUTDIR}/logs/vartrix/merge.log"
    threads: 4
    singularity:
        "bio_pipeline.sif"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})

        # BAM 병합
        samtools merge -@ {threads} {output.bam} {input.bam} 2>> {log}

        # 인덱스 생성
        samtools index {output.bam} 2>> {log}

        # barcode 병합
        cat {input.barcode} > {output.barcode}
        """



rule vartrix:
    input:
        bam=f"{OUTDIR}/vartrix/merged.bam",
        barcode=f"{OUTDIR}/vartrix/merged.barcodes.tsv",
        vcf=f"{OUTDIR}/vartrix/filtered_genotypeGVCFs.vcf"
    output:
        mtx=f"{OUTDIR}/vartrix/vartrix.mtx"
    log:
        f"{OUTDIR}/logs/vartrix/vartrix.log"
    singularity:
        "bio_pipeline.sif"
    threads: 10
    shell:
        """
        mkdir -p $(dirname {output.mtx})
        mkdir -p $(dirname {log})

        /opt/tools/vartrix \
            --bam {input.bam} \
            --cell-barcodes {input.barcode} \
            --fasta ./reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
            --vcf {input.vcf} \
            --out-matrix {output.mtx} \
            --scoring-method consensus \
            --threads {threads} \
            --umi \
            --mapq 30 \
            &> {log}
        """

rule filter_vartrix:
    input:
        barcodes     = f"{OUTDIR}/vartrix/merged.barcodes.tsv",
        vcf          = f"{OUTDIR}/vartrix/filtered_genotypeGVCFs.vcf",
        vartrix_mtx  = f"{OUTDIR}/vartrix/vartrix.mtx"
    output:
        filt_mtx = f"{OUTDIR}/vartrix/custom_filtered.mtx",
        var_idx  = f"{OUTDIR}/vartrix/filtered_variants.idx",
        cell_idx = f"{OUTDIR}/vartrix/filtered_cells.idx"
    params:
        top_n    = 300,
        cell_thr = 0.95
    log:
        f"{OUTDIR}/logs/vartrix/filtered_vartrix.log"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {log})
        python scripts/filter_vartrix.py \
            --barcodes {input.barcodes} \
            --vcf {input.vcf} \
            --vartrix-mtx {input.vartrix_mtx} \
            --top-n {params.top_n} \
            --cell-thr {params.cell_thr} \
            --out-dir {OUTDIR}/vartrix \
            > {log} 2>&1
        """


rule convert_mtx:
    input:
        mtx = f"{OUTDIR}/vartrix/custom_filtered.mtx"
    output:
        tsv = f"{OUTDIR}/autoencoder/filtered_matrix.tsv",
        log = f"{OUTDIR}/logs/autoencoder/convert_mtx.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mkdir -p $(dirname {output.log})

        python scripts/convert_mtx_to_tsv.py \
            --mtx {input.mtx} \
            --out {output.tsv} \
            > {output.log} 2>&1
        """

rule search_hyperparams:
    input:
        tsv = f"{OUTDIR}/autoencoder/filtered_matrix.tsv"
    output:
        best_params = f"{OUTDIR}/autoencoder/best_hyperparams.json",
        log = f"{OUTDIR}/logs/autoencoder/search_hyperparams.log"
    params:
        data_type = "real",           # 또는 "sim" 필요 시 변경
        num_clusters = 5              # sim일 경우에만 사용됨
    conda: 
        "envs/ae.yaml"
    shell:
        """
        mkdir -p $(dirname {output.log})

        python scripts/search_hyperparams.py \
            --data_type {params.data_type} \
            --num_clusters {params.num_clusters} \
            --input_tsv {input.tsv} \
            --output_json {output.best_params} \
            > {output.log} 2>&1

        """


rule train_model:
    input:
        hp = f"{OUTDIR}/autoencoder/best_hyperparams.json",
        filtered_matrix = f"{OUTDIR}/autoencoder/filtered_matrix.tsv"
    output:
        decoder_dir = directory(f"{OUTDIR}/autoencoder/decoder_outputs"),
        labels      = f"{OUTDIR}/autoencoder/real_pred_labels.csv"
    
    log:
        f"{OUTDIR}/logs/autoencoder/train_model.log"
    params:
        data_type    = "real",
        num_clusters = 7,
        out_dir      = f"{OUTDIR}/autoencoder"
    conda:
        "envs/ae.yaml"
    shell:
        """
        python scripts/train_wrapper.py \
            --hp {input.hp} \
            --data_type {params.data_type} \
            --num_clusters {params.num_clusters} \
            --data_path {input.filtered_matrix} \
            --out_path {params.out_dir}
        """



rule cleanup_temp_dirs:
    input:
        # train_model의 output을 명시 → train_model 이후에만 실행됨
        decoder_dir = directory(f"{OUTDIR}/autoencoder/decoder_outputs"),
        labels      = f"{OUTDIR}/autoencoder/real_pred_labels.csv"
    output:
        flag = f"{OUTDIR}/autoencoder/cleanup_done.flag"
    run:
        import shutil
        import os

        temp_dirs = [
            f"{OUTDIR}/tmp_norm",
            f"{OUTDIR}/genotypeGVCFs_tmp",
            f"{OUTDIR}/fastq",
            f"{OUTDIR}/genomics_DB/sample_genomics_DB",
            f"{OUTDIR}/genomics_DB/tmp"
        ]

        for d in temp_dirs:
            if os.path.exists(d):
                print(f"Deleting: {d}")
                shutil.rmtree(d)
            else:
                print(f"Directory not found, skipping: {d}")

        with open(output.flag, "w") as f:
            f.write("Cleanup completed.\n")

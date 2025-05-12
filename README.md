# CluVar Snakemake Pipeline

This pipeline performs RNA-seq-based variant calling using GATK. It includes:

- Cell barcode filtering
- BAM remapping via minimap2
- Read group and duplicate handling via Picard
- SplitNCigarReads and BQSR using GATK
- Variant calling via HaplotypeCaller

## How to Use

### Step 1: Create your config.yaml

```yaml
samples:
  - Sample1
  - Sample2
```

### Step 2: Run the pipeline

```bash
snakemake --cores 12 --use-conda --configfile config.yaml
```

Logs will be stored in the `logs/` directory. Temporary files like FASTQ will be automatically cleaned up.

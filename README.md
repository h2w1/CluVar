# CluVar overall pipeline
![image](https://github.com/user-attachments/assets/5915c4be-ac88-45a0-9a3b-7a3567ef3f0d)

# CluVar detail pipeline 
<img width="350" height="450" alt="image" src="https://github.com/user-attachments/assets/35db6b8c-ef74-40c3-831f-02d5ea1910c6" />




## 1. Setting Environment

Create a conda virtual environment for CluVar.

```bash
conda create -n CluVar-env -c conda-forge mamba snakemake python=3.10
conda activate CluVar-env 
```

clone the recent version of this repository.
```bash
git clone https://github.com/h2w1/CluVar.git
```
please download below directories and file on main CluVar directory and keep their name for Snakefile. 

1) reference : human reference genome, reference SNP (dbSNP) 
2) variants_DB : filtering normal, germline variants
3) results : there is test data input now. if you run your data, all results are saved as run_id directory in results.
4) bio_pipeline.sif : singularity image file
5) intervals.list :  chromosome number using for variants calling  
from https://drive.google.com/drive/folders/1bsQ4hoIcGslJMSlu_IYB5O3j5I-WFJ70?usp=drive_link

After download files, your CluVar directory is below. please check below tree.

```
├── bio_pipeline.sif
├── config
│   └── config_test.yaml
├── envs
│   ├── ae.yaml
│   ├── minimap2.yaml
│   └── python.yaml
├── intervals.list
├── README.md
├── reference
│   ├── refdata-gex-GRCh38-2024-A
│   │   ├── fasta
│   │   │   ├── genome.dict
│   │   │   ├── genome.fa
│   │   │   └── genome.fa.fai
│   │   └── reference.json
│   └── ref_vcf_human
│       ├── dbSNP_GRCh38_human.vcf.gz
│       ├── dbSNP_GRCh38_human.vcf.gz.tbi
│       └── editvcf.human.sh
├── results
│   └── test
│       └── raw
│           ├── test.bam
│           └── test.barcodes.tsv
├── scripts
│   ├── convert_mtx_to_tsv.py
│   ├── dataset.py
│   ├── filter_bam_by_barcode.py
│   ├── filter_vartrix.py
│   ├── main.py
│   ├── model.py
│   ├── renamer.py
│   ├── requirements.txt
│   ├── retag.py
│   ├── search_hyperparams.py
│   ├── train.py
│   ├── train_wrapper.py
│   └── utils.py
├── Snakefile
└── variants_DB
    ├── 1000g_pon.hg38.vcf
    ├── AllEditingSites.hg38.txt
    └── PoN.scRNAseq.hg38.tsv

``` 


Then you can run test dataset. 
With this code, conda environment setting is also done.

```bash
snakemake --use-conda --conda-frontend conda --use-singularity --configfile config/config_test.yaml  --cores 12
```




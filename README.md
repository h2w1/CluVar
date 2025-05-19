# CluVar pipeline
![image](https://github.com/user-attachments/assets/5915c4be-ac88-45a0-9a3b-7a3567ef3f0d)


## 1. Setting Environment

Create a conda virtual environment for FastTENET.

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





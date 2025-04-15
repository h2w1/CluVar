# CluVar: Clustering of Variants using autoencoder for inferring the phylogeny of cancer subclones from single cell RNA sequencing data

This repository contains PyTorch code for clustering single-cell RNA sequencing (scRNA-seq) data using CluVar.

## 📁 Project Structure

- `main.py` – Entry point to train model and perform clustering.
- `model.py` – Defines the AutoEncoder model.
- `train.py` – Training and clustering functions.
- `dataset.py` – Loads scRNA-seq data.
- `utils.py` – Utility functions like seed fixing, ARI evaluation, and result saving.
- `search_hyperparams.py` – (Optional) Searches for the best hyperparameter configuration based on loss and ARI.
- `requirements.txt` – Required Python packages.


## 💻 How to Run

1. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

2. Place data files in the `data/` folder.

3. Run:
    ```bash
   ## For Simulation Data, you can use promt below
    python main.py --data_type 'sim'--num_feature 300 --learning_rate 0.0005 --epochs 250 --num_clusters 7 --init_type Xavier
    ```
   ```bash
   ## For Real Data, you can use promt below
    python main.py --data_type 'merged' --num_feature 300 --learning_rate 0.0005 --epochs 250 --num_clusters 7 --init_type Xavier
    ```


Simulation data will output ARI score. Real data saves clustering results without evaluation.

## 🧪 Data Format

- Input `.tsv` with shape (samples, features)
- Optional `.tsv` target file (for simulation only)

## 🔍 Optional: Hyperparameter Search

CluVar supports automatic hyperparameter search using `search_hyperparams.py`.  
This script evaluates different combinations of initialization types, activation functions, and learning rates, and recommends the top 3 configurations based on reconstruction loss.

> ⚠️ Note: This search can take a long time.  
> You can skip it and directly run `main.py` using the recommended parameters shown ubove (`--learning_rate 0.0005`, `--init_type Xavier`, etc.)

### To run the hyperparameter search:

```bash
# For Simulation Data
python search_hyperparams.py --data_type 'sim' --num_clusters 7

# For Real Data
python search_hyperparams.py --data_type 'merged'

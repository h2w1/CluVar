# CluVar: Clustering~제목 반영하기

This repository contains PyTorch code for clustering single-cell RNA sequencing (scRNA-seq) data using CluVar.

## 📁 Project Structure

- `main.py` – Entry point to train model and perform clustering.
- `model.py` – Defines the AutoEncoder model.
- `train.py` – Training and clustering functions.
- `dataset.py` – Loads scRNA-seq data.
- `utils.py` – Utility functions like seed fixing, ARI evaluation, and result saving.
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
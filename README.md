# NvwaCE

This repository contains code for predicting regulatory effects of human genome variants from sequence with NvwaCE models and training new sequence-based regulatory model with any chromatin accessibility profile.

---

## Usage

Take Mouse(Female) as an example to run a standard model training procedures.

1. **Data Preprocessing** (Run `female_mus_pre.R` and `mouseF_5wcells_Generate_Dataset.ipynb` in the folder: `scripts/ResNextATAC/0_dataset_preprocessing/mouseF/`)  
   - Get peak matrix
   - Quality control (filter cells and features)
   - Organize file structure

2. **Training and Evaluating** (Run `bash run.sh` in the folder: `scripts/ResNextATAC/1_train/mouseF/`) 
   - Train and evaluate
   - Save the best model, metrics, etc.


3. **Whole genome predicting**  (Run `bash run.sh` in the folder: `scripts/ResNextATAC/6_predict_whole_genome/human_genome/`)
   - Scan the sliding windows in the human genome using the pre-trained model

---

## Experimental Environment  

- **Operating System:** Ubuntu 20.04.6 LTS (running on WSL2, GNU/Linux 5.15.0-134-generic x86_64) and Slurm cluster
- **Package Management and Environment Control:** Conda  
- **Tools and Versions Used:**  
  - torch==2.4.0+cu121
- **Training Resources:**  NVIDIA A100 Tensor Core GPUs (80GB VRAM) with multi-GPU configurations is recommended for distributed training paradigms.

---

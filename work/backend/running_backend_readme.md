# Pf-HaploAtlas Backend Process

This guide will cover steps to generate haplotype data files for 5,102 genes of the *P. falciparum* core genome. 

## Prerequsites:
You will need to have the following Python packages installed: 

```
numpy==1.21.5
pandas==2.0.3
PyVCF==0.6.8
biopython==1.75
scikit-allel==1.3.6
pyfasta==0.5.2
filelock==3.15.3
```

## Step 1: Run array job to generate new .pkl.xz files

### Step 1.1: Repo set-up
 - If you have not already done so, clone the pf-haploatlas repo to your local machine.
 - If you have a copy already, make sure to `git pull` to update your main branch
 - Create a new branch for the backend file regeneration work: `git checkout -b <branch_name>` and set up tracking with `git push -u origin main <branch_name>`

### Step 1.2: Navigate to correct location
The script is set up to operate from `/work/backend`. So `cd` to here. 

### Step 1.3: Submit job to LSF queue
The code to generate the .pkl.xz files is contained within `generate_gene_summary.py`. We split up the work of generating 5,102 .pkl.xz files across the same number of jobs - this means work can run in parallel and happens faster. 

Submit an array job to the LSF queue using the line of code contained within `generate_gene_summary.py`: 
```
bsub -J "gene_sums[1-5102]%500" -o output.%J%I -e error.%J%I -M6000 -R "select[mem>6000] rusage[mem=6000] span[hosts=1]" -G team342 'python3 generate_gene_summary.py ${LSB_JOBINDEX}'
```

This will trigger the generation of:
- A directory in which `.pkl.xz` files will be sent: `out_pkls`. 
- One `<gene_id>.pkl.xz` per gene: contains the gene-specific haplotype information
- One `.error` and `.output` file per gene: can be used to debug job if needed
- A `gene_logs.json` file: records basic run info per gene
- A `job_logs.json` file: records basic job info

## Step 2: Generate logs and clean up

### Step 2.1: Run job_log_parser.sh
- Cd to `pf-haploatlas/work/backend/`
- Run `./job_log_parser.sh` in the terminal. This will:
    - Archive the old `job_logs.json` in `job_logs_archive/`
    - Generate a new `job_logs.json` to replace the old one. This contains information on the job ID, run date, runtime, etc.
    - Clean up `.error` and `.output` files by removing them
    - Archive a copy of the `generate_gene_summary.py` file just run, dated and stored in `script_archive/`
    - Create a new directory with a date `YYYY-MM-DD_pkl_files` within `app/files`, into which files from `out_pkl` are copied. Files within `out_pkl` are then deleted when a successful copy is confirmed.


## Step 3: Remove old pkl files
Now new pkl.xzs have been generated, we want to remove the older set from the repo in `app/files/` to save on memory.
- Navigate to `app/files`
- Identify the older pkl_files directory by the date prefix
- Remove the files and then the directory

## Step 4: Point Pf-HaploAtlas frontend to new .pkl directory
- Navigate to `pf-haploatlas/app/src/`
- Open `utils.py`
- Update the `base_path` on Line 5 with the new directory name which contains the new .pkl.xz files: e.g. `base_path = "app/files/2024-05-29_pkl_files"``

## Step 5: Add, commit, push changes to repo
Make a PR to submit these changes for review. 
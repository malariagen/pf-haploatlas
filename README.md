# Pf-HaploAtlas
The _Plasmodium falciparum_ Haplotype Atlas (or Pf-HaploAtlas) allows anyone with an internet connection to study and track genetic mutations across any gene in the _P. falciparum_ genome! The app provides visualisations of haplotypes for all 4,952 core genes by using data from 40,074 samples, from 49 countries, and spread between the years 1944 and 2023, facilitating comprehensive spatial and temporal analyses of genes and variants of interest. 

Pf-HaploAtlas currently uses data generated using the MalariaGEN Pf9 whole genome sequencing data release, but will expand with each new MalariaGEN _Plasmodium_ data release. 

We encourage users to access and share the app using the following stable link to prevent outages in service: https://apps.malariagen.net/pf-haploatlas.




# 10-minute demo of the online app (also available in the [app](https://apps.malariagen.net/pf-haploatlas))

[![Pf-HaploAtlas Demo](https://img.youtube.com/vi/48f4r2frcdk/0.jpg)](https://www.youtube.com/watch?v=48f4r2frcdk)





# Pf-HaploAtlas manuscript on Bioinformatics
When publishing work that uses data and/or plots from the Pf-HaploAtlas, please cite the following: 

> Chiyun Lee, Eyyüb S Ünlü, Nina F D White, Jacob Almagro-Garcia, Cristina V Ariani, Richard D Pearson, Pf-HaploAtlas: an interactive web app for spatiotemporal analysis of Plasmodium falciparum genes, Bioinformatics, Volume 40, Issue 11, November 2024, btae673, https://doi.org/10.1093/bioinformatics/btae673





# How to run the app locally
This repository contains a Streamlit app. Follow the steps below to clone the repository and run the app locally.

### Prerequisites
- Python 3.10
- pip (Python package installer)

### 1. Clone the Repository
Open your terminal or command prompt and clone the repository using the following command:

```
git clone git@github.com:malariagen/pf-haploatlas.git
cd pf-haploatlas
```

### 2. Create a virtual environment (optional but recommended)
```
python -m venv .venv
```
Then, activate the virtual environment:
- On Windows: ```.venv\Scripts\activate```
- On macOS/Linux: ```source .venv/bin/activate```

### 3. Install dependencies
```
pip install -r requirements.txt
```

### 4. Run the app
```
python -m streamlit run app/main.py
```
The app should naturally open in your browser but if not, click on the appropriate URL that appears in the terminal. For further details, please refer to the [Streamlit documentation](https://streamlit.io/). 





# Contributing
We strongly encourage you to submit feature requests! Please open an issue [here](https://github.com/malariagen/pf-haploatlas/issues) where you will see a list of features already requested within our community, or contact support@malariagen.net!






# License
This project is licensed under the MIT License. See the LICENSE file for details.

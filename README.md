# scArches_reproducibility


Reproducing results from the [scArches](https://github.com/theislab/scarches).

## Getting Started
In the first step, you will need to download datasets, latent space representations, and pre-trained models to run the each notebook and reproduce the results.

All datasets are available in this drive [directory](https://drive.google.com/drive/folders/1QQXDuUjKG8CTnwWW_u83MDtdrBXr8Kpq?usp=sharing).
All latent space represenstations and pre-trained models are available in this drive [directory](https://drive.google.com/drive/folders/1QqNCBv19-iF9whYIMkVL6DyBRH0x_G57?usp=sharing). After you downloaded necessary datasets and latents, you can run the each notebook and reproduce the results.

## Table of Notebooks 


### Data Analysis
Study       | notebook path     
---------------| ---------------
| [*Toy (Splatter)*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/data_analysis/splatter.ipynb)| notebooks/data_analysis/pancreas.ipynb| 
| [*Pancreas*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/data_analysis/pancreas.ipynb)| notebooks/data_analysis/pancreas.ipynb| 
| [*Mouse Brain*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/data_analysis/MouseBrain.ipynb)|notebooks/data_analysis/MouseBrain.ipynb| 
| [*Tabula Senis Muris*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/data_analysis/TabulaSenis.ipynb)| notebooks/data_analysis/TabulaSenis.ipynb| 
| [*Mouse Cell Atlas*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/data_analysis/MouseCellAtlas.ipynb)| notebooks/data_analysis/MouseCellAtlas.ipynb| 
| [*Panorama*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/data_analysis/panorama.ipynb)| notebooks/data_analysis/panorama.ipynb| 
| [*Covid-19*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/data_analysis/Covid-19.ipynb)| notebooks/data_analysis/Covid-19.ipynb| 

### Notebooks (including paper plots) 
Notebook  | path     
---------------| ---------------
| [*Toy dataset*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/splatter.ipynb)| notebooks/sample_runs/splatter.ipynb| 
| [*Sample Pancreas training with classification*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/train_Pancreas.ipynb)| notebooks/sample_runs/train_Pancreas.ipynb| 
| [*Method Comparison - Pancreas*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/figures/methodComparison-pancreas.ipynb)| notebooks/figures/methodComparison-pancreas.ipynb| 
| [*Method Comparison - Mouse Brain*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/figures/methodComparison-mousebrain.ipynb)| notebooks/figures/methodComparison-mousebrain.ipynb| 
| [*Iterative Surgery - Pancreas (Alpha)*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/IterativeSurgery_with_OutOfSample/OoS+IS_Pancreas_Alpha+Gamma.ipynb)| notebooks/IterativeSurgery_with_OutOfSample/OoS+IS_Pancreas_Alpha.ipynb| 
| [*Iterative Surgery - Pancreas (Alpha + Gamma)*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/IterativeSurgery_with_OutOfSample/OoS+IS_Pancreas_Alpha+Gamma.ipynb)| notebooks/IterativeSurgery_with_OutOfSample/OoS+IS_Pancreas_Alpha+Gamma.ipynb| 
| [*Out of sample batch correction - Pancreas (Alpha + Gamma)*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/OutOfSample/OutOfSample_Pancreas_Alpha.ipynb)| notebooks/OutOfSample/OutOfSample_Pancreas_Alpha.ipynb| 
| [*Out of sample batch correction - Pancreas (Alpha + Gamma)*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/OutOfSample/OutOfSample_Pancreas_Alpha+Gamma.ipynb)| notebooks/OutOfSample/OutOfSample_Pancreas_Alpha+Gamma.ipynb| 
| [*Integration and Classification - TabulaSenisMuris*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/tabula_senis_muris.ipynb)| notebooks/sample_runs/tabula_senis_muris.ipynb| 
| [*Integration and Classification - MouseCellAtlas (with Tabula Senis Muris)*](https://nbviewer.jupyter.org/github/theislab/scnet_reproducibility/blob/master/notebooks/sample_runs/tabula_senis_mca.ipynb)| notebooks/sample_runs/tabula_senis_mca.ipynb| 

To run the notebooks and scripts you need following packages :

scnet, tensorflow, keras, scanpy, numpy, scikit-learn, matplotlib, scipy and splatter(R).

### Benchmarks 
Method  | Notebook Path     
---------------| ---------------
| [*mnnCorrect*](https://nbviewer.jupyter.org/github.com/theislab/scnet_reproducibility/blob/master/notebooks/Benchmarks/mnnCorrect.ipynb)| notebooks/Benchmarks/mnnCorrect.ipynb| 
| [*Conos*](https://nbviewer.jupyter.org/github.com/theislab/scnet_reproducibility/blob/master/notebooks/Benchmarks/Conos.ipynb)| notebooks/Benchmarks/Conos.ipynb| 
| [*Harmony*](https://nbviewer.jupyter.org/github.com/theislab/scnet_reproducibility/blob/master/notebooks/Benchmarks/Harmony.ipynb)| notebooks/Benchmarks/Harmony.ipynb| 
| [*Liger*](https://nbviewer.jupyter.org/github.com/theislab/scnet_reproducibility/blob/master/notebooks/Benchmarks/Liger.ipynb)| notebooks/Benchmarks/Liger.ipynb| 
| [*Scanorama*](https://nbviewer.jupyter.org/github.com/theislab/scnet_reproducibility/blob/master/notebooks/Benchmarks/Scanorama.ipynb)| notebooks/Benchmarks/Conos.ipynb| 
| [*PCA*](https://nbviewer.jupyter.org/github.com/theislab/scnet_reproducibility/blob/master/notebooks/Benchmarks/PCA.ipynb)| notebooks/Benchmarks/PCA.ipynb| 
| [*Seurat*](https://nbviewer.jupyter.org/github.com/theislab/scnet_reproducibility/blob/master/notebooks/Benchmarks/Seurat.ipynb)| notebooks/Benchmarks/Seurat.ipynb| 

In the notebooks, the data is assumed to be in a folder named "data" in the same directory as the notebook.

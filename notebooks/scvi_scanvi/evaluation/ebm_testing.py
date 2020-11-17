from evaluation.utils import entropy_batch_mixing
import scanpy as sc
import os

adata = sc.read(os.path.expanduser(f'~/Documents/git_repos/full_data.h5ad'))
ebm_s = entropy_batch_mixing(adata)
print(ebm_s)
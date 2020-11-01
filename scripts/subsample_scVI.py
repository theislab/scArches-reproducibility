import sys

sys.path.append("../")

import scanpy as sc
from sklearn.preprocessing import LabelEncoder
from _scVI import ADataset, scVI_Trainer
from scvi.models.vae import VAE
from scvi.dataset import AnnDatasetFromAnnData
import numpy as np
import torch
import os
import argparse
from scarches.utils import remove_sparsity

DATASETS = {
    "pancreas": {"name": "pancreas", "batch_key": "study", "cell_type_key": "cell_type",
                 "target": ["Pancreas SS2", "Pancreas CelSeq2"], "HV": False},
    "pbmc": {"name": "pbmc", "batch_key": "study", "cell_type_key": "cell_type", "target": ["Drop-seq", "inDrops"],
             "HV": False},
    "toy": {"name": "toy", "batch_key": "batch", "cell_type_key": "celltype", "target": ["Batch8", "Batch9"],
            "HV": False},
    "mouse_brain": {"name": "mouse_brain", "batch_key": "study", "cell_type_key": "cell_type",
                    "target": ["Tabula_muris", "Zeisel"], "HV": False},
}

parser = argparse.ArgumentParser(description='scNet')
arguments_group = parser.add_argument_group("Parameters")
arguments_group.add_argument('-d', '--data', type=str, required=True,
                             help='data name')
args = vars(parser.parse_args())

data_dict = DATASETS[args['data']]
data_name = data_dict['name']
batch_key = data_dict['batch_key']
cell_type_key = data_dict['cell_type_key']
target_batches = data_dict['target']
highly_variable = data_dict['HV']

adata = sc.read(f"./data/{data_name}/{data_name}_normalized.h5ad")
adata.X = adata.raw.X
adata = remove_sparsity(adata)

n_epochs = 1000
lr = 1e-3
early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
use_batches = True
n_samples = adata.shape[0]

for i in range(5):
    for subsample_frac in [0.1, 0.2, 0.4, 0.6, 0.8, 1.0]:
        final_adata = None
        for target in target_batches:
            adata_sampled = adata[adata.obs[batch_key] == target, :]
            keep_idx = np.loadtxt(f'./data/subsample/{data_name}/{target}/{subsample_frac}/{i}.csv', dtype='int32')
            adata_sampled = adata_sampled[keep_idx, :]

            if final_adata is None:
                final_adata = adata_sampled
            else:
                final_adata = final_adata.concatenate(adata_sampled)

        final_adata.obs['labels'] = LabelEncoder().fit_transform(final_adata.obs[cell_type_key])
        final_adata.obs['batch_indices'] = LabelEncoder().fit_transform(final_adata.obs[batch_key])

        scvi_dataset = AnnDatasetFromAnnData(final_adata)

        print(scvi_dataset.n_batches, scvi_dataset.n_labels)

        vae = VAE(scvi_dataset.nb_genes, n_batch=scvi_dataset.n_batches * use_batches, n_layers=2, n_hidden=128)

        model = scVI_Trainer(vae, scvi_dataset,
                             train_size=0.80,
                             frequency=n_epochs,
                             early_stopping_kwargs=early_stopping_kwargs)
        os.makedirs(f"./results/subsample/{data_name}", exist_ok=True)
        model.train(f"./results/subsample/{data_name}/scVI_frac={subsample_frac}-{i}.csv", n_epochs=n_epochs, lr=lr)

        os.makedirs("./models/scVI/subsample/", exist_ok=True)
        torch.save(model.model.state_dict(), f"./models/scVI/subsample/{data_name}-{subsample_frac}-{i}.pt")

        posterior = model.create_posterior(model.model, scvi_dataset, indices=np.arange(len(scvi_dataset)))
        latent, batch_indices, labels = posterior.sequential().get_latent()

        latent_adata = sc.AnnData(latent)
        latent_adata.obs[cell_type_key] = final_adata.obs[cell_type_key].values
        latent_adata.obs[batch_key] = final_adata.obs[batch_key].values
        
        latent_adata.write_h5ad(f"./results/subsample/{data_name}/scVI/{subsample_frac}/results_adata_{i}.h5ad")

import os
import matplotlib.pyplot as plt
import scanpy as sc
import time
import json
import scarches as sca
import numpy as np

datasets = ['pancreas', 'brain', 'pbmc']
n_epochs_vae = 500
n_epochs_scanvi = 300
early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
early_stopping_kwargs_scanvi = {
    "early_stopping_metric": "accuracy",
    "save_best_state_metric": "accuracy",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

for dataset in datasets:
    dir_path = os.path.expanduser(
        f'~/Documents/scanvi_vanilla_bench/{dataset}/')
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    control_path = f'{dir_path}controlling/'
    if not os.path.exists(control_path):
        os.makedirs(control_path)

    if dataset == 'brain':
        adata_all = sc.read(
            os.path.expanduser(f'~/Documents/benchmarking_datasets/mouse_brain_subsampled_normalized_hvg.h5ad'))
        batch_key = "study"
        cell_type_key = "cell_type"
    elif dataset == 'pancreas':
        adata_all = sc.read(
            os.path.expanduser(f'~/Documents/benchmarking_datasets/pancreas_normalized.h5ad'))
        batch_key = "study"
        cell_type_key = "cell_type"
    elif dataset == 'pbmc':
        adata_all = sc.read(
            os.path.expanduser(f'~/Documents/benchmarking_datasets/Immune_ALL_human_wo_villani_rqr_normalized_hvg.h5ad'))
        batch_key = "condition"
        cell_type_key = "final_annotation"

    adata = adata_all.raw.to_adata()
    sca.dataset.setup_anndata(adata, batch_key=batch_key, labels_key=cell_type_key)
    vae = sca.models.SCANVI(
        adata,
        "Unknown",
        n_layers=2,
        use_cuda=True,
    )
    print("Labelled Indices: ", len(vae._labeled_indices))
    print("Unlabelled Indices: ", len(vae._unlabeled_indices))
    full_time = time.time()
    vae.train(
        n_epochs_unsupervised=n_epochs_vae,
        n_epochs_semisupervised=n_epochs_scanvi,
        unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
        semisupervised_trainer_kwargs=dict(metrics_to_monitor=["elbo", "accuracy"],
                                           early_stopping_kwargs=early_stopping_kwargs_scanvi),
        frequency=1
    )
    full_time = time.time() - full_time

    adata_latent = sc.AnnData(vae.get_latent_representation())
    adata_latent.obs["predictions"] = vae.predict()
    adata_latent.obs["celltype"] = adata.obs[cell_type_key].tolist()
    adata_latent.obs["batch"] = adata.obs[batch_key].tolist()

    full_acc = np.mean(adata_latent.obs["predictions"].tolist() == adata.obs[cell_type_key])

    sc.pp.neighbors(adata_latent)
    sc.tl.leiden(adata_latent)
    sc.tl.umap(adata_latent)
    adata_latent.write_h5ad(filename=f'{dir_path}full_data.h5ad')

    plt.figure()
    sc.pl.umap(
        adata_latent,
        color=["batch", "celltype"],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}umap_full.png', bbox_inches='tight')

    ref_path = f'{dir_path}model/'
    if not os.path.exists(ref_path):
        os.makedirs(ref_path)
    vae.save(ref_path, overwrite=True)

    times = dict()
    times["full_time"] = full_time
    with open(f'{dir_path}results_times.txt', 'w') as filehandle:
        json.dump(times, filehandle)

    accs = dict()
    accs["full_acc"] = full_acc
    with open(f'{dir_path}results_accs.txt', 'w') as filehandle:
        json.dump(accs, filehandle)

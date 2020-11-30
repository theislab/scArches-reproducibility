import os
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import time
import json
import scvi
import numpy as np

sc.set_figure_params(figsize=(4, 4))

dataset = 'brain'
versions = ['vanilla', 'scarches']
deep_injects_main = [True, False]

n_epochs_vae = 500
early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

for version in versions:
    deep_injects = deep_injects_main
    if version == "vanilla":
        deep_injects = [True]
    for deep_inject in deep_injects:
        # Save right dir path
        if version == "scarches" and deep_inject == True:
            deep_label = "/deep_cond"
        elif version == "scarches" and deep_inject == False:
            deep_label = "/first_cond"
        else:
            deep_label = ""
        dir_path = os.path.expanduser(f'~/Documents/benchmarking_results/full_integration_{version}/scvi/{dataset}{deep_label}/')
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        control_path = f'{dir_path}controlling/'
        if not os.path.exists(control_path):
            os.makedirs(control_path)

        if dataset == 'pancreas':
            adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/'
                                                   f'pancreas_normalized.h5ad'))
            batch_key = "study"
            cell_type_key = "cell_type"
        if dataset == 'brain':
            adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/'
                                                   f'mouse_brain_subsampled_normalized_hvg.h5ad'))
            batch_key = "study"
            cell_type_key = "cell_type"
        if dataset == 'immune_all_human_fig6':
            adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/'
                                                   f'Immune_ALL_human_wo_villani_rqr_normalized_hvg.h5ad'))
            batch_key = "condition"
            cell_type_key = "final_annotation"

        adata = adata_all.raw.to_adata()

        scvi.data.setup_anndata(adata, batch_key=batch_key)

        if version == "vanilla":
            vae = scvi.model.SCVI(
                adata,
                n_layers=2,
                use_cuda=True,
            )
        elif version == "scarches":
            vae = scvi.model.SCVI(
                adata,
                n_layers=2,
                use_cuda=True,
                encode_covariates=True,
                deeply_inject_covariates=deep_inject,
                use_layer_norm="both",
                use_batch_norm="none",
                use_observed_lib_size=True
            )

        full_time = time.time()
        vae.train(n_epochs=n_epochs_vae, frequency=1, early_stopping_kwargs=early_stopping_kwargs)
        full_time = time.time() - full_time

        plt.figure()
        plt.plot(vae.trainer.history["elbo_train_set"][2:], label="train")
        plt.plot(vae.trainer.history["elbo_test_set"][2:], label="test")
        plt.title("Negative ELBO over training epochs")
        plt.legend()
        plt.savefig(f'{control_path}reference_elbo.png', bbox_inches='tight')

        adata.obsm["X_scVI"] = vae.get_latent_representation()

        ref_cropped = sc.AnnData(adata.obsm["X_scVI"])
        ref_cropped.obs["celltype"] = adata.obs[cell_type_key].tolist()
        ref_cropped.obs["batch"] = adata.obs[batch_key].tolist()

        sc.pp.neighbors(ref_cropped)
        sc.tl.leiden(ref_cropped)
        sc.tl.umap(ref_cropped)
        ref_cropped.write_h5ad(filename=f'{dir_path}data.h5ad')

        plt.figure()
        sc.pl.umap(
            ref_cropped,
            color=["batch", "celltype"],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_reference.png', bbox_inches='tight')

        torch.save(vae.model.state_dict(), f'{dir_path}model_state_dict')
        path = f'{dir_path}model/'
        if not os.path.exists(path):
            os.makedirs(path)
        vae.save(path, overwrite=True)

        times = dict()
        times["full_time"] = full_time
        with open(f'{dir_path}results_times.txt', 'w') as filehandle:
            json.dump(times, filehandle)

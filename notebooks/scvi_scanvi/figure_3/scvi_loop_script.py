import os
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import time
import json
import scvi
import numpy as np

sc.set_figure_params(figsize=(4, 4))

deep_injects = [True, False]
dataset = 'brain'
test_nrs = [1, 2, 3]

batch_key = "study"
cell_type_key = "cell_type"

n_epochs_vae = 500
n_epochs_surgery = 500
early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
for deep_inject in deep_injects:
    for test_nr in test_nrs:
        print(f"\n\n TEST NR {test_nr} with conditions deep {deep_inject}--------------------------------------------")
        if deep_inject:
            dir_path = os.path.expanduser(
                f'~/Documents/benchmarking_results/figure_3/scvi/brain/test_{test_nr}_deep_cond/')
        else:
            dir_path = os.path.expanduser(
                f'~/Documents/benchmarking_results/figure_3/scvi/brain/test_{test_nr}_first_cond/')
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        control_path = f'{dir_path}controlling/'
        if not os.path.exists(control_path):
            os.makedirs(control_path)
        if dataset == 'pancreas':
            adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/pancreas_normalized.h5ad'))
            if test_nr == 1:
                reference = ['Pancreas inDrop']
                query = ['Pancreas SS2', 'Pancreas CelSeq2', 'Pancreas CelSeq', 'Pancreas Fluidigm C1']
            elif test_nr == 2:
                reference = ['Pancreas inDrop', 'Pancreas SS2']
                query = ['Pancreas CelSeq2', 'Pancreas CelSeq', 'Pancreas Fluidigm C1']
            elif test_nr == 3:
                reference = ['Pancreas inDrop', 'Pancreas SS2', 'Pancreas CelSeq2']
                query = ['Pancreas CelSeq', 'Pancreas Fluidigm C1']
            elif test_nr == 4:
                reference = ['Pancreas inDrop', 'Pancreas SS2', 'Pancreas CelSeq2', 'Pancreas CelSeq']
                query = ['Pancreas Fluidigm C1']

        if dataset == 'brain':
            adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/mouse_brain_subsampled_normalized_hvg.h5ad'))
            if test_nr == 1:
                reference = ['Rosenberg']
                query = ['Saunders', 'Zeisel', 'Tabula_muris']
            elif test_nr == 2:
                reference = ['Rosenberg', 'Saunders']
                query = ['Zeisel', 'Tabula_muris']
            elif test_nr == 3:
                reference = ['Rosenberg', 'Saunders', 'Zeisel']
                query = ['Tabula_muris']

        adata = adata_all.raw.to_adata()
        print(adata)

        ref_ind = np.array([s in reference for s in adata.obs.study])
        query_ind = np.array([s in query for s in adata.obs.study])
        adata_ref = adata[ref_ind].copy()
        adata_query = adata[query_ind].copy()
        print(adata_ref)
        print(adata_query)

        scvi.data.setup_anndata(adata_ref, batch_key=batch_key)

        vae = scvi.model.SCVI(
            adata_ref,
            n_layers=2,
            use_cuda=True,
            encode_covariates=True,
            deeply_inject_covariates=deep_inject,
            use_layer_norm="both",
            use_batch_norm="none",
            use_observed_lib_size=True
        )

        ref_time = time.time()
        vae.train(n_epochs=n_epochs_vae, frequency=1, early_stopping_kwargs=early_stopping_kwargs)
        ref_time = time.time() - ref_time

        plt.plot(vae.trainer.history["elbo_train_set"][2:], label="train")
        plt.plot(vae.trainer.history["elbo_test_set"][2:], label="test")
        plt.title("Negative ELBO over training epochs")
        plt.legend()
        plt.savefig(f'{control_path}reference_elbo.png', bbox_inches='tight')

        adata_ref.obsm["X_scVI"] = vae.get_latent_representation()

        ref_cropped = sc.AnnData(adata_ref.obsm["X_scVI"])
        ref_cropped.obs["celltype"] = adata_ref.obs[cell_type_key].tolist()
        ref_cropped.obs["batch"] = adata_ref.obs[batch_key].tolist()
        sc.pp.neighbors(ref_cropped)
        sc.tl.leiden(ref_cropped)
        sc.tl.umap(ref_cropped)
        ref_cropped.write_h5ad(filename=f'{dir_path}reference_data.h5ad')

        plt.figure()
        sc.pl.umap(
            ref_cropped,
            color=["batch", "celltype"],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_reference.png', bbox_inches='tight')

        torch.save(vae.model.state_dict(), f'{dir_path}reference_model_state_dict')
        ref_path = f'{dir_path}ref_model/'
        if not os.path.exists(ref_path):
            os.makedirs(ref_path)
        vae.save(ref_path, overwrite=True)

        model = scvi.model.SCVI.load_query_data(
            adata_query,
            ref_path,
            use_cuda=True,
            freeze_dropout=True,
            freeze_expression=True,
            freeze_decoder_first_layer=True,
            freeze_batchnorm_encoder=True,
            freeze_batchnorm_decoder=True,
            freeze_classifier=False,
        )

        query_time = time.time()
        model.train(n_epochs=n_epochs_surgery, frequency=1, early_stopping_kwargs=early_stopping_kwargs, weight_decay=0)
        query_time = time.time() - query_time

        plt.figure()
        plt.plot(model.trainer.history["elbo_train_set"][2:], label="train")
        plt.plot(model.trainer.history["elbo_test_set"][2:], label="test")
        plt.title("Negative ELBO over training epochs")
        plt.legend()
        plt.savefig(f'{control_path}surgery_elbo.png', bbox_inches='tight')

        adata_query.obsm["X_scVI"] = model.get_latent_representation()

        q_cropped = sc.AnnData(adata_query.obsm["X_scVI"])
        q_cropped.obs["celltype"] = adata_query.obs[cell_type_key].tolist()
        q_cropped.obs["batch"] = adata_query.obs[batch_key].tolist()
        sc.pp.neighbors(q_cropped)
        sc.tl.leiden(q_cropped)
        sc.tl.umap(q_cropped)
        q_cropped.write_h5ad(filename=f'{dir_path}query_data.h5ad')

        plt.figure()
        sc.pl.umap(
            q_cropped,
            color=["batch", "celltype"],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_query.png', bbox_inches='tight')

        adata_full = adata_ref.concatenate(adata_query)
        adata_full.uns["_scvi"] = adata_query.uns["_scvi"]
        print(adata_full.obs[batch_key].unique().tolist())
        print(adata_full.obs["_scvi_batch"].unique().tolist())
        adata_full.obsm["X_scVI"] = model.get_latent_representation(adata=adata_full)

        f_cropped = sc.AnnData(adata_full.obsm["X_scVI"])
        f_cropped.obs["celltype"] = adata_full.obs[cell_type_key].tolist()
        f_cropped.obs["batch"] = adata_full.obs[batch_key].tolist()
        sc.pp.neighbors(f_cropped)
        sc.tl.leiden(f_cropped)
        sc.tl.umap(f_cropped)
        f_cropped.write_h5ad(filename=f'{dir_path}full_data.h5ad')
        plt.figure()
        sc.pl.umap(
            f_cropped,
            color=["batch", "celltype"],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_full.png', bbox_inches='tight')

        torch.save(model.model.state_dict(), f'{dir_path}surgery_model_state_dict')
        surgery_path = f'{dir_path}surg_model/'
        if not os.path.exists(surgery_path):
            os.makedirs(surgery_path)
        model.save(surgery_path, overwrite=True)

        times = dict()
        times["ref_time"] = ref_time
        times["query_time"] = query_time
        times["full_time"] = ref_time + query_time
        with open(f'{dir_path}results_times.txt', 'w') as filehandle:
            json.dump(times, filehandle)
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
test_nrs = [1,2,3]    # For Brain only 1-3 possible


batch_key = "study"
cell_type_key = "cell_type"

n_epochs_vae = 500
n_epochs_scanvi = 500
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
early_stopping_kwargs_surgery = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
for deep_inject in deep_injects:
    for test_nr in test_nrs:
        print(f"\n\n TEST NR {test_nr} with conditions deep {deep_inject}--------------------------------------------")
        # Save right dir path
        if deep_inject:
            dir_path = os.path.expanduser(f'~/Documents/benchmarking_results/figure_3/scanvi/{dataset}/test_{test_nr}_deep_cond/')
        else:
            dir_path = os.path.expanduser(f'~/Documents/benchmarking_results/figure_3/scanvi/{dataset}/test_{test_nr}_first_cond/')
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

        scvi.data.setup_anndata(adata_ref, batch_key=batch_key, labels_key=cell_type_key)

        vae = scvi.model.SCANVI(
            adata_ref,
            "Unknown",
            n_layers=2,
            use_cuda=True,
            encode_covariates=True,
            deeply_inject_covariates=deep_inject,
            use_layer_norm="both",
            use_batch_norm="none",
            use_observed_lib_size=True,
        )
        print(adata_ref.obs[cell_type_key].unique())
        print(adata_ref.obs["_scvi_labels"].unique())
        print(adata_ref.obs[batch_key].unique())
        print(adata_ref.obs["_scvi_batch"].unique())
        print("Labelled Indices: ", vae._labeled_indices.shape[0])
        print("Unlabelled Indices: ", vae._unlabeled_indices.shape[0])

        ref_time = time.time()
        vae.train(
            n_epochs_unsupervised=n_epochs_vae,
            n_epochs_semisupervised=n_epochs_scanvi,
            unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
            semisupervised_trainer_kwargs=dict(metrics_to_monitor=["elbo", "accuracy"],
                                               early_stopping_kwargs=early_stopping_kwargs_scanvi),
            frequency=1
        )
        ref_time = time.time() - ref_time
        ref_predictions = vae.predict(adata_ref)
        adata_ref.obsm["X_scANVI"] = vae.get_latent_representation()
        adata_ref.obs["predictions"] = vae.predict()
        print("Acc: {}".format(np.mean(ref_predictions == adata_ref.obs[cell_type_key])))

        plt.figure()
        plt.plot(vae.trainer.history['accuracy_full_dataset'][2:], label="ACC")
        plt.title("ACC")
        plt.legend()
        plt.savefig(f'{control_path}reference_acc.png', bbox_inches='tight')

        plt.figure()
        plt.plot(vae.trainer.history['elbo_full_dataset'][2:], label="ELBO")
        plt.title("ELBO")
        plt.legend()
        plt.savefig(f'{control_path}reference_elbo.png', bbox_inches='tight')

        sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
        sc.tl.leiden(adata_ref)
        sc.tl.umap(adata_ref)
        plt.figure()
        sc.pl.umap(
            adata_ref,
            color=[batch_key, cell_type_key],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_reference.png', bbox_inches='tight')

        adata_ref.write_h5ad(filename=f'{dir_path}reference_data.h5ad')
        torch.save(vae.model.state_dict(), f'{dir_path}reference_model_state_dict')
        ref_path = f'{dir_path}ref_model/'
        if not os.path.exists(ref_path):
            os.makedirs(ref_path)
        vae.save(ref_path, overwrite=True)

        model = scvi.model.SCANVI.load_query_data(
            adata_query,
            ref_path,
            use_cuda=True,
            freeze_batchnorm_encoder=True,
            freeze_batchnorm_decoder=True,
            freeze_expression=True,
        )
        print(adata_query.obs[cell_type_key].unique())
        print(adata_query.obs["_scvi_labels"].unique())
        print(adata_query.obs[batch_key].unique())
        print(adata_query.obs["_scvi_batch"].unique())
        model._unlabeled_indices = np.arange(adata_query.n_obs)
        model._labeled_indices = []
        print("Labelled Indices: ", len(model._labeled_indices))
        print("Unlabelled Indices: ", model._unlabeled_indices.shape[0])

        query_time = time.time()
        model.train(
            n_epochs_semisupervised=n_epochs_surgery,
            train_base_model=False,
            semisupervised_trainer_kwargs=dict(metrics_to_monitor=["accuracy", "elbo"],
                                               weight_decay=0,
                                               early_stopping_kwargs=early_stopping_kwargs_surgery
                                               ),
            frequency=1
        )
        query_time = time.time() - query_time

        adata_query.obsm["X_scANVI"] = model.get_latent_representation()
        adata_query.obs["predictions"] = model.predict()
        query_predictions = model.predict()
        print("Acc: {}".format(np.mean(query_predictions == adata_query.obs[cell_type_key])))

        plt.figure()
        plt.plot(model.trainer.history['accuracy_full_dataset'][2:], label="ACC")
        plt.title("ACC")
        plt.legend()
        plt.savefig(f'{control_path}surgery_acc.png', bbox_inches='tight')

        plt.figure()
        plt.plot(model.trainer.history['elbo_full_dataset'][2:], label="ELBO")
        plt.title("ELBO")
        plt.legend()
        plt.savefig(f'{control_path}surgery_elbo.png', bbox_inches='tight')

        sc.pp.neighbors(adata_query, use_rep="X_scANVI")
        sc.tl.leiden(adata_query)
        sc.tl.umap(adata_query)
        plt.figure()
        sc.pl.umap(
            adata_query,
            color=[batch_key, "cell_type"],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_query.png', bbox_inches='tight')
        adata_query.write_h5ad(filename=f'{dir_path}query_data.h5ad')

        adata_full = adata_ref.concatenate(adata_query)
        adata_full.uns["_scvi"] = adata_query.uns["_scvi"]
        print(adata_full.obs[cell_type_key].unique())
        print(adata_full.obs["_scvi_labels"].unique())
        print(adata_full.obs[batch_key].unique())
        print(adata_full.obs["_scvi_batch"].unique())
        adata_full.obsm["X_scANVI"] = model.get_latent_representation(adata=adata_full)

        full_predictions = model.predict(adata_full)
        print("Acc: {}".format(np.mean(full_predictions == adata_full.obs[cell_type_key])))

        sc.pp.neighbors(adata_full, use_rep="X_scANVI")
        sc.tl.leiden(adata_full)
        sc.tl.umap(adata_full)
        plt.figure()
        sc.pl.umap(
            adata_full,
            color=[batch_key, cell_type_key],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_full.png', bbox_inches='tight')

        sc.pl.umap(
            adata_full[adata_ref.n_obs:],
            color=["predictions", cell_type_key],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}umap_full_pred.png', bbox_inches='tight')
        adata_full.write_h5ad(filename=f'{dir_path}full_data.h5ad')
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
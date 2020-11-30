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
n_epochs_surgery = 300
test_nrs = [1,2,3,4]
batch_key = "study"
cell_type_key = "celltype"

n_epochs_vae = 500
n_epochs_scanvi = 500
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

for deep_inject in deep_injects:
    for test_nr in test_nrs:
        # Save right dir path
        if deep_inject:
            dir_path = os.path.expanduser(f'~/Documents/benchmarking_results/figure_4/scanvi/test_{test_nr}_deep_cond/')
        else:
            dir_path = os.path.expanduser(
                f'~/Documents/benchmarking_results/figure_4/scanvi/test_{test_nr}_first_cond/')
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        control_path = f'{dir_path}controlling/'
        if not os.path.exists(control_path):
            os.makedirs(control_path)

        adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/toy_normalized.h5ad'))
        adata = adata_all.raw.to_adata()

        if test_nr == 1:
            reference = ['Batch1']
            query = ['Batch2', 'Batch3', 'Batch4', 'Batch5']
        elif test_nr == 2:
            reference = ['Batch1', 'Batch2']
            query = ['Batch3', 'Batch4', 'Batch5']
        elif test_nr == 3:
            reference = ['Batch1', 'Batch2', 'Batch3']
            query = ['Batch4', 'Batch5']
        elif test_nr == 4:
            reference = ['Batch1', 'Batch2', 'Batch3', 'Batch4']
            query = ['Batch5']

        ref_ind = np.array([s in reference for s in adata.obs[batch_key]])
        query_ind = np.array([s in query for s in adata.obs[batch_key]])
        adata_ref = adata[ref_ind].copy()
        adata_query = adata[query_ind].copy()

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

        print(adata_ref.obs[cell_type_key].unique().tolist())
        print(adata_ref.obs["_scvi_labels"].unique().tolist())
        print(adata_ref.obs[batch_key].unique().tolist())
        print(adata_ref.obs["_scvi_batch"].unique().tolist())
        print("Labelled Indices: ", len(vae._labeled_indices))
        print("Unlabelled Indices: ", len(vae._unlabeled_indices))

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

        ref_cropped = sc.AnnData(adata_ref.obsm["X_scANVI"])
        ref_cropped.obs["celltype"] = adata_ref.obs[cell_type_key].tolist()
        ref_cropped.obs["batch"] = adata_ref.obs[batch_key].tolist()
        ref_cropped.obs["predictions"] = adata_ref.obs["predictions"].tolist()

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

        model = scvi.model.SCANVI.load_query_data(
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

        print(adata_query.obs[cell_type_key].unique().tolist())
        print(adata_query.obs["_scvi_labels"].unique().tolist())
        print(adata_query.obs[batch_key].unique().tolist())
        print(adata_query.obs["_scvi_batch"].unique().tolist())
        model._unlabeled_indices = np.arange(adata_query.n_obs)
        model._labeled_indices = []
        print("Labelled Indices: ", len(model._labeled_indices))
        print("Unlabelled Indices: ", len(model._unlabeled_indices))

        query_time = time.time()
        model.train(
            n_epochs_semisupervised=n_epochs_surgery,
            train_base_model=False,
            semisupervised_trainer_kwargs=dict(metrics_to_monitor=["accuracy"],
                                               weight_decay=0
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

        q_cropped = sc.AnnData(adata_query.obsm["X_scANVI"])
        q_cropped.obs["celltype"] = adata_query.obs[cell_type_key].tolist()
        q_cropped.obs["batch"] = adata_query.obs[batch_key].tolist()
        q_cropped.obs["predictions"] = adata_query.obs["predictions"].tolist()

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
        print(adata_full.obs[cell_type_key].unique().tolist())
        print(adata_full.obs["_scvi_labels"].unique().tolist())
        print(adata_full.obs[batch_key].unique().tolist())
        print(adata_full.obs["_scvi_batch"].unique().tolist())
        adata_full.obsm["X_scANVI"] = model.get_latent_representation(adata=adata_full)

        adata_full.obs["predictions"] = model.predict(adata_full)
        full_predictions = model.predict(adata_full)
        print("Acc: {}".format(np.mean(full_predictions == adata_full.obs[cell_type_key])))

        f_cropped = sc.AnnData(adata_full.obsm["X_scANVI"])
        f_cropped.obs["celltype"] = adata_full.obs[cell_type_key].tolist()
        f_cropped.obs["batch"] = adata_full.obs[batch_key].tolist()
        f_cropped.obs["predictions"] = adata_full.obs["predictions"].tolist()

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

        plt.figure()
        sc.pl.umap(
            f_cropped,
            color=["predictions"],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{control_path}pred_full.png', bbox_inches='tight')

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

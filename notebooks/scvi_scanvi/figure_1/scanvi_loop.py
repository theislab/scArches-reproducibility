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
leave_out_cell_types = ['Pancreas Alpha', 'Pancreas Gamma']

target_batches = ["Pancreas SS2", "Pancreas CelSeq2"]
batch_key = "study"
cell_type_key = "cell_type"

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

for deep_inject in deep_injects:
    # Save right dir path
    if deep_inject:
        dir_path = os.path.expanduser(f'~/Documents/benchmarking_results/figure_1/scanvi/deep_cond_ood2/')
    else:
        dir_path = os.path.expanduser(f'~/Documents/benchmarking_results/figure_1/scanvi/first_cond_ood2/')

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    control_path = f'{dir_path}controlling/'
    if not os.path.exists(control_path):
        os.makedirs(control_path)
    adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/pancreas_normalized.h5ad'))
    adata = adata_all.raw.to_adata()

    query = np.array([s in target_batches for s in adata.obs[batch_key]])
    query_1 = np.array([s in [target_batches[0]] for s in adata.obs[batch_key]])
    query_2 = np.array([s in [target_batches[1]] for s in adata.obs[batch_key]])
    adata_ref_full = adata[~query].copy()
    adata_ref = adata_ref_full[~adata_ref_full.obs[cell_type_key].isin(leave_out_cell_types)].copy()
    adata_query_1 = adata[query_1].copy()
    adata_query_2 = adata[query_2].copy()

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

    adata_ref.obsm["X_scANVI"] = vae.get_latent_representation()
    adata_ref.obs["predictions"] = vae.predict()

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

    model_1 = scvi.model.SCANVI.load_query_data(
        adata_query_1,
        ref_path,
        use_cuda=True,
        freeze_dropout=True,
        freeze_expression=True,
        freeze_decoder_first_layer=True,
        freeze_batchnorm_encoder=True,
        freeze_batchnorm_decoder=True,
        freeze_classifier=False,
    )
    model_1._unlabeled_indices = np.arange(adata_query_1.n_obs)
    model_1._labeled_indices = []
    query_1_time = time.time()
    model_1.train(
        n_epochs_semisupervised=n_epochs_surgery,
        train_base_model=False,
        semisupervised_trainer_kwargs=dict(metrics_to_monitor=["accuracy"],
                                           weight_decay=0
                                           ),
        frequency=1
    )
    query_1_time = time.time() - query_1_time

    adata_query_1.obsm["X_scANVI"] = model_1.get_latent_representation()
    adata_query_1.obs["predictions"] = model_1.predict()

    plt.figure()
    plt.plot(model_1.trainer.history['accuracy_full_dataset'][2:], label="ACC")
    plt.title("ACC")
    plt.legend()
    plt.savefig(f'{control_path}surgery_acc.png', bbox_inches='tight')

    q1_cropped = sc.AnnData(adata_query_1.obsm["X_scANVI"])
    q1_cropped.obs["celltype"] = adata_query_1.obs[cell_type_key].tolist()
    q1_cropped.obs["batch"] = adata_query_1.obs[batch_key].tolist()
    q1_cropped.obs["predictions"] = adata_query_1.obs["predictions"].tolist()

    sc.pp.neighbors(q1_cropped)
    sc.tl.leiden(q1_cropped)
    sc.tl.umap(q1_cropped)
    q1_cropped.write_h5ad(filename=f'{dir_path}query_1_data.h5ad')

    plt.figure()
    sc.pl.umap(
        q1_cropped,
        color=["batch", "celltype"],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}umap_query_1.png', bbox_inches='tight')

    adata_full_1 = adata_ref.concatenate(adata_query_1)
    adata_full_1.uns["_scvi"] = adata_query_1.uns["_scvi"]
    adata_full_1.obsm["X_scANVI"] = model_1.get_latent_representation(adata=adata_full_1)
    adata_full_1.obs["predictions"] = model_1.predict(adata_full_1)

    f1_cropped = sc.AnnData(adata_full_1.obsm["X_scANVI"])
    f1_cropped.obs["celltype"] = adata_full_1.obs[cell_type_key].tolist()
    f1_cropped.obs["batch"] = adata_full_1.obs[batch_key].tolist()
    f1_cropped.obs["predictions"] = adata_full_1.obs["predictions"].tolist()

    sc.pp.neighbors(f1_cropped)
    sc.tl.leiden(f1_cropped)
    sc.tl.umap(f1_cropped)
    f1_cropped.write_h5ad(filename=f'{dir_path}full_1_data.h5ad')

    plt.figure()
    sc.pl.umap(
        f1_cropped,
        color=["batch", "celltype"],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}umap_full_1.png', bbox_inches='tight')
    sc.pl.umap(
        f1_cropped,
        color=["predictions"],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}pred_full_1.png', bbox_inches='tight')

    torch.save(model_1.model.state_dict(), f'{dir_path}surgery_1_model_state_dict')
    surgery_1_path = f'{dir_path}surg_1_model/'
    if not os.path.exists(surgery_1_path):
        os.makedirs(surgery_1_path)
    model_1.save(surgery_1_path, overwrite=True)

    model_2 = scvi.model.SCANVI.load_query_data(
        adata_query_2,
        surgery_1_path,
        use_cuda=True,
        freeze_dropout=True,
        freeze_expression=True,
        freeze_decoder_first_layer=True,
        freeze_batchnorm_encoder=True,
        freeze_batchnorm_decoder=True,
        freeze_classifier=False,
    )

    model_2._unlabeled_indices = np.arange(adata_query_2.n_obs)
    model_2._labeled_indices = []

    query_2_time = time.time()
    model_2.train(
        n_epochs_semisupervised=n_epochs_surgery,
        train_base_model=False,
        semisupervised_trainer_kwargs=dict(metrics_to_monitor=["accuracy"],
                                           weight_decay=0),
        frequency=1
    )
    query_2_time = time.time() - query_2_time

    plt.figure()
    plt.plot(model_2.trainer.history['accuracy_full_dataset'][2:], label="ACC")
    plt.title("ACC")
    plt.legend()
    plt.savefig(f'{control_path}surgery_2_acc.png', bbox_inches='tight')

    adata_query_2.obsm["X_scANVI"] = model_2.get_latent_representation()
    adata_query_2.obs["predictions"] = model_2.predict()

    q2_cropped = sc.AnnData(adata_query_2.obsm["X_scANVI"])
    q2_cropped.obs["celltype"] = adata_query_2.obs[cell_type_key].tolist()
    q2_cropped.obs["batch"] = adata_query_2.obs[batch_key].tolist()
    q2_cropped.obs["predictions"] = adata_query_2.obs["predictions"].tolist()

    sc.pp.neighbors(q2_cropped)
    sc.tl.leiden(q2_cropped)
    sc.tl.umap(q2_cropped)
    q2_cropped.write_h5ad(filename=f'{dir_path}query_2_data.h5ad')

    plt.figure()
    sc.pl.umap(
        q2_cropped,
        color=["batch", "celltype"],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}umap_query_2.png', bbox_inches='tight')

    adata_full_2 = adata_full_1.concatenate(adata_query_2)
    adata_full_2.uns["_scvi"] = adata_query_2.uns["_scvi"]

    adata_full_2.obsm["X_scANVI"] = model_2.get_latent_representation(adata=adata_full_2)
    adata_full_2.obs["predicitions"] = model_2.predict(adata_full_2)

    f2_cropped = sc.AnnData(adata_full_2.obsm["X_scANVI"])
    f2_cropped.obs["celltype"] = adata_full_2.obs[cell_type_key].tolist()
    f2_cropped.obs["batch"] = adata_full_2.obs[batch_key].tolist()
    f2_cropped.obs["predictions"] = adata_full_2.obs["predictions"].tolist()
    sc.pp.neighbors(f2_cropped)
    sc.tl.leiden(f2_cropped)
    sc.tl.umap(f2_cropped)
    f2_cropped.write_h5ad(filename=f'{dir_path}full_2_data.h5ad')

    plt.figure()
    sc.pl.umap(
        f2_cropped,
        color=["batch", "celltype"],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}umap_full_2.png', bbox_inches='tight')
    plt.figure()
    sc.pl.umap(
        f2_cropped,
        color=["predictions"],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}pred_full_2.png', bbox_inches='tight')

    torch.save(model_2.model.state_dict(), f'{dir_path}surgery_2_model_state_dict')
    surgery_2_path = f'{dir_path}surg_2_model/'
    if not os.path.exists(surgery_2_path):
        os.makedirs(surgery_2_path)
    model_2.save(surgery_2_path, overwrite=True)
    times = dict()
    times["ref_time"] = ref_time
    times["query_1_time"] = query_1_time
    times["query_2_time"] = query_2_time
    times["full_time"] = ref_time + query_1_time + query_2_time
    with open(f'{dir_path}results_times.txt', 'w') as filehandle:
        json.dump(times, filehandle)
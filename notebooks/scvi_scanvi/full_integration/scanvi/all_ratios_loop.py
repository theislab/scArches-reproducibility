import os
import matplotlib.pyplot as plt
import scanpy as sc
import torch
import time
import json
import scvi
import numpy as np

sc.set_figure_params(figsize=(4, 4))

dataset = 'immune_all_human_fig6'
test_nrs = [1,2,3,4]

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
for test_nr in test_nrs:
    print(f"\n\n\nSTART WITH RATIO {test_nr}")
    # Save right dir path
    dir_path = os.path.expanduser(
        f'~/Documents/benchmarking_results/full_integration/scanvi/{dataset}/label_ratio_{test_nr}/')
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
        if test_nr == 1:
            labeled = ['Pancreas inDrop']
            unlabeled = ['Pancreas SS2', 'Pancreas CelSeq2', 'Pancreas CelSeq', 'Pancreas Fluidigm C1']
        elif test_nr == 2:
            labeled = ['Pancreas inDrop', 'Pancreas SS2']
            unlabeled = ['Pancreas CelSeq2', 'Pancreas CelSeq', 'Pancreas Fluidigm C1']
        elif test_nr == 3:
            labeled = ['Pancreas inDrop', 'Pancreas SS2', 'Pancreas CelSeq2']
            unlabeled = ['Pancreas CelSeq', 'Pancreas Fluidigm C1']
        elif test_nr == 4:
            labeled = ['Pancreas inDrop', 'Pancreas SS2', 'Pancreas CelSeq2', 'Pancreas CelSeq']
            unlabeled = ['Pancreas Fluidigm C1']
        elif test_nr == 5:
            labeled = ['Pancreas inDrop', 'Pancreas SS2', 'Pancreas CelSeq2', 'Pancreas CelSeq', 'Pancreas Fluidigm C1']
            unlabeled = []
    if dataset == 'brain':
        adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/'
                                               f'mouse_brain_subsampled_normalized_hvg.h5ad'))
        batch_key = "study"
        cell_type_key = "cell_type"
        if test_nr == 1:
            labeled = ['Rosenberg']
            unlabeled = ['Saunders', 'Zeisel', 'Tabula_muris']
        elif test_nr == 2:
            labeled = ['Rosenberg', 'Saunders']
            unlabeled = ['Zeisel', 'Tabula_muris']
        elif test_nr == 3:
            labeled = ['Rosenberg', 'Saunders', 'Zeisel']
            unlabeled = ['Tabula_muris']
        elif test_nr == 4:
            labeled = ['Rosenberg', 'Saunders', 'Zeisel', 'Tabula_muris']
            unlabeled = []
    if dataset == 'immune_all_human_fig6':
        adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/'
                                               f'Immune_ALL_human_wo_villani_rqr_normalized_hvg.h5ad'))
        batch_key = "condition"
        cell_type_key = "final_annotation"
        if test_nr == 1:
            labeled = ['10X']
            unlabeled = ['Oetjen', 'Sun', 'Freytag']
        elif test_nr == 2:
            labeled = ['10X', 'Oetjen']
            unlabeled = ['Sun', 'Freytag']
        elif test_nr == 3:
            labeled = ['10X', 'Oetjen', 'Sun']
            unlabeled = ['Freytag']
        elif test_nr == 4:
            labeled = ['10X', 'Oetjen', 'Sun', 'Freytag']
            unlabeled = []

    adata = adata_all.raw.to_adata()
    scvi.data.setup_anndata(adata, batch_key=batch_key, labels_key=cell_type_key)

    vae = scvi.model.SCANVI(
        adata,
        "Unknown",
        n_layers=2,
        use_cuda=True,
    )

    labeled_ind = np.array([s in labeled for s in adata.obs.study])
    unlabeled_ind = np.array([s in unlabeled for s in adata.obs.study])
    vae._unlabeled_indices = np.arange(adata.n_obs)[unlabeled_ind]
    vae._labeled_indices = np.arange(adata.n_obs)[labeled_ind]
    print("Labeled Conditions:", adata.obs[batch_key][labeled_ind].unique().tolist())
    print("Unlabeled Conditions:", adata.obs[batch_key][unlabeled_ind].unique().tolist())
    print("Unlabeled Indices: ", vae._unlabeled_indices.shape[0])
    print("Labeled Indices: ", vae._labeled_indices.shape[0])

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

    ref_predictions = vae.predict(adata)
    adata.obsm["X_scANVI"] = vae.get_latent_representation()
    adata.obs["predictions"] = vae.predict()
    print("Acc: {}".format(np.mean(ref_predictions == adata.obs[cell_type_key])))

    plt.figure()
    plt.plot(vae.trainer.history['accuracy_full_dataset'][2:], label="ACC")
    plt.title("ACC")
    plt.legend()
    plt.savefig(f'{control_path}acc.png', bbox_inches='tight')

    plt.figure()
    plt.plot(vae.trainer.history['elbo_full_dataset'][2:], label="ELBO")
    plt.title("ELBO")
    plt.legend()
    plt.savefig(f'{control_path}elbo.png', bbox_inches='tight')

    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.leiden(adata)
    sc.tl.umap(adata)
    plt.figure()
    sc.pl.umap(
        adata,
        color=[batch_key, cell_type_key],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}umap.png', bbox_inches='tight')

    sc.pl.umap(
        adata,
        color=["predictions", cell_type_key],
        frameon=False,
        ncols=1,
        show=False
    )
    plt.savefig(f'{control_path}umap_pred.png', bbox_inches='tight')

    adata.write_h5ad(filename=f'{dir_path}data.h5ad')
    torch.save(vae.model.state_dict(), f'{dir_path}model_state_dict')
    ref_path = f'{dir_path}model/'
    if not os.path.exists(ref_path):
        os.makedirs(ref_path)
    vae.save(ref_path, overwrite=True)

    times = dict()
    times["full_time"] = full_time
    with open(f'{dir_path}results_times.txt', 'w') as filehandle:
        json.dump(times, filehandle)
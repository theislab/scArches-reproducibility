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
surgery_options = ['freezed_expr', 'freezed', 'unfreezed']

batch_key = "study"
cell_type_key = "cell_type"
reference = ['Rosenberg', 'Saunders']
query = ['Zeisel', 'Tabula_muris']

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
adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/mouse_brain_subsampled_normalized_hvg.h5ad'))
adata = adata_all.raw.to_adata()
print(adata)

for deep_inject in deep_injects:
    if deep_inject:
        dir_path = os.path.expanduser(
            f'~/Documents/benchmarking_results/figure_2/scvi/deep_cond/')
    else:
        dir_path = os.path.expanduser(
            f'~/Documents/benchmarking_results/figure_2/scvi/first_cond/')

    ref_model_path = f'{dir_path}reference/'
    if not os.path.exists(ref_model_path):
        os.makedirs(ref_model_path)
    ref_control_path = f'{ref_model_path}controlling/'
    if not os.path.exists(ref_control_path):
        os.makedirs(ref_control_path)

    ref_ind = np.array([s in reference for s in adata.obs.study])
    adata_ref = adata[ref_ind].copy()
    scvi.data.setup_anndata(adata_ref, batch_key=batch_key)
    print(adata_ref)

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

    plt.figure()
    plt.plot(vae.trainer.history["elbo_train_set"][2:], label="train")
    plt.plot(vae.trainer.history["elbo_test_set"][2:], label="test")
    plt.title("Negative ELBO over training epochs")
    plt.legend()
    plt.savefig(f'{ref_control_path}reference_elbo.png', bbox_inches='tight')

    adata_ref.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata_ref, use_rep="X_scVI")
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
    plt.savefig(f'{ref_control_path}umap_reference.png', bbox_inches='tight')

    adata_ref.write_h5ad(filename=f'{ref_model_path}reference_data.h5ad')
    torch.save(vae.model.state_dict(), f'{ref_model_path}reference_model_state_dict')
    ref_path = f'{ref_model_path}ref_model/'
    if not os.path.exists(ref_path):
        os.makedirs(ref_path)
    vae.save(ref_path, overwrite=True)

    for surgery_option in surgery_options:
        print(f"\n\n Surgery {surgery_option} with conditions deep {deep_inject}---------------------------------------")
        surg_model_path = f'{dir_path}{surgery_option}/'
        if not os.path.exists(surg_model_path):
            os.makedirs(surg_model_path)
        surg_control_path = f'{surg_model_path}controlling/'
        if not os.path.exists(surg_control_path):
            os.makedirs(surg_control_path)

        if surgery_option == 'freezed_expr':
            freeze_exp = True
            full_retrain = False
        if surgery_option == 'freezed':
            freeze_exp = False
            full_retrain = False
        if surgery_option == 'unfreezed':
            freeze_exp = False
            full_retrain = True

        query_ind = np.array([s in query for s in adata.obs.study])
        adata_query = adata[query_ind].copy()
        print(adata_query)

        model = scvi.model.SCVI.load_query_data(
            adata_query,
            ref_path,
            use_cuda=True,
            freeze_expression=freeze_exp,
            full_retrain=full_retrain
        )
        for name, p in model.model.named_parameters():
            if ("z_encoder.encoder.fc_layers" in name or "decoder.px_decoder.fc_layers" in name) and "weight" in name:
                print(name)
                print(p[0])

        query_time = time.time()
        model.train(n_epochs=n_epochs_surgery, frequency=1, early_stopping_kwargs=early_stopping_kwargs, weight_decay=0)
        query_time = time.time() - query_time

        plt.figure()
        plt.plot(model.trainer.history["elbo_train_set"][2:], label="train")
        plt.plot(model.trainer.history["elbo_test_set"][2:], label="test")
        plt.title("Negative ELBO over training epochs")
        plt.legend()
        plt.savefig(f'{surg_control_path}surgery_elbo.png', bbox_inches='tight')

        adata_query.obsm["X_scVI"] = model.get_latent_representation()
        sc.pp.neighbors(adata_query, use_rep="X_scVI")
        sc.tl.leiden(adata_query)
        sc.tl.umap(adata_query)
        plt.figure()
        sc.pl.umap(
            adata_query,
            color=[batch_key, cell_type_key],
            frameon=False,
            ncols=1,
            show=False
        )
        plt.savefig(f'{surg_control_path}umap_query.png', bbox_inches='tight')
        adata_query.write_h5ad(filename=f'{surg_model_path}query_data.h5ad')

        adata_full = adata_ref.concatenate(adata_query)
        adata_full.uns["_scvi"] = adata_query.uns["_scvi"]
        print(adata_full.obs[batch_key].unique())
        print(adata_full.obs["_scvi_batch"].unique())
        adata_full.obsm["X_scVI"] = model.get_latent_representation(adata=adata_full)

        sc.pp.neighbors(adata_full, use_rep="X_scVI")
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
        plt.savefig(f'{surg_control_path}umap_full.png', bbox_inches='tight')

        adata_full.write_h5ad(filename=f'{surg_model_path}full_data.h5ad')
        torch.save(model.model.state_dict(), f'{surg_model_path}surgery_model_state_dict')
        surgery_path = f'{surg_model_path}surg_model/'
        if not os.path.exists(surgery_path):
            os.makedirs(surgery_path)
        model.save(surgery_path, overwrite=True)

        times = dict()
        times["ref_time"] = ref_time
        times["query_time"] = query_time
        times["full_time"] = ref_time + query_time
        with open(f'{surg_model_path}results_times.txt', 'w') as filehandle:
            json.dump(times, filehandle)

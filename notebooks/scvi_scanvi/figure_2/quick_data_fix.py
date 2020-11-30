import os
import scvi
import numpy as np
import scanpy as sc
import torch
import matplotlib.pyplot as plt

surgery_option = 'freezed_expr'
#surgery_option = 'freezed'
#surgery_option = 'unfreezed'

batch_key = "study"
cell_type_key = "cell_type"
reference = ['Rosenberg', 'Saunders']
query = ['Zeisel', 'Tabula_muris']
adata_all = sc.read(os.path.expanduser(f'~/Documents/benchmarking_datasets/mouse_brain_subsampled_normalized_hvg.h5ad'))
adata = adata_all.raw.to_adata()
ref_ind = np.array([s in reference for s in adata.obs.study])
adata_ref = adata[ref_ind].copy()
scvi.data.setup_anndata(adata_ref, batch_key=batch_key)
dir_path = os.path.expanduser(
            f'~/Documents/benchmarking_results/figure_2/scvi/first_cond/')
ref_model_path = f'{dir_path}reference/'
ref_path = f'{ref_model_path}ref_model/'
query_ind = np.array([s in query for s in adata.obs.study])
adata_query = adata[query_ind].copy()
print(adata_query)
model = scvi.model.SCVI.load_query_data(
            adata_query,
            ref_path,
            use_cuda=True,
            unfrozen=False,
            freeze_expression=True,
            freeze_decoder_first_layer=True,
            freeze_batchnorm_decoder=True,
            freeze_batchnorm_encoder=True,
            freeze_dropout=True,
            freeze_classifier=False,
        )

surg_model_path = f'{dir_path}{surgery_option}/'
if not os.path.exists(surg_model_path):
    os.makedirs(surg_model_path)
surg_control_path = f'{surg_model_path}controlling/'
if not os.path.exists(surg_control_path):
    os.makedirs(surg_control_path)

model.model.load_state_dict(torch.load(f'{surg_model_path}surgery_model_state_dict'), strict=True)
model.is_trained_ = True

adata_query.obsm["X_scVI"] = model.get_latent_representation()

q_cropped = sc.AnnData(adata_query.obsm["X_scVI"])
q_cropped.obs["celltype"] = adata_query.obs[cell_type_key].tolist()
q_cropped.obs["batch"] = adata_query.obs[batch_key].tolist()
sc.pp.neighbors(q_cropped)
sc.tl.leiden(q_cropped)
sc.tl.umap(q_cropped)
q_cropped.write_h5ad(filename=f'{surg_model_path}query_data.h5ad')
plt.figure()
sc.pl.umap(
    q_cropped,
    color=["batch", "celltype"],
    frameon=False,
    ncols=1,
    show=False
)
plt.savefig(f'{surg_control_path}umap_query_check.png', bbox_inches='tight')

adata_full = adata_ref.concatenate(adata_query)
adata_full.uns["_scvi"] = adata_query.uns["_scvi"]
print(adata_full.obs[batch_key].unique())
print(adata_full.obs["_scvi_batch"].unique())
adata_full.obsm["X_scVI"] = model.get_latent_representation(adata=adata_full)

f_cropped = sc.AnnData(adata_full.obsm["X_scVI"])
f_cropped.obs["celltype"] = adata_full.obs[cell_type_key].tolist()
f_cropped.obs["batch"] = adata_full.obs[batch_key].tolist()
sc.pp.neighbors(f_cropped)
sc.tl.leiden(f_cropped)
sc.tl.umap(f_cropped)
f_cropped.write_h5ad(filename=f'{surg_model_path}full_data.h5ad')
plt.figure()
sc.pl.umap(
    f_cropped,
    color=["batch", "celltype"],
    frameon=False,
    ncols=1,
    show=False
)
plt.savefig(f'{surg_control_path}umap_full_check.png', bbox_inches='tight')
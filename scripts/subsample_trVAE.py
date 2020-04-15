import argparse
import os

import numpy as np
import scanpy as sc

import scnet as sn

DATASETS = {
    "pancreas": {"name": "pancreas", "batch_key": "study", "cell_type_key": "cell_type",
                 "target": ["Pancreas SS2", "Pancreas CelSeq2"]},
    "toy": {"name": "toy", "batch_key": "batch", "cell_type_key": "celltype", "target": ["Batch8", "Batch9"]},
    "pbmc": {"name": "pbmc_subset", "batch_key": "study", "cell_type_key": "cell_type",
             "target": ["inDrops", "Drop-seq"]},
    "mouse_brain": {"name": "mouse_brain", "batch_key": "study", "cell_type_key": "cell_type",
                    "target": ["Tabula_muris", "Zeisel"]},
}


def train_and_evaluate(data_dict, loss_fn="mse"):
    data_name = data_dict['name']
    cell_type_key = data_dict['cell_type_key']
    condition_key = data_dict['batch_key']
    target_conditions = data_dict['target']

    adata = sc.read(f"./data/{data_name}/{data_name}_normalized.h5ad")

    path_to_save = f"./results/subsample/{data_name}/"
    os.makedirs(path_to_save, exist_ok=True)

    if loss_fn == "mse":
        clip_value = 1000.0
    else:
        clip_value = 3.0

    for i in range(5):
        scores = []
        for subsample_frac in [1.0, 0.8, 0.6, 0.4, 0.2, 0.1]:
            final_adata, raw_out_of_sample = None, None
            for condition in target_conditions:
                condition_adata = adata[adata.obs[condition_key] == condition]
                keep_idx = np.loadtxt(f'./data/subsample/{data_name}/{condition}/{subsample_frac}/{i}.csv',
                                      dtype='int32')
                condition_adata_subsampled = condition_adata[keep_idx, :]
                final_adata = condition_adata_subsampled if final_adata is None \
                    else final_adata.concatenate(condition_adata_subsampled)
                raw_out_of_sample = sc.AnnData(
                    condition_adata_subsampled.raw.X) if raw_out_of_sample is None else raw_out_of_sample.concatenate(
                    sc.AnnData(condition_adata_subsampled.raw.X))
            final_adata.raw = raw_out_of_sample

            train_adata, valid_adata = sn.utils.train_test_split(final_adata, 0.80)
            n_conditions = len(train_adata.obs[condition_key].unique().tolist())

            z_dim = 10
            architecture = [128, 64, 32]

            network = sn.archs.scNet(x_dimension=train_adata.shape[1],
                                     z_dimension=z_dim,
                                     architecture=architecture,
                                     use_batchnorm=False,
                                     n_conditions=n_conditions,
                                     lr=0.001,
                                     alpha=0.00005,
                                     beta=1000.0,
                                     eta=1.0,
                                     clip_value=clip_value,
                                     loss_fn=loss_fn,
                                     model_path=f"./models/CVAE/subsample/{data_name}/scratch/",
                                     dropout_rate=0.05,
                                     output_activation='relu')

            conditions = final_adata.obs[condition_key].unique().tolist()
            condition_encoder = sn.utils.create_dictionary(conditions, [])

            network.train(train_adata,
                          valid_adata,
                          condition_key=condition_key,
                          cell_type_key=cell_type_key,
                          le=condition_encoder,
                          n_epochs=10000,
                          batch_size=1024,
                          early_stop_limit=50,
                          lr_reducer=40,
                          n_per_epoch=-1,
                          save=True,
                          retrain=False,
                          verbose=5)

            encoder_labels, _ = sn.utils.label_encoder(
                final_adata, label_encoder=network.condition_encoder, condition_key=condition_key)

            latent_adata = network.to_mmd_layer(final_adata, encoder_labels, encoder_labels)

            latent_adata.write_h5ad(os.path.join(path_to_save, f'trVAE/{subsample_frac}/results_adata_{i}.h5ad'))

            encoder_labels, _ = sn.utils.label_encoder(final_adata, label_encoder=network.condition_encoder,
                                                          condition_key=condition_key)
            latent_adata = network.to_latent(final_adata, encoder_labels)

            asw = sn.metrics.asw(latent_adata, label_key=condition_key)
            ari = sn.metrics.ari(latent_adata, label_key=cell_type_key)
            nmi = sn.metrics.nmi(latent_adata, label_key=cell_type_key)
            knn_15 = sn.metrics.knn_purity(latent_adata, label_key=cell_type_key, n_neighbors=15)
            knn_25 = sn.metrics.knn_purity(latent_adata, label_key=cell_type_key, n_neighbors=25)
            knn_50 = sn.metrics.knn_purity(latent_adata, label_key=cell_type_key, n_neighbors=50)
            knn_100 = sn.metrics.knn_purity(latent_adata, label_key=cell_type_key, n_neighbors=100)
            knn_200 = sn.metrics.knn_purity(latent_adata, label_key=cell_type_key, n_neighbors=200)
            knn_300 = sn.metrics.knn_purity(latent_adata, label_key=cell_type_key, n_neighbors=300)
            ebm_15 = sn.metrics.entropy_batch_mixing(latent_adata, label_key=condition_key, n_pools=1,
                                                        n_neighbors=15)
            ebm_25 = sn.metrics.entropy_batch_mixing(latent_adata, label_key=condition_key, n_pools=1,
                                                        n_neighbors=25)
            ebm_50 = sn.metrics.entropy_batch_mixing(latent_adata, label_key=condition_key, n_pools=1,
                                                        n_neighbors=50)
            ebm_100 = sn.metrics.entropy_batch_mixing(latent_adata, label_key=condition_key, n_pools=1,
                                                         n_neighbors=100)
            ebm_200 = sn.metrics.entropy_batch_mixing(latent_adata, label_key=condition_key, n_pools=1,
                                                         n_neighbors=200)
            ebm_300 = sn.metrics.entropy_batch_mixing(latent_adata, label_key=condition_key, n_pools=1,
                                                         n_neighbors=300)

            scores.append(
                [subsample_frac, asw, ari, nmi, knn_15, knn_25, knn_50, knn_100, knn_200, knn_300, ebm_15, ebm_25,
                 ebm_50, ebm_100, ebm_200, ebm_300])
            print([subsample_frac, asw, ari, nmi, knn_15, knn_25, knn_50, knn_100, knn_200, knn_300, ebm_15, ebm_25,
                   ebm_50, ebm_100, ebm_200, ebm_300])

        scores = np.array(scores)

        filename = "scores_scratch_scNet"
        filename += "_count" if loss_fn == 'nb' else "_normalized"
        filename += f"_{i}.log"

        np.savetxt(os.path.join(path_to_save, filename), X=scores, delimiter=",")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='scNet')
    arguments_group = parser.add_argument_group("Parameters")
    arguments_group.add_argument('-d', '--data', type=str, required=True,
                                 help='data name')
    arguments_group.add_argument('-c', '--loss', type=int, default=0, required=False,
                                 help='if 1 will use nb else will use mse')
    args = vars(parser.parse_args())

    data_name = args['data']
    loss_fn = "nb" if args['loss'] > 0 else "mse"
    data_dict = DATASETS[data_name]
    train_and_evaluate(data_dict=data_dict, loss_fn=loss_fn)

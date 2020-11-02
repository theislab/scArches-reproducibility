import os
import time
import scarches as sca
import scanpy as sc


datasets = {
    'pancreas': {
        'path': '/home/mohsen/data/pancreas/pancreas_normalized.h5ad',
        'condition_key': 'study',
        'cell_type_key': 'cell_type',
        'runs': [
            {'reference': ['Pancreas inDrop']},
            {'reference': ['Pancreas inDrop', 'Pancreas SS2']},
            {'reference': ['Pancreas inDrop', 'Pancreas SS2', 'Pancreas CelSeq2']},
            {'reference': ['Pancreas inDrop', 'Pancreas SS2', 'Pancreas CelSeq2', 'Pancreas CelSeq']},
        ]
    },
    'mouse_brain': {
        'path': '/home/mohsen/data/mouse_brain/mouse_brain_subsampled_normalized_hvg.h5ad',
        'condition_key': 'study',
        'cell_type_key': 'cell_type',
        'runs': [
            {'reference': ['Rosenberg']},
            {'reference': ['Rosenberg', 'Saunders']},
            {'reference': ['Rosenberg', 'Saunders', 'Zeisel']},
        ]
    },
    'pbmc': {
        'path': '/home/mohsen/data/PBMC/Immune_ALL_human_wo_villani_rqr_normalized_hvg.h5ad',
        'condition_key': 'condition',
        'cell_type_key': 'final_annotation',
        'runs': [
            {'reference': ['10X']},
            {'reference': ['10X', 'Oetjen_A', 'Oetjen_P', 'Oetjen_U']},
            {'reference': ['10X', 'Oetjen_A', 'Oetjen_P', 'Oetjen_U', 'Sun_sample1_CS', 'Sun_sample2_KC', 'Sun_sample3_TB', 'Sun_sample4_TC']},
        ]
    },
    'toy': {
        'path': '/home/mohsen/data/toy/toy_normalized.h5ad',
        'condition_key': 'study',
        'cell_type_key': 'celltype',
        'runs': [
            {'reference': ['Batch1']},
            {'reference': ['Batch1', 'Batch2']},
            {'reference': ['Batch1', 'Batch2', 'Batch3']},
            {'reference': ['Batch1', 'Batch2', 'Batch3', 'Batch4']},
        ]
    },
}

def train_model(data_name, path, condition_key, cell_type_key, reference_conditions):
    results_path = f'/home/mohsen/scarches_review_results/results/ref_query_ratio/CVAE_MSE/{data_name}/{len(reference_conditions)}/'
    os.makedirs(results_path, 
                exist_ok=True)
    adata = sc.read(path)

    source_adata = adata[adata.obs[condition_key].isin(reference_conditions)]
    target_adata = adata[~adata.obs[condition_key].isin(reference_conditions)]
    target_conditions = target_adata.obs[condition_key].unique().tolist()

    network = sca.models.scArches(task_name=f"ref_query_ratio_{data_name}_{len(reference_conditions)}_reference",
                                  x_dimension=source_adata.shape[1], 
                                  z_dimension=10,
                                  architecture=[128],
                                  conditions=reference_conditions,
                                  gene_names=source_adata.var_names.tolist(),
                                  alpha=0.01,
                                  beta=0.0,
                                  eta=100.0,
                                  use_batchnorm=False,
                                  clip_value=3.0,
                                  loss_fn='mse',
                                  model_path=f"./models/CVAE/",
                                  dropout_rate=0.05,
                                  output_activation='relu')

    start_time = time.time()
    network.train(source_adata,
                  train_size=0.9,
                  condition_key=condition_key,
                  n_epochs=300,
                  batch_size=32,
                  early_stop_limit=15,
                  lr_reducer=10,
                  save=True,
                  retrain=True,
                  )
    end_time = time.time()
    with open(os.path.join(results_path, 'time.txt'), 'w') as f:
        f.write(f'reference_time: {end_time - start_time}')

    latent_adata = network.get_latent(source_adata, condition_key)
    latent_adata.write_h5ad(os.path.join(results_path, 'before.h5ad'))

    new_network = sca.operate(network, 
                             new_task_name=f'ref_query_ratio_{data_name}_{len(reference_conditions)}_after',
                             new_conditions=target_conditions,
                             version='scArches v1',
                             new_training_kwargs={'alpha': 0.023, 'eta': 1000.0},
                             remove_dropout=False)

    start_time = time.time()
    new_network.train(target_adata,
                      train_size=0.9, 
                      condition_key=condition_key,
                      n_epochs=100,
                      batch_size=32, 
                      early_stop_limit=15,
                      lr_reducer=10,
                      save=True, 
                      retrain=True,)
    end_time = time.time()

    latent_adata = new_network.get_latent(adata, condition_key)
    latent_adata.write_h5ad(os.path.join(results_path, 'all.h5ad'))

    with open(os.path.join(results_path, 'time.txt'), 'w') as f:
        f.write(f'query_time: {end_time - start_time}')



def run_all_tasks(datasets, data_name):
    dataset = datasets[data_name]

    for run in dataset['runs']:
        reference_conditions = run['reference']
        
        data_dict = dataset.copy()
        data_dict.pop('runs')
        data_dict['reference_conditions'] = reference_conditions
        
        train_model(data_name, **data_dict)

# run_all_tasks(datasets, 'pancreas')
run_all_tasks(datasets, 'mouse_brain')
# run_all_tasks(datasets, 'pbmc')
# run_all_tasks(datasets, 'toy')
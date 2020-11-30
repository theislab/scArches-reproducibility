import os
import torch
from scvi.core.modules import VAE

torch.set_printoptions(edgeitems=10, sci_mode=False, precision=4)

deep_cond = False
surgery_opt = 'freezed_expr'
#surgery_opt = 'freezed'
#surgery_opt = 'unfreezed'

check_for_retrained_weights = False
check_for_all = True
if deep_cond:
    dir_path = os.path.expanduser(
        f'~/Documents/benchmarking_results/figure_2/scvi/deep_cond/')
else:
    dir_path = os.path.expanduser(
        f'~/Documents/benchmarking_results/figure_2/scvi/first_cond/')

ref_model_path = f'{dir_path}reference/reference_model_state_dict'
ref_model = VAE(
    n_input=4000,
    n_batch=2,
    n_layers=2,
    use_batch_norm="none",
    use_layer_norm="both",
    encode_covariates=True,
    deeply_inject_covariates=deep_cond
)
ref_model.load_state_dict(torch.load(ref_model_path), strict=True)

surg_model_path = f'{dir_path}{surgery_opt}/surgery_model_state_dict'
surg_model = VAE(
    n_input=4000,
    n_batch=4,
    n_layers=2,
    use_batch_norm="none",
    use_layer_norm="both",
    encode_covariates=True,
    deeply_inject_covariates=deep_cond
)
surg_model.load_state_dict(torch.load(surg_model_path), strict=True)
if check_for_retrained_weights:
    for name, p in surg_model.named_parameters():
        for name_r, p_r in ref_model.named_parameters():
            if name == name_r:
                if not torch.equal(p, p_r):
                    print("\n\n")
                    print(name)
                    print(p.size(0), p.size(-1))
                    print(p[0])
                    print(p_r.size(0), p_r.size(-1))
                    print(p_r[0])
if check_for_all:
    for name, p in surg_model.named_parameters():
        for name_r, p_r in ref_model.named_parameters():
            if name == name_r:
                print("\n")
                print(name)
                if not torch.equal(p, p_r):
                    print("RETRAINED")

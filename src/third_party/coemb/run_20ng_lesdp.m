function run_20ng_lesdp(topics,dim_list,n_bst,procmode,ndocs)

wts = 0.001;

alg_params = wts;

alg_fun = 'run_sdp_allw';
alg_str = 'sdple';

for di=1:length(dim_list)
     dim = dim_list(di);
     main_all(topics,dim,n_bst,alg_str,alg_fun,alg_params,procmode,ndocs)
end

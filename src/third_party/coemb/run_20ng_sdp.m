function run_20ng_sdp(topics,dim_list,n_bst,procmode,ndocs)


alg_fun = 'run_code_sdp';
alg_str = 'sdp';

for di=1:length(dim_list)
     dim = dim_list(di);
     main_all(topics,dim,n_bst,alg_str,alg_fun,alg_params,procmode,ndocs)
end

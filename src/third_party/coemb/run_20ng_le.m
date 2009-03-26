function run_20ng_le(topics,dim_list,n_bst,procmode,wmarg,ndocs)

max_repeats = 1;
if ~wmarg
   alg_str = 'le';
   alg_fun = 'logemb_grad_dim';
else
   alg_str = 'mrgle';
   alg_fun = 'mrgle_grad_dim';
end

alg_params = max_repeats;

for di=1:length(dim_list)
     dim = dim_list(di);
     main_all(topics,dim,n_bst,alg_str,alg_fun,alg_params,procmode,ndocs)
end

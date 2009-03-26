function run_20ng_ca(topics,dim_list,n_bst,procmode,ndocs)

max_repeats = 10;
alg_str = 'log_cca';
alg_fun = 'log_cca';
%alg_str = 'ca_cca';
%alg_fun = 'ca_cca';

alg_params = '';

main_all(topics,dim_list,n_bst,alg_str,alg_fun,alg_params,procmode,ndocs);

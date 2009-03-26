function make_purity(tmpfn,group,dim,nw,ndocs,proc,newfn)

[topics,setname] = get_set(group);
[lik_ds,w,phi,psi] = read_sdp_tmp(tmpfn,dim);
tmp_psi = psi;
tmp_phi = phi;
psi = tmp_phi;
phi = tmp_psi;

options ='';
[pwd,fn,class_ids,chosen_words] = main_all(topics,dim,nw,'','',options,proc,ndocs);


[purity{dim}, purity_nrm{dim}]=test_embed(psi{dim},class_ids);  


if nargin>6
   save(newfn,'phi','psi','purity');
end

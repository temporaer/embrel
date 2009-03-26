function [phi,psi,c,lik] = run_lrl2(pxy0,dim_list,options)

max_restarts0 = 1;
[options,n_restarts]         = take_from_struct(options,'n_restarts',max_restarts0);
[options,prt_level]          = take_from_struct(options,'prt_level',1);
[options,cg_opts]            = take_from_struct(options,'cg_opts',-1*ones(1,9));
[options,b_log]            = take_from_struct(options,'b_log',1);
[options,zero_thr]            = take_from_struct(options,'zero_thr',0);
t0=cputime;


[NX,NY] = size(pxy0);
cg_params.NX = NX;
cg_params.NY = NY;

% If we are going to take the log, must threshold before
if b_log
  pxy0 = pxy0+zero_thr;
  pxy0 = pxy0/sum(pxy0(:));
end

for curr_dim = dim_list
    for(i_repeat=1:n_restarts)
        if(prt_level>0)
            fprintf('curr_dim=%d i_repeat=%3d \n', curr_dim,i_repeat);    
        end

	phi0 = rand(NX,curr_dim);
        psi0 = rand(NY,curr_dim);
	c0 = rand;
	
        x0 = [phi0(:);psi0(:);c0];
        cg_params.dim = curr_dim;
	px0 = sum(pxy0,2);
	py0 = sum(pxy0,1);
	% This is what should be modeled as a distance
	if b_log
	  cg_params.lrat = -log(pxy0./(px0*py0));
	else
	  cg_params.lrat = pxy0./(px0*py0);
	end
	
        cg_params.NX = NX;
	cg_params.NY = NY;	
	cg_params.b_log = b_log;
	
        x = conj_grad('lrl2_grad',cg_params,x0,cg_opts);
        [rep_phis{i_repeat},rep_psis{i_repeat},rep_cs{i_repeat}] = read_lrl2_params(x,cg_params);
        liks(i_repeat) = -lrl2_grad(x,cg_params);    
        if(prt_level>0)
            fprintf('Likelihood=%g\n',liks(i_repeat));
        end
    end
    liks = full(liks)
    [max_lik,bst_repeat] = max(liks);
    bst_repeat = bst_repeat(1);
    fprintf('choose best score out of repeats\n');  
    fprintf('max_lik=%f bst_repeat=%d\n',max_lik,bst_repeat);    
    
    
    lik(curr_dim) = max_lik;
    phi{curr_dim} = rep_phis{bst_repeat};
    psi{curr_dim} = rep_psis{bst_repeat};
    c{curr_dim} = rep_cs{bst_repeat};
end

return

% ===============================================================
function [out_options,val]=take_from_struct(options,fieldname,default)
%
% Take values from the options structure
%
out_options = options;    
if(isfield(options,fieldname));
    val = getfield(options,fieldname);
else
    val=default;
    out_options=setfield(out_options,fieldname,val);  
end

return

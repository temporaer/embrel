function  [phi,psi,sdp_liks,evs,wts,phis,psis] = run_code_sdp(nxy,curr_dim,params)

% Assume the field params.tmp_fname has been filled

NX = size(nxy,1);
NY = size(nxy,2);
ver = params.version;
wts = params.wts;

force_marg = 1;

for wi = 1:length(wts)
    w = wts(wi);
    if ver==1
      gram = code_sdp_solve(nxy,NX,NY,w);      
    elseif ver==2
      force_marg = 1;
      gram = code_sdp_solve5(nxy,NX,NY,w,force_marg,params);
    elseif ver==3
      gram = le_trace_gd_str(nxy,NX,NY,w);    
    elseif ver==4
      gram = code_sdp_solve_grad(nxy,NX,NY,w);
    elseif ver==5
      gram = gramcode_solve(nxy,NX,NY,w);      
    elseif ver==6
      precond = 0;      
      gram = code_v6(nxy,NX,NY,w,precond,params);      
    elseif ver==7
      precond = 1;
      gram = code_v6(nxy,NX,NY,w,precond,params);      
    elseif ver==8
      force_marg = 0;
      gram = code_sdp_solve5(nxy,NX,NY,w,force_marg,params);
    elseif ver==9
      if params.b_ycond
          gram = real(cond_sdp_code(nxy',NY,NX,w,params));      
          gram = gram';
      else
         gram = real(cond_sdp_code(nxy,NX,NY,w,params));      
      end
    elseif ver==10
         gram = real(code_nwtn_cond(nxy,NX,NY,w,params));      
    end


    % Generate reduced Gram matrix
    
    [v,d] = eig(gram);
    [v,d] = sortem(v,d);
    evs{wi} = diag(d);
    
    % Generate reduced Gram matrix
    newpoints = v(:,1:curr_dim)*sqrt(d(1:curr_dim,1:curr_dim));
    phis{wi} = real(newpoints(1:NX,:));
    psis{wi} = real(newpoints(NX+1:end,:));

    % Get its likelihood on the current data
    if ver==9
      sdp_liks(wi) =full(get_le_cond_lik(phis{wi},psis{wi},nxy,params.b_ycond));
    else
      sdp_liks(wi) =full(get_le_lik(phis{wi},psis{wi},nxy,force_marg));
    end
    
    save tmp_res
end


[v,i] = max(sdp_liks);
phi{curr_dim} = phis{i};
psi{curr_dim} = psis{i};
  

%plot([1 length(wts)],[base_lik base_lik],1:length(wts),sdp_liks);

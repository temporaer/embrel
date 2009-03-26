function [phi,psis,as,bs,liks] = run_mway_code(pxy0s,dim_list,options)

%
% function [PHI,PSIS,LIK]=CODE(PXY0S,DIM_LIST,OPTIONS)
%
% Implement euclidean embedding of co-occurence data, using
% conjugate gradient.
%
% Inputs: 
% 
% PXY0S      - Cell array of co-occurrence matrices (normalized to
% one). All matrices are assumed to have a common X (same number of rows).
% DIM_LIST - list of dimensions, for the low-dimensional embeddings
%            default is dim_list=[2];
% OPTIONS  - 
%     a structure with the following fields.
%     n_restarts   - repeat the algorithm n_restarts times and choose
%                    the one that maximizes the likelihood
%                    default = 20;
%     conjgrad     - options vector for conjgrad. see CONJGRAD below.
%     prt_level    - verbosity. default=1;
%
%
% Outputs: 
% PHIX - embeding of x objects. a cell array with entries for the
%     relevant dimensions only.
% PSIY - embeding of y objects. a cell array with entries for the
%     relevant dimensions only.
% L    - likelihood of solutions.
%
% 
% See
% Globerson A, Chechik G, Pereira F and Tishby N,  
% Euclidean embedding of co-occurence data, 
% Adavances in Neural Information Processing systems 17, (NIPS*2004)
%

%     n_restarts   - repeat the algorithm n_restarts times and choose
%                    the one that maximizes the likelihood
%                    default = 20;
%     conjgrad     - options vector for conjgrad. see CONJGRAD below.
%     prt_level    - verbosity. default=1;
%     margin_init -  How to set A(x),B(y) in the case when they're
%                    unchanged by the algorithm
%                    1 - constant, A(x)=B(y)=const(0) forall x,y
%                    2 - empirical, A(x)=p(x), B(y)=p(y);
%                    Default 2
%     w_nxy        - weight of each co-occurrence matrix. Determines its
%                    relative weight in the likelihood
%                    maximization. 
%                    Default: Uniform
%     phi0         - A fixed value to set phi(x) to. If phi0 does
%                    not exist or is empty, phi(x) is optimized
%                    over
%     b_keep_a     - If true, initializes A(x) according to
%                    margin_init and keeps it constant througout
%                    the optimization. 
%                    Default: 1
%     b_keep_b     - As b_keep_a, but for B(y).
%     
%     b_cond       - Make model conditional on x

% (C) G. Chechik, A. Globerson 2004


if ~iscell(pxy0s)
   npxy0s = {}
   npxy0s{1} = pxy0s;
   clear pxy0s;
   pxy0s = npxy0s;
   clear npxy0s;
end



max_restarts0 = 1;
[options,phi0] = take_from_struct(options,'phi0','');
cg_params.b_keep_phi = ~isempty(phi0);
cg_params.b_keep_psi = 0;
[options,cg_params.b_keep_a] = take_from_struct(options,'b_keep_a',1);
[options,cg_params.b_keep_b] = take_from_struct(options,'b_keep_b',1);
[options,n_restarts]         = take_from_struct(options,'n_restarts',max_restarts0);
[options,prt_level]          = take_from_struct(options,'prt_level',1);
[options,cg_params.w_nxy]    = take_from_struct(options,'w_nxy',[]);
[options,margin_init]        = take_from_struct(options,'margin_init',2);
[options,cg_opts]            = take_from_struct(options,'cg_opts',-1*ones(1,9));
[options,prt_level]          = take_from_struct(options,'prt_level',1);
[options,cg_params.b_cond]   = take_from_struct(options,'b_cond',0);

b_take_weights = isempty(cg_params.w_nxy);

t0=cputime;

if cg_params.b_cond & ~cg_params.b_keep_b
   if length(pxy0s)>1
      error('Unsupported');
   end
   pxy0s{1} = pxy0s{1}';
end

[NX,NY] = size(pxy0s{1});
cg_params.NX = NX;

for curr_dim = dim_list
    for(i_repeat=1:n_restarts)
        if(prt_level>0)
            fprintf('curr_dim=%d i_repeat=%3d \n', curr_dim,i_repeat);    
        end
        
        
        cg_params.NYs = NY;
        if ~cg_params.b_keep_phi
            phi0 = rand(NX,curr_dim);
        else
            phi0 = options.phi0;
        end
        
        x0 = [phi0(:)];
        cg_params.dim = curr_dim;
        cg_params.pxy0s = pxy0s;
        % Prepare marginals
        for i=1:length(pxy0s)
            cg_params.px0s{i} = sum(cg_params.pxy0s{i},2);
            cg_params.py0s{i} = sum(cg_params.pxy0s{i},1);      
            cg_params.NYs(i) = length(cg_params.py0s{i});
            if b_take_weights
                cg_params.w_nxy(i) = 1;
            end
            psi0 = rand(cg_params.NYs(i),curr_dim);
            x0 = [x0;psi0(:)];            
            
            if ~cg_params.b_keep_a
                a0 = rand(NX,1);
                x0 = [x0;a0(:)];
            else
                switch margin_init
                case 1
                    a0 = zeros(NX,1);
                case 2
                    a0 = log(cg_params.px0s{i});
                end
                cg_params.a0{i} = a0;
            end
            if ~cg_params.b_keep_b
                b0 = rand(cg_params.NYs(i),1);
                x0 = [x0;b0(:)];                
            else
                switch margin_init
                case 1
                    b0 = zeros(cg_params.NYs(i),1);
                case 2
                    b0 = log(cg_params.py0s{i});
                end
                cg_params.b0{i} = b0;                
            end
        end
        
        x = conj_grad('mway_logemb_grad',cg_params,x0,cg_opts);
        [rep_phis{i_repeat},rep_psis{i_repeat},rep_as{i_repeat},rep_bs{i_repeat}] = read_mway_params(x,cg_params);
        liks(i_repeat) = -mway_logemb_grad(x,cg_params);    
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
    if length(pxy0s)==1
       psis{curr_dim} = rep_psis{bst_repeat}{1};
    else
       psis{curr_dim} = rep_psis{bst_repeat};
    end
     
    as = rep_as{bst_repeat};
    bs = rep_bs{bst_repeat};
    
    % Normalize features in SDR case
    if ~cg_params.b_keep_a & ~cg_params.b_keep_b & length(psis)==1
        [phi,psis{1}] = norm_sdr_features(phi,psis{1});
    end
    
end

if cg_params.b_cond & ~cg_params.b_keep_b
   tmp_phi = phi;
   tmp_psi = psis;
   for curr_dim = dim_list
       phi{curr_dim} = tmp_psi{curr_dim};
       psis{curr_dim} = tmp_phi{curr_dim};
   end 
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

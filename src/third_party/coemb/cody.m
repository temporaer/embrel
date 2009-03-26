function [phi,psis,as,bs,liks] = cody(pxy0s,dim_list,options)

%
% function [PHI,PSIS,LIK]=CODE(PXY0S,dim,OPTIONS)
%
% Implement euclidean embedding of co-occurence data, using
% conjugate gradient.
%
% Inputs: 
% =======
% PXY0S    - Cell array of co-occurrence joint probability matrices
%          . All matrices are assumed to have a common X (same
%            number of rows).
% DIM      - Dimension of embedding, is dim_list=[2];
% OPTIONS  - 
%     a structure with the following fields.
%     n_restarts   - repeat the algorithm n_restarts times and choose
%                    the one that maximizes the likelihood
%                    default = 20;
%     conjgrad     - options vector for conjgrad. see CONJGRAD below.
%     prt_level    - verbosity. default=1;
%     w_nxy        - weight of each co-occurrence matrix. Determines its
%                    relative weight in the likelihood
%                    maximization. 
%                    Default: Uniform
%     phi0         - A fixed value to set phi(x) to. If phi0 does
%                    not exist or is empty, phi(x) is optimized
%                    over
%     b_cond       - Make model conditional on x
%     x_marg       - How to treat marginal of x. One of:  
%                    'F' - free. Optimized over
%                    'U' - fixed as Uniform
%                    'M' - fixed as the empirical marginals
%     y_marg       - How to treat marginal of y. see x_marg.
%
% Outputs: 
% ========
% PHIX - embeding of x objects. a cell array with entries for the
%        relevant dimensions only.
% PSIY - embeding of y objects. a cell array with entries for the
%        relevant dimensions only.
% L    - likelihood of solutions.
%
% 
% See
% Globerson A, Chechik G, Pereira F and Tishby N,  
% Euclidean embedding of co-occurence data, 
% Adavances in Neural Information Processing systems 17, (NIPS*2004)
%
% (C) G. Chechik, A. Globerson 2004


%     margin_init -  How to set A(x),B(y) in the case when they're
%                    unchanged by the algorithm
%                    1 - constant, A(x)=B(y)=const(0) forall x,y
%                    2 - empirical, A(x)=p(x), B(y)=p(y);
%                    Default 2
%     b_keep_a     - If true, initializes A(x) according to
%                    margin_init and keeps it constant througout
%                    the optimization. 
%                    Default: 1
%     b_keep_b     - As b_keep_a, but for B(y).
%     





if ~iscell(pxy0s)
   npxy0s = {}
   npxy0s{1} = pxy0s;
   clear pxy0s;
   pxy0s = npxy0s;
   clear npxy0s;
end

b_reverse = 0;
% Reverse x and y
if options.b_cond & strcmp(options.y_marg,'F')
  error('Unsupported');
end

max_restarts0 = 10;
[options,phi_init] = take_from_struct(options,'phi_init','');
[options,tmp_fname] = take_from_struct(options,'tmp_fname','');
options.tmp_fname = tmp_fname;
[options,psi_init] = take_from_struct(options,'psi_init','');
[options,a_init] = take_from_struct(options,'a_init','');
[options,b_init] = take_from_struct(options,'b_init','');
[options,pxx0] = take_from_struct(options,'pxx0','');
if ~isempty(pxx0);  pxx0 = pxx0/sum(pxx0(:)); end
[options,pyy0] = take_from_struct(options,'pyy0','');
if ~isempty(pyy0);  pyy0 = pyy0/sum(pyy0(:)); end
[options,w_pxx0] = take_from_struct(options,'w_pxx0',0);
[options,w_pyy0] = take_from_struct(options,'w_pyy0',0);
[options,cg_params.b_keep_phi] = take_from_struct(options,'b_keep_phi',0);
[options,cg_params.b_keep_psi] = take_from_struct(options,'b_keep_psi',0);
[options,n_restarts]       = take_from_struct(options,'n_restarts',max_restarts0);
[options,prt_level]        = take_from_struct(options,'prt_level',1);
[options,cg_params.w_nxy]  = take_from_struct(options,'w_nxy',[]);
[options,margin_init]      = take_from_struct(options,'margin_init',2);
[options,cg_opts]          = take_from_struct(options,'cg_opts',-1*ones(1,9));
[options,prt_level]        = take_from_struct(options,'prt_level',1);
[options,cg_params.b_cond] = take_from_struct(options,'b_cond',0);
[options,x_marg]           = take_from_struct(options,'x_marg','F');
[options,y_marg]           = take_from_struct(options,'y_marg','F');
[options,noise]           = take_from_struct(options,'noise',1);


switch(x_marg)
 case 'U', cg_params.b_keep_a = 1; margin_init_x=1;
 case 'M', cg_params.b_keep_a = 1; margin_init_x=2;  
 case 'F', cg_params.b_keep_a = 0; margin_init_x=0;
end

switch(y_marg)
 case 'U', cg_params.b_keep_b = 1; margin_init_y=1;
 case 'M', cg_params.b_keep_b = 1; margin_init_y=2;  
 case 'F', cg_params.b_keep_b = 0; margin_init_y=0;  
end

if options.b_cond
  cg_params.b_keep_a = 1;
  cg_params.b_keep_b = 1;
  margin_init_x = 1;
end

%[options,cg_params.b_keep_a] = take_from_struct(options,'b_keep_a',1);
%[options,cg_params.b_keep_b] = take_from_struct(options,'b_keep_b',1);


b_take_weights = isempty(cg_params.w_nxy);

t0=cputime;



[NX,NY] = size(pxy0s{1});
cg_params.NX = NX;


x0=[];

for curr_dim = dim_list
    for(i_repeat=1:n_restarts)
        x0=[];
        if(prt_level>0)
            fprintf('curr_dim=%d i_repeat=%3d \n', curr_dim,i_repeat);    
        end
        
        
        cg_params.NY = NY;
        cg_params.dim = curr_dim;
        cg_params.pxy0s = pxy0s;
	cg_params.pxx0 = pxx0;
	cg_params.pyy0 = pyy0;
	cg_params.w_pxx0 = w_pxx0;
	cg_params.w_pyy0 = w_pyy0;

	psi0 = 2*(0.5-rand(NY,curr_dim))*noise;
	x0 = [x0;psi0(:)];
	
        % Prepare marginals
        for i=1:length(pxy0s)
	  cg_params.px0s{i} = sum(cg_params.pxy0s{i},2);
	  cg_params.py0s{i} = sum(cg_params.pxy0s{i},1);      
	  cg_params.NXs(i) = length(cg_params.px0s{i});
	  if b_take_weights
	    cg_params.w_nxy(i) = 1 %w_nxy(i);
	  end
	  
	  phi0 = 2*(0.5-rand(cg_params.NXs(i),curr_dim))*noise; 
	  
	  x0 = [x0;phi0(:)];            
	  
	  if isempty(a_init)
	    a0 = 2*(0.5-rand(cg_params.NXs(i),1))*noise;
	  else
	    a0 = a_init;
	  end
	  
	  if ~cg_params.b_keep_a
	    x0 = [x0;a0(:)];
	  else
	    switch margin_init_x
	      case 1
		a0 = zeros(cg_params.NXs(i),1);
	      case 2
		a0 = log(cg_params.px0s{i});
	    end
	    cg_params.a0{i} = a0;
	  end
	  
	  if isempty(b_init)
	    b0 = 2*(0.5-rand(cg_params.NY,1))*noise;
	  else
	    b0 = b_init{i}
	  end
	  
	  if ~cg_params.b_keep_b
	    x0 = [x0;b0(:)];                
	  else
	    switch margin_init_y
	      case 1
		b0 = zeros(cg_params.NY,1);
	      case 2
		b0 = log(cg_params.py0s{i});
	    end
	    cg_params.b0{i} = b0;                
	  end
	end

	cg_params.w_nxy
	x = conj_grad('mway_logemb_grad_y',cg_params,x0,cg_opts);
%        x = bfgs_orig(x0,'mway_logemb_grad',1e-6,1e5,0,cg_params); 
        [curr_phis,curr_psis,curr_as,curr_bs] = read_mway_params_y(x,cg_params);

        liks(i_repeat) = -mway_logemb_grad_y(x,cg_params);    
        if(prt_level>0)
            fprintf('Likelihood=%g\n',liks(i_repeat));
	end
        if liks(i_repeat)>=max(liks)
           fprintf('Found max %g\n',liks(i_repeat));
           bst_phi = curr_phis;
           bst_psi = curr_psis;
           bst_as = curr_as;
           bst_bs = curr_bs;
           save(options.tmp_fname,'bst_phi','bst_psi','bst_as','bst_bs','i_repeat','liks');
        end
     end
    liks = full(liks)
    [max_lik,bst_repeat] = max(liks);
    bst_repeat = bst_repeat(1);
    fprintf('choose best score out of repeats\n');  
    fprintf('max_lik=%f bst_repeat=%d\n',max_lik,bst_repeat);    
    
    
    lik(curr_dim) = max_lik;
    phi{curr_dim} = bst_phi;
    if length(pxy0s)==1
       psis{curr_dim} = bst_psi{1};
    else
       psis{curr_dim} = bst_psi;
    end
     
    as = bst_as;
    bs = bst_bs;
    
    % Normalize features in SDR case
    if ~cg_params.b_keep_a & ~cg_params.b_keep_b & length(psis)==1
        [phi,psis{1}] = norm_sdr_features(phi,psis{1});
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

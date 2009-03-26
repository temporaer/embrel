function [best_x_inds,ns] = filter_terms(Nxy,mode,n_words,Nxw,w_vec);
%
% filter relevant terms according to one of several modes.  Leave only
% n_words that meat the criterion; Assume words are rows and documents
% are cols. 
%
% TF - Term frequency
% TC - Term count
% NTF - Term frequecnies - take low weights
% I1dH - One-symbol info by entropy
% I1dP - One-symbol info by distribution
% DIdH - difference of I1dH (for IB-SI)
% DIdP - difference of I1dP (for IB-SI)
%
%
%

if(nargin<5)
  w_vec = [1 -1];
end


% init
disp(sprintf('filter terms: using mode %s ',mode));
mode(end+1:4) = '_';
if(n_words>=size(Nxy,1))
  best_x_inds = [1:size(Nxy,1)];
  return
end


switch(mode)
 case 'TF__', % 
    Nx = sum(Nxy');
    [ns,is] = sort(Nx);
    best_x_inds = is(end-n_words+1:end);
    disp(sprintf('              min TF allowed is %d',ns(end-n_words+1)));
 
 case 'NTF_', % 
    Nx = sum(Nxy');
    [ns,is] = sort(Nx);
    best_x_inds = is(1:n_words);
    disp(sprintf('              max TF allowed is %d',ns(n_words)));    
    
  case 'TC__', % 
    appear = Nxy>0;
    Nx = sum(appear');
    [ns,is] = sort(Nx);
    best_x_inds = is(end-n_words+1:end);
    disp(sprintf('              min TC allowed is %d',ns(end-n_words+1)));
               
  case 'I1dH', % Information from one symbol: delta Entropies
    Ix = I1dH(Nxy);
    [ns,is] = sort(Ix);
    best_x_inds = is(end-n_words+1:end);
    disp(sprintf('              min I1dH allowed is %f',ns(end-n_words+1)));
    
  case 'I1dP', % Information from one symbol: delta Entropies
    Ix = I1dP(Nxy);
    size(Ix)
    [ns,is] = sort(Ix);
    best_x_inds = is(end-n_words+1:end);
    disp(sprintf('              min I1dP allowed is %f',ns(end-n_words+1)));
  case 'DIdP', % Delta-Information from one symbol    
    if(nargin<4)
      error('Nxw is undefined');
    end
    Ix_ony = I1dP(Nxy);
    Ix_onw = I1dP(Nxw);
    dIyw = w_vec(1)*Ix_ony+w_vec(2)*Ix_onw;        
    [ns,is] = sort(dIyw);
    lst = [1:9 10:10:50 100:100:length(ns)];
    lst = lst(find(lst<length(ns)));
    for(i=1+length(ns)-lst)
	disp(sprintf('%d %6.4f',i,ns(i)));
    end
    n_pos = length(find(ns>0));
    n_pos_words = min(n_words,n_pos);
    best_x_inds = is(end-n_pos_words+1:end);
    disp(sprintf('min DIdP allowed is %f',ns(end-n_words+1)));
    disp(sprintf('    DIdP n_words is %d',n_pos_words));    
    % dIyw(best_x_inds)'    
    
  case 'IG__',
    Ix = IG(Nxy);
    [ns,is] = sort(Ix);
    best_x_inds = is(end-n_words+1:end);
    disp(sprintf('min IG allowed is %f',ns(end-n_words+1)));
    
  case 'DG__',
    Ix_ony = IG(Nxy);
    Ix_onw = IG(Nxw);
    dIyw = Ix_ony-Ix_onw;
    [ns,is] = sort(dIyw);
    best_x_inds = is(end-n_words+1:end);
    disp(sprintf('min DG allowed is %f',ns(end-n_words+1)));
  
  otherwise,
    mode
    error('illegal mode in filter_terms');
end



return

% ======================
function Ix = I1dH(Nxy)

    % Init
    Ny = sum(Nxy);
    Nx = sum(Nxy');    
    N = sum(Ny);  
    Pxy = Nxy/N;
    Py = Ny/N;
    Px = Nx/N;    
    Ix = zeros(length(Nx),1);
    % H(Y)
    inds = find(Py>0);    
    Hy =  - sum(Py(inds).*log2(Py(inds)));
    % H(Y|X=x)    
    for(ix=1:length(Nx))
      Pxx = Px(ix);
      Pygx = Pxy(ix,:);
      inds = find(Pygx>0);      
      Pygx = Pygx(inds)/Pxx;
      Hygx = - sum(Pygx.*log2(Pygx));
      Ix(ix) = Pxx*(Hy - Hygx);
    end    
    
return

% ======================
function Ix = I1dP(Nxy)

    % Init 
    Pxy = sparse(Nxy/sum(Nxy(:)));
    Py = sum(Pxy);
    Px = sum(Pxy');
    nx = size(Pxy,1);
    clear Nxy;
    Ix = zeros(nx,1);
    p2=Py; 
    
    Pxy2 = sp_make_cond_dist(Pxy);
    clear Pxy;
    Pxy2 = Pxy2+eps;

    Pxy2 = Pxy2.*log(Pxy2./repmat(p2,nx,1));
    Pxy2 = diag(Px)*Pxy2;
    Ix = sum(Pxy2,2);
    
%    for(i_x=1:nx)      
%      i_x
%      p1=Pxy(i_x,:);
%      p1 = p1/sum(p1);
%      inds = find(p1>0);
%      Ix(i_x) = Px(i_x) * sum(p1(inds).*log(p1(inds)./p2(inds)));
%    end
%    keyboard
return

% ======================
function Ix = IG(Nxy)

    % Init
    Ny = sum(Nxy);
    Nx = sum(Nxy');    
    N = sum(Ny);  
    Pxy = Nxy/N;
    Py = Ny/N;
    Px = Nx/N;    
    % H(Y)
    inds = find(Py>0);    
    Hy =  - sum(Py(inds).*log2(Py(inds)));
    
    for(ix=1:length(Nx));    
      Pxx  = Px(ix);      
      Pygx = Pxy(ix,:);
      i_x  = find(Pygx>0);      
      Pygx = Pygx(i_x)/Pxx;
      Hygx = - sum(Pygx.*log2(Pygx));
      
      Pxnx  = 1-Px(ix);      
      Pygnx = Py-Pxy(ix,:);
      i_nx  = find(Pygnx>0);
      Pygnx = Pygnx(i_nx)/Pxnx;
      Hygnx = - sum(Pygnx.*log2(Pygnx));      
      Ix(ix) = Hy - Pxx* Hygx - Pxnx* Hygnx;    
    end    

    
return

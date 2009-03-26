function [pxy0,filename,class_ids,chosen_words] = main_all(topics,dim_list,n_bst,alg_str,alg_fun,alg_params,procmode,n_docs)
%
% function main_all(topics,dim_list,n_bst,alg_str,alg_fun,alg_params,procmode)
%
%rand('seed',0);

% Init
t0=cputime;
addpath('../../../Sources');

% Load data
% =========
[tags,mat,words] = load_20ng; % M,L,words;
dirlist = load_dirlist;

% Extract the relevant topic
% ============================
% alt = 1; cmp = [2:6] sci = [12:15]

n_topics = length(topics)
{dirlist{topics}}
infix = sprintf('-%d', topics)


disp(sprintf('main: seperate relevant topic'));
[nxy,class_ids]= extract_topics(mat,tags,topics);
clear mat tags;

if 1 
docinds = []; 
% Choose ndocs from each subject
for i=1:max(class_ids)
  class_docs = find(class_ids==i);
  docinds = [docinds;class_docs(1:n_docs)];
end

nxy  = nxy(:,docinds);
class_ids = class_ids(docinds);

%pxy0 = sparse(pxy0);
end

% Filter words
% ============
inds1 = filter_terms(nxy ,'TF',n_bst+100);nxy1= nxy (inds1,:);
inds2 = filter_terms(nxy1,'NTF',n_bst);nxy2= nxy1(inds2,:);
clear nxy1 nxy2;
all_inds = inds1(inds2);


chosen_words = words(all_inds);
eval(sprintf('save words%d_%s chosen_words',length(all_inds),infix));


if(length(all_inds)<2) disp('All filtered out!!!');return; end;

% Prepare Pxy 
% ===========
rel_nxy = nxy(all_inds,:);
clear nxy;
empty_cols=find(sum(rel_nxy)==0);
fill = rand(size(rel_nxy,1),1)/100;
%rel_nxy(:,empty_cols) = rand(size(rel_nxy(:,empty_cols)))/100;
rel_nxy(:,empty_cols) =fill(:,ones(1,length(empty_cols)));

rel_words = {words{all_inds}};clear words;
n_rel = length(rel_words);

rel_pxy = rel_nxy/sum(rel_nxy(:));
clear rel_nxy;
[pxy0,procstr] = preprocess(rel_pxy,procmode);
pxy0 = pxy0/sum(pxy0(:));
pxy0 = sparse(pxy0);


% Just generate the data
if isempty(alg_fun)
  filename='';
  return;
end

% Test alg
max_repeats=10;
[nx,ny] = size(pxy0);
fprintf('\nStart %s size(pxy0)=%dx%d cputime=%4.2fmin\n',...
	alg_str,nx,ny,(cputime-t0)/60);
if ~isempty(alg_params)
  % Create filename for backup. 
  filename = sprintf('tmp_CODE_%s_%s_topics%s_dim%d_nx_%d_ndocs_%d.mat',...
	alg_str,procstr,infix,dim_list(1),nx,n_docs);
  alg_params.tmp_fname = filename;
  [phi,psi,lik] = feval(alg_fun,pxy0,dim_list,alg_params);
else
  [phi,psi,lik] = feval(alg_fun,pxy0,dim_list);
end


% Save results to files 
n_dims = length(dim_list);
for(i_dim = 1:n_dims)
  dim = dim_list(i_dim);

  if iscell(psi)
     [purity{dim}, purity_nrm{dim}]=test_embed(psi{dim},class_ids);  
  else
     [purity{dim}, purity_nrm{dim}]=test_embed(psi,class_ids);  
  end

  
  if n_docs~=1000
    filename = sprintf('CODE_%s_%s_topics%s_dim%d_nx_%d_ndocs_%d.mat',...
	alg_str,procstr,infix,dim,nx,n_docs);
  else
    filename = sprintf('CODE_%s_%s_topics%s_dim%d_nx_%d.mat',...
	alg_str,procstr,infix,dim,nx);
  end

  save(filename,...
       'topics','dim_list',...
       'psi','phi', ...
       'purity','alg_params','alg_fun','procmode','n_bst','n_docs');
end

return

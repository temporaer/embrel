
procmode=0;

% sets={'-1-2-12','-2-3','-4-3' '-12-13','-12-13-14','-12-13-14-15','-18-19'};
% n_bst=1000;
% topics = [1 2 12];
% topics = [4 3  ];

dim_list=[2 3 4 5 6 8 10];

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

n_topics = length(topics);
{dirlist{topics}};
infix = sprintf('-%d', topics);
disp(infix);


disp(sprintf('main: seperate relevant topic'));
[nxy,class_ids]= extract_topics(mat,tags,topics);
clear mat tags;

% Filter words
% ============
%inds1 = filter_terms(nxy ,'TF',5000);nxy1= nxy (inds1,:);
%inds2 = filter_terms(nxy1,'TC',2000);nxy2= nxy1(inds2,:);
%inds3 = filter_terms(nxy2,'I1dP',200);nxy3= nxy2(inds3,:);

inds1 = filter_terms(nxy ,'TF',n_bst+100);nxy1= nxy (inds1,:);
inds2 = filter_terms(nxy1,'NTF',n_bst);nxy2= nxy1(inds2,:);
inds3 = filter_terms(nxy2,'I1dP',n_bst);nxy3= nxy2(inds3,:);
clear nxy1 nxy2 nxy3;

all_inds = inds1(inds2(inds3));
if(length(all_inds)<2) disp('All filtered out!!!');return; end;

% Prepare Nxy for SVD
rel_nxy = nxy(all_inds,:);
clear nxy;
empty_cols=find(sum(rel_nxy)==0);
fill = rand(size(rel_nxy,1),1)/100;
rel_nxy(:,empty_cols) = fill(:,ones(1,length(empty_cols)));
[nx,ny] = size(rel_nxy);
rel_words = {words{all_inds}};clear words;
n_rel = length(rel_words);
[pxy0,procstr] = preprocess(rel_pxy,procmode);


% Test SVD 
fprintf('Exec SVD\n');
[U,S,V] = svds(log(rel_nxy+1),max(dim_list));
n_dims = length(dim_list);

for(i_dim = 1:n_dims)    
  dim = dim_list(i_dim);  
  svd_psi = V(:,1:dim);
  svd_phi = U(:,1:dim);  
  % plot_embedding([],svd_psi,rel_nxy,n_topics,rel_words,0,1);
  [svd_purity{i_dim}, svd_purity_nrm{i_dim}]=test_embed(svd_psi,class_ids);
end
  
filename = sprintf('CODE_svd%s_topics%s_dims_nx_%d',procstr,infix,nx);
save(filename,...
     'topics','dim_list',...
     'svd_purity', 'svd_psi','svd_phi')
fig_save(filename);

return
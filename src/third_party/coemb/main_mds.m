
procmode=0;

%n_bst=2000;

metric='euclidean'
metric='cosine'
topics=[2 3 ];
dim_list=[2 3 4 5 6 8 10];

% Init
addpath('../../../Sources');

% Load data
% =========
[tags,mat,words] = load_20ng; % M,L,words;
dirlist = load_dirlist;

% Extract the relevant topic
% ============================
% alt = 1; cmp = [2:6] sci = [12:15]

n_topics = length(topics);
fprintf(1,'%s ',dirlist{topics})
fprintf('\n');
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

% Prepare Nxy for CA
rel_nxy = nxy(all_inds,:);
clear nxy;
empty_cols=find(sum(rel_nxy)==0);
fill = rand(size(rel_nxy,1),1)/100;
rel_nxy(:,empty_cols) = fill(:,ones(1,length(empty_cols)));

rel_words = {words{all_inds}};clear words;
n_rel = length(rel_words);
[nx,ny]=size(rel_nxy);
[rel_nxy,procstr] = preprocess(rel_pxy,procmode);

% Test MDS
fprintf('Exec MDS. size(nxy)=[%d %d]\n',size(rel_nxy));

D =pdist(rel_nxy',metric);
if(min(D)<0 & min(D)>-10^(-10)) D = D-min(D); end;
[mds_psi e] = cmdscale(D);
n_dims = length(dim_list);
for(i_dim = 1:n_dims)
  dim = dim_list(i_dim);  
  % plot_embedding([],mds_psi,rel_nxy,n_topics,rel_words,0,1);
  [mds_purity{i_dim}, mds_purity_nrm{i_dim}]=...
      test_embed(mds_psi(:,1:dim),class_ids);
end
mds_phi=[];
mds_psi = mds_psi(:,1:20);

filename = sprintf('CODE_mds%c%s_topics%s_dims_nx_%d',...
		   metric(1),procstr,infix,nx);
save(filename,...
     'n_bst', 'topics','dim_list',...
     'mds_psi','mds_phi', ...     
     'mds_purity');
% fig_save(filename);

return

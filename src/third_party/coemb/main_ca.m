procmode=0;

% 
% Noam's topics
% =============
% sci crypt, electronics, med, space [12 13 14 15]
% talk politicsmideast politics misc [18 19]
%
% comp.graphincs rec.motorcycles rec.sport.baseball
% sci.space talk.politics.mideast
% 



%topics= [2 3 4 5 6];
% n_bst=1000;
topics = [12 13 14 15];
topics = [1 2 12];
dim_list=[2 3 4 5 6 8 10];

% Init
t0=cputime;
addpath('../../Sources');
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
[nx,ny] = size(rel_nxy);
[rel_nxy,procstr] = preprocess(rel_pxy,procmode);

% Test CA
fprintf('Exec CA\n');
n_dims = length(dim_list);
for(i_dim = 1:n_dims)
  dim = dim_list(i_dim);  
  fprintf('dim=%d\n',dim);    
  [ca_phi,ca_psi,Du] = corr_anl(rel_nxy,dim+1);
  [ca_purity{i_dim}, ca_purity_nrm{i_dim}]=test_embed(ca_psi,class_ids);
end

filename = sprintf('CODE_ca%s_topics%s_dims_nx_%d',procstr,infix,nx);
save(filename,...
     'topics','dim_list',...
     'ca_psi','ca_phi',...
     'ca_purity');

return
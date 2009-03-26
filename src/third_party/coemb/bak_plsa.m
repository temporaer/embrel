function [pwgz,pdgz,pz] = plsa(train,dim)

niter = 100;

[d_tr,w_tr,v_tr] = find(train);

NW = size(train,2);
ND = size(train,1);
NTRAIN = length(v_tr);

%Initalize everything uniformly
pz = ones(dim,1)/dim;
pwgz = rand(dim,NW)/NW;
pdgz = rand(dim,ND)/ND;
pwgz = sp_make_cond_dist(pwgz);
pdgz = sp_make_cond_dist(pdgz);

del = 1e-6;
w_mat = sparse(NTRAIN,NW);

for wi=1:NW
   where_w{wi} = find(w_tr==wi);
   w_mat(:,wi) = sparse(w_tr==wi);
end
 
d_mat = sparse(NTRAIN,ND);
for di=1:ND
   where_d{di} = find(d_tr==di);
   d_mat(:,di) = sparse(d_tr==di);   
end

lik_hist = 0;

for ni=1:niter
   lik = 0 ;
   
   pzgdw = repmat(pz,1,NTRAIN).*pdgz(:,d_tr).*pwgz(:,w_tr);
   lik = dot(v_tr,log(sum(pzgdw)));
   pzgdw = sp_make_cond_dist(pzgdw')';

   % Now multiply each instance by its count
   pzgdw = pzgdw*spdiags(v_tr,0,NTRAIN,NTRAIN);

   tmp_pwgz = pzgdw*w_mat;
   tmp_pdgz = pzgdw*d_mat;
   
   if 0
   for wi=1:NW
      tmp_pwgz(:,wi) = sum(pzgdw(:,where_w{wi}),2);
   end
   for di=1:ND
      tmp_pdgz(:,di) = sum(pzgdw(:,where_d{di}),2);
    end
  end
   tmp_pz = sum(pzgdw,2);

if 0
   tmp_pwgz = zeros(size(pwgz));
   tmp_pdgz = zeros(size(pdgz));
   tmp_pz = zeros(size(pz)); 

   % Calculate p(z|d,w) for training data
   for i=1:length(v_tr)
     v = pz.*pdgz(:,d_tr(i)).*pwgz(:,w_tr(i));
     lik = lik+v_tr(i)*log(sum(v));
     pzgdw = v/sum(v);
     tmp_pwgz(:,w_tr(i)) = tmp_pwgz(:,w_tr(i))+ v_tr(i)*pzgdw;
     tmp_pdgz(:,d_tr(i)) = tmp_pdgz(:,d_tr(i))+ v_tr(i)*pzgdw;
     tmp_pz = tmp_pz + v_tr(i)*pzgdw;
   end
end
   % Normalize
   pz = tmp_pz /sum(tmp_pz);
   pwgz = sp_make_cond_dist(tmp_pwgz);
   pdgz = sp_make_cond_dist(tmp_pdgz);

   fprintf(1,'Iter=%d Lik=%g\n',ni,lik);
   lik_hist(ni) = lik;
   if ni>1 & lik_hist(end)-lik_hist(end-1)<del
      break;
   end
end

function [lambdas,liks,perp] = delinterp(train_dat,heldout,new_model,lam0)

% Find optimal parameters for
% p(w|h) = \lambda_1 p_TD(w|h) + \lambda_2 p_TD(w) + \lambda_3 c
% where p_TD are the empirical counts.
% \lambda's are optimized to maximize likelihood over heldout data
% new_model gives the probability of another model on the heldout
% data. When 0, it is ignored

trained_bi = sp_make_cond_dist(train_dat);
trained_uni = sum(train_dat,2);
trained_uni = trained_uni/sum(trained_uni(:));
NW = size(train_dat,2);

if new_model==0
  n_lambda = 3;
else
  n_lambda = 4;
end

if lam0==0
  lambdas = ones(1,n_lambda)/n_lambda;
  niter = 100;
else
  lambdas = lam0;
  niter = 1;
end


[iho,jho,vho] = find(heldout);
k = 1/NW;
NHO = length(iho);
liks = [];
tol = 1e-2;

% Get all bigram probabilities
bi_fact = trained_bi(find(heldout));
uni_fact = trained_uni(jho);
kv = ones(size(uni_fact))*k;

N = sum(vho);

for ni=1:niter
  if (n_lambda==3)
    v = [lambdas(1)*bi_fact lambdas(2)*uni_fact lambdas(3)*kv];
  else
    v = [lambdas(1)*bi_fact lambdas(2)*uni_fact lambdas(3)*kv lambdas(4)*new_model];    
  end

  p = sum(v,2);
  v = sp_make_cond_dist(v);
  % Multiply by counts
  cnt = spdiags(vho,0,NHO,NHO)*v;
  lik = dot(vho,log(p));
  cnt = sum(cnt);
  last_lam = lambdas;
  lambdas = cnt/sum(sum(cnt));

  if sum(lambdas==0)==0
    rat = last_lam./lambdas;
    rat = max(max(rat,1./rat));
  else
    rat = 1e9;
  end
  liks(end+1) = lik;
  perp = exp(-lik/N);
  fprintf(1,'Perp=%g Lik=%g. Rat=%g\n',full(perp),full(lik),full(rat));
  if (rat<1+tol)
    break;
  end
  
end














if 0  
  lik = 0;
  cnt = zeros(size(lambdas));
  for i=1:NHO
    v = [lambdas(1)*trained_bi(iho(i),jho(i)) ...
	  lambdas(2)*trained_uni(jho(i)) lambdas(3)*k];
    % This is the probability of the current point according to the
    % mixture model
    p = sum(v);
    v = v/p;
    curr_count = vho(i);
    cnt = cnt+curr_count*v;
    lik = lik+curr_count*log(p);
  end
  lambdas = cnt/sum(cnt);
end
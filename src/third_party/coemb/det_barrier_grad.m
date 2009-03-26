function [bar_grad,bar_hess] = det_barrier_grad(x,N)

% Derivative and hessian of -log det(X)
fx = vec_to_smat_fast(x,N);

%ifx = inv(fx);

N = length(fx);

c = chol(fx);
sqifx = inv(c);
ifx = sqifx*sqifx';
bar_hess = zeros(length(x),length(x));


ind1 = 0; ind2=0;
for i1=1:N
  for j1=i1:N
    ind1 = ind1+1;
    ind2 = 0;
    for i2=1:N
      for j2=i2:N
	ind2 = ind2+1;
	fact =2;
	if i1==j1 & i2==j2 
	  fact = 0.5;
	elseif i1==j1 | i2==j2
	  fact = 1;
	end
	 
	bar_hess(ind1,ind2) = fact*(ifx(i1,j2)*ifx(j1,i2) + ifx(i1,i2)*ifx(j1,j2));
      end
    end
  end
end


%bar_grad = bar_grad';

N = size(fx,1);
mf = ones(N)*2;
mf = mf - eye(N);
g2 = -ifx;
bar_grad = smat_to_vec_c(mf.*g2)';



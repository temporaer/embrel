function bar_grad = det_barrier_grad_gd(fx)

%N = sqrt(size(F,1)); %size(F{1},1);
N = size(fx,1);
% Derivative and hessian of -log det(X)
%fx = F*x; %vec_to_smat(x,F);
%fx = reshape(fx,N,N);

%ifx = inv(fx);
mf = ones(N)*2;
mf = mf - eye(N);
bar_grad = -inv(fx); %
%bar_grad = uptri_to_x(mf.*bar_grad)';
bar_grad = smat_to_vec_c(mf.*bar_grad)';

return;

uind = triu(ones(N,N));
uind = find(uind');  % Indices indicating upper triangle

ndiag_inds = triu(ones(N,N),1);
ndiag_inds = find(ndiag_inds');

diag_inds = eye(N);
diag_inds = find(diag_inds);

rowind = 0;
bar_hess = zeros(length(F),length(F));

for i1=1:N
    for i2=i1:N
        rowind = rowind+1;
        curr_row = kron(ifx(i1,:),ifx(i2,:))+kron(ifx(i2,:),ifx(i1,:));
        if i1==i2
            curr_row(diag_inds) = curr_row(diag_inds)/2;
            curr_row = curr_row(uind);            
            bar_hess(rowind,:) = curr_row;
        else
            curr_row(ndiag_inds) = curr_row(ndiag_inds)*2;
            curr_row = curr_row(uind);            
            bar_hess(rowind,:) = curr_row;
        end
    end
end
a=1;














if 0 
    N = size(F{1},1);
    uind = triu(ones(N,N));
    uind = find(uind');  % Indices indicating upper triangle
    kr = kron(ifx,ifx);
    kr = kr(uind,uind);
    kr = kr*2;
    % mf = ones(size(kr))*2;
    % mf = mf - eye(length(kr));
    % kr = kr.*mf;
    
    krmf = kron(ifx,ifx);
    krmf = krmf(uind,uind);
    a=mf(uind);
    mff=a*a';  %repmat(a,1,45);
    krmff = krmf.*mff;
    
    
    for i=1:length(F)
        %    bar_grad2(i) = -trace(F{i}*ifx);
        for j=1:length(F)
            bar_hess2(i,j) = trace(ifx*F{i}*ifx*F{j});      
        end
    end
    
end
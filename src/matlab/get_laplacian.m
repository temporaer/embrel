function L = get_laplacian(M)
	n = size(M,1);
	D = diag(sum(M)); % degree matrix
	L = D - M;        % laplacian
	L = eye(n,n) - D^-.5 * M * D^-.5;   % normalized laplacian
	L = M;  % adjmat
end


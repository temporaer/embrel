function [final_gram,final_d,final_z] = run_csdp(alld,NX,NY,b_zmax)
% A = mat_completion(Ainit,indices,C);
% dxy - Given distances between x and y. Must be kept up to a constant. All others are free

% There are NX^2+NY^2+1 variables, for the Gram of X,Y and the shift variable

dxy = alld(1:NX,NX+1:end);


N = NX+NY;
PSDN = N + (NX*NY*4)+1;
% Z is the last element of P

if b_zmax
  m = sparse(PSDN,PSDN);
  m(PSDN,PSDN) = -1;
else
  % Minimize trace
  m = speye(PSDN,PSDN);
  m(PSDN,PSDN) = 0;
end

% This means we minimize Z only

C(:,1) = m(:);
A=[];
b=[];
blk={};
A = sparse(PSDN^2,8*NX*NY+1);
ind = 0;
curr_base = (NX+NY+1);

for xi=1:NX
    for yi=1:NY
        fprintf('xi=%d yi=%d\n',xi,yi);
        % Set the constraint P(xi,xi)+P(NX+yi,NX+yi)-2*P(xi,NX+yi)==Z*dxy(xi,yi)-Z);	
        % D =  [Z -1; dxy0 dxy] 
        % Set d(1,1) = c;    
        m = sparse(PSDN,PSDN);        
        m(curr_base,curr_base) = 1;
        m(PSDN,PSDN) = -1;
        ind = ind +1;
    	A(:,ind) = m(:);
    	b(end+1) = 0;
        
        % Set d(2,2) = dxy_model
        m = sparse(PSDN,PSDN);
        m(curr_base+1,curr_base+1) = -1;
    	m(xi,xi) = 1;        
		m(NX+yi,NX+yi) = 1;
		m(xi,NX+yi) = -1;
		m(NX+yi,xi) = -1;	
    	b(end+1) = 0;        
        ind = ind+1;
    	A(:,ind) = m(:);
        
        
        % Set d(2,1) = sqrt(dxy0)
        m = sparse(PSDN,PSDN);
        m(curr_base,curr_base+1) = sqrt(dxy(xi,yi));
        ind = ind + 1;
    	A(:,ind) = m(:);
        b(end+1) = 1;

        % Set d(2,1) = sqrt(dxy0)        
        m = sparse(PSDN,PSDN);
        m(curr_base+1,curr_base) = sqrt(dxy(xi,yi));
        ind = ind + 1;
    	A(:,ind) = m(:);
        b(end+1) = 1;
        
        curr_base = curr_base + 2;

        m = sparse(PSDN,PSDN);        
        m(curr_base,curr_base) = 1;
        m(PSDN,PSDN) = -1;
        ind = ind +1;
    	A(:,ind) = m(:);
    	b(end+1) = 0;
        
        % Set d(2,2) = dxy_model
        m = sparse(PSDN,PSDN);
        m(curr_base+1,curr_base+1) = -1;
    	m(xi,xi) = 1;        
		m(NX+yi,NX+yi) = 1;
		m(xi,NX+yi) = -1;
		m(NX+yi,xi) = -1;	
    	b(end+1) = 0;        
        ind = ind+1;
    	A(:,ind) = m(:);
        
        
        % Set d(2,1) = sqrt(dxy0)
        m = sparse(PSDN,PSDN);
        m(curr_base,curr_base+1) = sqrt(dxy(xi,yi));
        ind = ind + 1;
    	A(:,ind) = m(:);
        b(end+1) = 1;

        % Set d(2,1) = sqrt(dxy0)        
        m = sparse(PSDN,PSDN);
        m(curr_base+1,curr_base) = sqrt(dxy(xi,yi));
        ind = ind + 1;
    	A(:,ind) = m(:);
        b(end+1) = 1;
    end
end

% Add centering requirement
m = ones(PSDN,PSDN);
m(NX+NY+1:end,NX+NY+1:end) = 0;
ind = ind +1;
A(:,ind) = m(:);
b(end+1) = 0;
%x0 = allg;
%x0(N+1,N+1) = 0;
%x0 = x0+eye(N+1)*1e-8;
%y0 = zeros(1,length(A));
%z0 = x0;

%sqlparameters;
K.s = PSDN;
[final_gram] = csdp(A,b',C,K);
final_gram = reshape(final_gram,PSDN,PSDN);
final_z = final_gram(PSDN,PSDN);
final_gram = final_gram(1:(NX+NY),1:(NX+NY));

final_d = g_to_d(final_gram);

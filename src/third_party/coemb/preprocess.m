function [processed,str] = preprocess(data,mode,quiet);
%
%
%

if(nargin<3)
  quiet=0;
end

switch(mode)
 case 0, 
  str = 'none';
  if(quiet>0),  fprintf('No processing\n');end;
  processed = data;
 case 1, 
  str = 'margy';
  if(quiet>0),  fprintf('Normalized Y marginals\n');end;
  processed = make_cond_dist(data,1);
 case 2, 
  str = 'copul';
  if(quiet>0),  fprintf('Normalized BOTH X and Y marginals\n');end;
  processed = make_cond_dist(data,2);  
 case 3, 
  str = 'tfidf';
  if(quiet>0),  fprintf('TFIDF normalization\n');end;
  processed = tfidf(data);  
 case 4, 
  str = 'ltctd';
  if(quiet>0),  fprintf('LTC-TFIDF normalization\n');end;
  error('Not implemented yet\n');
 case 5, 
  str = 'logtfidf';
  if(quiet>0),  fprintf('LogTFIDF normalization\n');end;
  data = make_cond_dist(data,1);
  processed = logtfidf(data);  
 otherwise, 
  mode
  error('Illegal mode\n');
end
return



% =============================================
function pxy_cond = make_cond_dist(pxy,bcond_x)

% Assume x is given in the rows
pxy_cond = zeros(size(pxy));
switch(bcond_x)
 case 0,
  px = sum(pxy,2);
  for i=1:size(pxy,1)
    pxy_cond(i,:)=pxy(i,:)/px(i);
  end  
  
 case 1,
  py = sum(pxy,1);
  for i=1:size(pxy,2)
    pxy_cond(:,i)=pxy(:,i)/py(i);
  end
  
 case 2, 
  px = sum(pxy,2);
  py = sum(pxy,1);  
  for i=1:size(pxy,1)
    for j=1:size(pxy,2)    
      pxy_cond(i,j)=pxy(i,j)/(px(i)*py(j));
    end
  end    
  
end
return


% =============================================
function processed = tfidf(data)
%
%
%

processed = zeros(size(data));
df = sum(data>0,2);
for(i_row=[1:size(data,1)])
  processed(i_row,:) = data(i_row,:)/df(i_row);
end

return

% =============================================
function processed = logtfidf(data)
%
%
%

processed = zeros(size(data));
df = sum(data>0,2);
df = size(data,2)./df;

for(i_row=[1:size(data,1)])
  processed(i_row,:) = data(i_row,:)*log2(df(i_row));
end

return


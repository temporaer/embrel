function m = vec_to_smat(v,F)

m  = zeros(size(F{1}));
for i=1:length(F)
    m = m+v(i)*F{i};
end
% ind = 1;
% for k=0:N-1
%     m = m+diag(v(ind:ind+N-k),k);
%     if k>0
%         m = m+diag(v(ind:ind+N-k),k);        
%     end
%     ind = ind+N-k;
% end

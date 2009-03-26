function op = get_gram_op_inds(n,xinds,yinds)

ind = 0;
op=[];

% Only go over half the matrix
for indi=1:length(xinds)
  i = xinds(indi);
  j = yinds(indi);
  d = sparse(n,n);
  if i~=j
    d(i,j) = -2;
    d(j,i) = -2;
    d(i,i) = 1;
    d(j,j) = 1;
  end      
  op(end+1,:) = smat_to_vec_c(full(d));
end

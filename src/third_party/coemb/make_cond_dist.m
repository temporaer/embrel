function pxy_cond = make_cond_dist(pxy,bcond_x)

% Assume x is given in the rows
pxy_cond = zeros(size(pxy));
if bcond_x
  py = sum(pxy,1);
  for i=1:size(pxy,2)
      pxy_cond(:,i)=pxy(:,i)/py(i);
  end
else
  px = sum(pxy,2);
  for i=1:size(pxy,1)
      pxy_cond(i,:)=pxy(i,:)/px(i);
  end
end


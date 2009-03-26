function [x,y] = joint_to_indic(nxy)

sx = size(nxy,1);
sy = size(nxy,2);
n = sum(nxy(:));

x = zeros(sx,n);
y = zeros(sy,n);

ind = 1;
for xi = 1:sx
  for yi = 1:sy
    x(xi,ind:ind+nxy(xi,yi)-1) = 1;
    y(yi,ind:ind+nxy(xi,yi)-1) = 1;    

    ind = ind + nxy(xi,yi);
  end
end


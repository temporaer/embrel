function [x,y] = joint_to_indic_red(nxy)

sx = size(nxy,1);
sy = size(nxy,2);
n = sum(nxy(:));

x = zeros(sx-1,n);
y = zeros(sy-1,n);

ind = 1;
for xi = 1:sx
  for yi = 1:sy
    if xi<sx
      x(xi,ind:ind+nxy(xi,yi)-1) = 1;
    end
    if yi<sy
      y(yi,ind:ind+nxy(xi,yi)-1) = 1;    
    end
    ind = ind + nxy(xi,yi);
  end
end


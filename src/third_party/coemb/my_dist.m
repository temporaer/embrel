function d = my_dist(x,y)

Xs = size(x,1);
Ys = size(y,2);

 x_norm = sum(x.^2,2);
 y_norm = sum(y.^2,1);
 
 d = x*y;
 d = repmat(x_norm,1,Ys) + repmat(y_norm,Xs,1) -  2*d ;

%for xi=1:size(x,1)
%    for yi=1:size(y,2)
%%        d(xi,yi) = sum((x(xi,:)-y(:,yi)').^2);
%    end
%end

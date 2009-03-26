function pn=draw_from_p(p,n)
% p is assumed one dimensional

u=rand(1,n);
pcs=cumsum(p);

for i=1:n
    ind=max(find(u(i)>pcs));
    if isempty(ind)
        ind=0;
    end
    pn(i)=ind+1;
end




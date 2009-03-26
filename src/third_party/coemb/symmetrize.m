function newm = symmetrize(m)

N = size(m,1);
newm = triu(m);
newm = newm + triu(m,1)';

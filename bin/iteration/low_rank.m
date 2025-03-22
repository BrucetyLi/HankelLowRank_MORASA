function Ik = low_rank(I0,opt)
[u,s,v] = svds(I0,opt.J);

Ik = u*s*v';

end
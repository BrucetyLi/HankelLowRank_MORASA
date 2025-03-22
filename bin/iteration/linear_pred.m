function Ik = linear_pred(I0,opt)
row = opt.ny;
col = opt.nx;
echo = opt.nt;
Ik = zeros(size(I0));
niu = opt.niu;

for i=1:row*col
temp = I0(i,:);
H1 = hankel(temp(1:floor(echo/2)),temp(floor(echo/2):end));
[u,s,v] = svd(H1);
s = max(s-niu,0);
H2 = u*s*v';
H2 = H2(end:-1:1,:);

for j= 1:echo
Ik(i,j) = mean(diag(H2,j-floor(echo/2)));
end

end

end

function Ik = Sparsity_prior(I0,opt)
row = opt.ny;
col = opt.nx;
echo = opt.nt;

Ik = zeros(size(I0));
[WC,WS] = wavedec2(reshape(I0(:,1),[row,col]),4,'db4');
alpha = zeros(size(WC,2),echo);
opt.waveletS = WS;

for i=1:echo
[WC,~] = wavedec2(reshape(I0(:,i),[row,col]),4,'db4');
alpha(:,i) = permute(WC,[2,1]);
end

alpha = JointSoftThresh(alpha,opt.tau);

for i=1:echo
Ik(:,i) = conj(reshape(waverec2(alpha(:,i),opt.waveletS,'db4'),[row*col,1]));
end


end
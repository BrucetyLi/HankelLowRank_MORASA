function Ik = data_consistency(I0,mat,opt)
row = opt.ny;
col = opt.nx;
echo = opt.nt;
channel = opt.nc;

C = reshape(opt.maps,[row*col,1,channel]);
U = reshape(opt.mask,[row,col,echo,1]);
vec1 = C.*reshape(I0,[row*col,echo,1]);
vec2 =(1-U).*fftc(fftc(reshape(vec1,[row,col,echo,channel]),1),2)+ mat;
Ik = sum(reshape(ifftc(ifftc(vec2,1),2),[row*col,echo,channel]).*conj(C),3);

end
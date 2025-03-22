function arr = EPGarray(Tm,opt)
row = opt.ny;
col = opt.nx;
echo = opt.nt;
arr = zeros(row*col,echo);
for i = 1:(row*col)
arr(i,:) = FSEsig(Tm(i,1),Tm(i,2),opt);
end

end
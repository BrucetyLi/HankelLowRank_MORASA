function maps = Sense(mat,num_low)
windows = hann(num_low)*hann(num_low)';
mat = sum(mat,3);
[row,col,channel] = size(mat);
data = zeros(row,col,channel);
data(ceil(0.5*row-0.5*num_low+1):floor(0.5*row+0.5*num_low), ...
    ceil(0.5*col-0.5*num_low+1):floor(0.5*col+0.5*num_low), ...
    :)= mat(ceil(0.5*row-0.5*num_low+1):floor(0.5*row+0.5*num_low), ...
    ceil(0.5*col-0.5*num_low+1):floor(0.5*col+0.5*num_low),:).*windows;

data = ifftc(ifftc(data,2),1);
maps = data./sqrt(sum(abs(data).^2,3));

end
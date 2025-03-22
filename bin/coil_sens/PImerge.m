function result = PImerge(C,temp)
[row,col,frame,channel] = size(temp);
C1 = repmat(reshape(C,[row,col,1,channel]),[1,1,frame,1]);
result = sum(conj(C1).*temp,4);%./repmat(sum(abs(C).^2,3),[1,1,frame]);
end



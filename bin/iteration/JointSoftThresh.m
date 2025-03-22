function y = JointSoftThresh(x,tau)
temp1 = sqrt(sum(abs(x).^2,2));
temp1(temp1==0) = 1e-5;
y = x./temp1.*max(0,temp1-tau);
end
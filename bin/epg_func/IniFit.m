function [T2, B1, amp] = IniFit(S,opt)

    [X,~,~,~,~] = lsqnonlin(@diff_sig,opt.lsq.Icomp.X0,...
        opt.lsq.Icomp.XL,opt.lsq.Icomp.XU,opt.lsq.fopt);
    T2 = X(1); amp = X(2); B1 = X(3);
      
function delta = diff_sig(X)
    
     Sfit = X(2).*FSEsig(X(1),X(3),opt);
     delta = abs(S(:) - Sfit(:));
    
end
end

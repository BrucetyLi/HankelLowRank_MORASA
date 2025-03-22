function S = FSEsig(T2,B1,opt)
%   Variable to replicate flip angles
if isfield(opt,'FA_array')
    replicate = opt.FA_array';
else
    replicate = ones(opt.etl,1);
end


switch lower(opt.mode(1))
    
    case 'n'
        FA = B1*pi/180*opt.RFr.angle;
        FA = replicate*FA;
        S = epgMEX(opt.T1(T2),T2,opt.esp,FA);
%         S = S(2:2:end); %
    case 's'

        Me = sin(B1*opt.RFe.alpha);
        FA = replicate*(B1*opt.RFr.alpha);
        Mr = epgMEX(opt.T1(T2),T2,opt.esp,FA);
        Mnet = (ones(opt.etl,1)*Me) .* Mr;
        S = sum(Mnet,2);
        S = S./opt.Nz;
        
end
        




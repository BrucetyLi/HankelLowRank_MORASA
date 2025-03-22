function opt = optset(mode)
if nargin < 1 || (~strcmp(mode,'s') && ~strcmp(mode,'n'))
    disp('Assuming non-selective refocusing');
    opt.mode = 'n';
else
    opt.mode = mode;
end

opt.rho1 = 1;
opt.rho2 = 0.01;
opt.rho3 = 1;
opt.lambda1 = 1e-5;
opt.iter = 60;

opt.esp = 10e-3;            %   Inter-echo spacing (s)
opt.etl = 20;               %   Echo train length
opt.RFe.angle = 90;         %   Prescribed excitation pulse angle (degrees)
opt.RFr.angle = 180;        %   Prescribed refocusing pulse angle (degrees)

opt.T1  = @(T2)3.0;     %   T1 value (s)

if strcmpi(opt.mode,'s')
    opt.Dz  = [-0.5 0.5];   %   Spatial bounds (cm) (should exceed slice thickness, can start at 0 for half profile)
    opt.Nz  = 51;           %   Number of positions to simulate across slice
    opt.Nrf = 64;           %   Number of resampled points in the RF waveform

    opt.RFe.RF = [];        %   Placeholder for RF waveform
    opt.RFe.tau = 2.000e-3; %   Duration of excitation pulse (s)
    opt.RFe.G = 0.5000;     %   Slice select gradient during excitation (G/cm)
    opt.RFe.phase = 0;      %   Relative phase (0 in CPMG) (degrees)
    opt.RFe.ref = 1.00;     %   Refocusing fraction (x2, i.e. near unity for excite; zero for refocus)
    opt.RFe.alpha = [];     %   Actual tip angle distribution (degrees) (computed later...)

    opt.RFr.RF = [];        %   Placeholder for RF wavform
    opt.RFr.tau = 2.000e-3; %   Duration of excitation pulse (s)
    opt.RFr.G = 0.5000;     %   Slice select gradient during excitation (G/cm)
    opt.RFr.phase = 90;     %   Relative phase (90 in CPMG) (degrees)
    opt.RFr.ref = 0.00;     %   Refocusing fraction (x2, i.e. near unity for excite; zero for refocus)
    opt.RFr.alpha = [];     %   Actual tip angle distribution (degrees) (computed later...)

end

%   Define boundaries for various tissue compartments and fit models
opt.lsq.Icomp.X0   = [0.015 1e-3 1.00];      
opt.lsq.Icomp.XU   = [3.000 1e+1 1.60];      
opt.lsq.Icomp.XL   = [0.000 0.00 0.40];
%   Numeric fitting options
opt.lsq.fopt = optimset('lsqnonlin');
opt.lsq.fopt.TolX = 5e-4;     %   Fitting accuracy: 0.5 ms
opt.lsq.fopt.TolFun = 1.0e-9;
opt.lsq.fopt.MaxIter = 100;
opt.lsq.fopt.Display = 'off';




clear;clc;close;

addpath(genpath('D:\理综\重建\T2map_proj\T2map_proj2\data'))
addpath(genpath('bin'))

%% Loading data
data = '180NIST38791_12x';
data_all = load(data);
kspace = squeeze(data_all.(char(fieldnames(data_all))));

%% Basic Settings
% iteration parameters
opt = optset('n');
opt.th = 0;
opt.tau = 0.01;    % regularization parameter for wavelet transform 
opt.J = 3;         % low rank parameter
opt.niu = 0.1;     % low rank parameter for Hankel matrix


opt.T1=@(T2)4;
opt.esp = 10/1000;
opt.etl = 32;
opt.RFr.angle = 180;
opt.FA_array = [180/180,180/180,ones(1,30)];


opt.num_low_phase = 6;

% ESPIRIT parameters
opt.es.ksize = [max(ceil(opt.num_low_phase/4),2),max(ceil(opt.num_low_phase/4),2)];
opt.es.eigThresh_1 = 0.02;
opt.es.eigThresh_2 = 0.95;

%% pretreatment
if length(size(kspace)) ==4
    [opt.ny,opt.nx,opt.nt,opt.nc] = size(kspace);
else
    [opt.ny,opt.nx,opt.sl,opt.nt,opt.nc] = size(kspace);
    selection = 1;  
    kspace = squeeze(kspace(:,:,selection,:,:));
end
opt.mask = abs(squeeze(kspace(:,:,:,1)))>0;

% scaling  
scale = max(max(max(sqrt(sum(abs(ifftc(ifftc(kspace,2),1)).^2,4)))));
mat = kspace/scale;

% espirit
calib = crop(squeeze(sum(mat,3)/size(mat,3)),[opt.num_low_phase,opt.num_low_phase,opt.nc]);
[k,SE] = dat2Kernel(calib,opt.es.ksize);
idx = find(SE >= SE(1)*opt.es.eigThresh_1, 1, 'last' );
[ME,WE] = kernelEig(k(:,:,:,1:idx),[opt.ny,opt.nx]);
opt.maps = ME(:,:,:,end).*repmat(WE(:,:,end)>opt.es.eigThresh_2,[1,1,opt.nc]);

%% Iteration
row = opt.ny;
col = opt.nx;
echo = opt.nt;
channel = opt.nc;

I0 = reshape(PImerge(opt.maps,ifftc(ifftc(mat,2),1)),[row*col,echo]);

for i = 1:40
disp(['iteration:  ',num2str(i),'   '])
old = I0;
I0 = Sparsity_prior(I0,opt);
I0 = data_consistency(I0,mat,opt);
I0 = low_rank(I0,opt);
I0 = data_consistency(I0,mat,opt);
I0 = linear_pred(I0,opt);
I0 = data_consistency(I0,mat,opt);

disp(['compare:   ',num2str(mat_compare(I0,old))])
end

recon_results = reshape(I0*scale,[row,col,echo]);

%%   EPG fitting
% Z = abs(reshape(I0,[row*col,echo]));
% Z(abs(Z(:,1))<=opt.th,:) = 0;
% T2 = zeros(row*col,1);
% B1 = zeros(row*col,1);
% I = zeros(row*col,1);
% parfor j= 1:row*col
% [T2(j),B1(j),I(j)] = IniFit(squeeze(Z(j,:)),opt);
% end
% Tm = [T2,B1];
% 
% TIZ = [Tm,I*scale,Z*scale];
% save(fullfile('results',['prospective_',data,'_',num2str(opt.J),'_',num2str(opt.tau),'_',num2str(opt.niu),'.mat']),'TIZ')












%% Preliminaries
close all; clear all; clc
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','DOPBoxV1_7','DOPBoxV1-7','DOPbox'))
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis'))
c = 299.792458;
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
wn2T = @(wn) 1./((c*wn)./(1E7)); %wave number to period
T2wn = @(tau) 1./((c*tau)./(1E7));
%% Load the data, compress it, save it.
if (~exist('all_data_mat','var') || ~exist('T2vec','var'));
    load(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','data_stack.mat'));
end
normv = @(x) x/max(x);
mdata = squeeze(mean(all_data_mat,3));
p2ex = normv((mean(mdata,1)));
p2det = normv((mean(mdata,2)));
Wx = eye(size(mdata,2));%diag(p2ex);
Wy = eye(size(mdata,1));%diag(p2det);
[~,Wx_p] = chol(Wx); assert(Wx_p==0);
[~,Wy_p] = chol(Wy); assert(Wy_p==0);
Nbf = 60;
Trange = 1:185;
rDm = zeros(Nbf,Nbf,length(Trange));
[~,rDm(:,:,1),Bxw,Byw] = dop2DApproxWeighted(all_data_mat(:,:,1),Nbf,Nbf,Wx,Wy);
% for k=2:length(Trange);
%     [~,rDm(:,:,k),~,~] = dop2DApproxWeighted(all_data_mat(:,:,k),Nbf,Nbf,Wx,Wy);
% end
T = -T2vec(Trange);
% save(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','reduced_stack.mat'),'rDm','Bx','By','T');
imagesc(UncompressData(squeeze(rDm(:,:,1)),Bxw,Byw))
figure; imagesc((UncompressData(squeeze(rDm(:,:,1)),Bxw,Byw) - all_data_mat(:,:,1))./max(max(all_data_mat(:,:,1))))
%% Weird shit yo
x = mdata(300,:)';
bx = dop(length(x),100);
xr = zeros(672,100);
for k=1:100;
    xr(:,k) = bx*rand(100)*bx'*x;
end

clear all; close all; clc
% load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/PostProcessingData/' 'data_stack']);
load(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-07-06\PostProcessingData\' 'data_stack']);
c = 299.792458;
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
flatten_cube = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
flatten_matrix = @(m) reshape(m,[size(m,1)*size(m,2) 1]);

unflatten_matrix = @(m,X,Y,Z) reshape(m,[X Y Z]);
unflatten_vector = @(m,X,Y) reshape(m,[X Y]);
bin_pixels = @(m,L) squeeze(sum(reshape(m,[L size(m,1)/L size(m,2) size(m,3)]),1));
bin_axis = @(x,L) squeeze(mean(reshape(x,[L length(x)/L]),1)); %bin 4 consecutive pixels
data.data = bin_pixels(all_data_mat,5);
clear all_data_mat;
yaxis = bin_axis(y_axis_w_roi,5);
xaxis = x_axis_w_roi;

data.mean_data = squeeze(mean(data.data,3));
data.data_no_mean = data.data - repmat(squeeze(mean(data.data,3)),[1 1 size(data.data,3)]);
data.dataf = flatten_cube(data.data_no_mean);

%%

X = size(data.data_no_mean,1);
Y = size(data.data_no_mean,2);
f1 = tukeywin(X,0.1);
f2 = tukeywin(Y,0.1);
F = f1*f2';
N = 40;
Z = 1;%size(data.data_no_mean,3);
B = GenerateAssociatedLegendreBasis(X,Y,N,2);
Bf = flatten_cube(B);
Amat = zeros(N^2,Z);
D = flatten_matrix(data.data(:,:,1).*F);
for k=1:Z;
    Amat(:,k) = Bf\D;
end

imagesc(unflatten_vector(Bf*Amat(:,Z),X,Y)); figure; imagesc(unflatten_vector(D,X,Y));
figure; imagesc((unflatten_vector(D,X,Y) - unflatten_vector(Bf*Amat(:,Z),X,Y))./max(max(D)));
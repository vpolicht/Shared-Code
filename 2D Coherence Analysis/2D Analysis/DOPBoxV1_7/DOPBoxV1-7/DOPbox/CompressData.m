clear all; close all; clc;
load(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'data_stack']);
rDm = zeros(60,60,size(all_data_mat,3));
[~,rDm(:,:,1),Bx,By] = dop2DApprox(all_data_mat(:,:,1),60,60);
for k=2:size(all_data_mat,3);
    [~,rDm(:,:,k),~,~] = dop2DApprox(all_data_mat(:,:,k),60,60);
end
T = T2vec;
save(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'reduced_stack'],'rDm','Bx','By','T');
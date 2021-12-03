close all; clear all; clc;
load(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'reduced_bkg_stack']);
init_cube = UncompressData( squeeze(F(:,:,1)), Bx, By);
bkg_cube = zeros(size(init_cube,1),size(init_cube,2),length(T));
for k=1:length(T);
    bkg_cube(:,:,k) = UncompressData( squeeze(F(:,:,k)), Bx, By);
end
save(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'bkg_cube'],'bkg_cube','T');
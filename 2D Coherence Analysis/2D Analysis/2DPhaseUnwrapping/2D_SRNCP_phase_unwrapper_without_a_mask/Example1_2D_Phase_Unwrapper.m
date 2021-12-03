clear; clc; close all
%This Matlab program calls a 2D Phase unwrapper written in C language
%The wrapped phase image is floating point data type. 
%Also, the unwrapped phase image is floating data type.

img_name = 'complex_map_at_244.8129.mat';
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
P2image = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','Figures','Beat Analysis',img_name);
if(~exist('complex_map','var'))
    if(exist(P2image,'file'))
      load(P2image);                      %Load complex image
    end
end

%call the 2D phase unwrapper from C language
%To compile the C code: in Matlab Command Window type
%          mex Miguel_2D_unwrapper.cpp
%The wrapped phase should have the single (float in C) data type
WrappedPhase = single(angle(complex_map(200:900,100:500)));
UnwrappedPhase = Miguel_2D_unwrapper(WrappedPhase);
figure(2)
colormap('jet')
contourf(UnwrappedPhase); set(gca,'Ydir','normal');
hold all;
contour(abs(complex_map(200:900,100:500)),5,'-k');

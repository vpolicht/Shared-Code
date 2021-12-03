% Spectrometer calibration script
% June 2020 - V R Policht

clc; close all; clear all;
WhichOS; %Determines system syntax and base folder

%% Define spectrum to be checked against which calibration standard spectrum

calibspec = 'neon_calib_lamp.mat';

proj = 'TMD HS';
subproj = 'MWSe2';
date = '12-02-2021';

testspec = 'calib-760nm'; %Excel file saved directly from LabView GUI (.xlsx)
filext = 1; % 0 - .xlsx, 1 - .txt
%% Load files

load(calibspec);

sfold = [basefolder sla proj sla subproj sla date sla];
if ~filext
    tspec = xlsread([sfold testspec '.xlsx']);
else
    tspec = importdata([sfold testspec '.txt']);
    tspec = tspec.data(:,1:2);
end

if size(tspec,2) > 2
    tmp = tspec(:,1:2);
    clear tspec
    tspec = tmp;
end

nmr = @(x) (x-min(x))./max(x-min(x));
tspec(:,2) = nmr(tspec(:,2));

if length(unique(diff(tspec(:,1)))) > 1
    tmpwl = linspace(tspec(1,1),tspec(end,1),size(tspec,1));
    tmpsp = interp1(tspec(:,1),tspec(:,2),tmpwl);
    
    tspec(:,1) = tmpwl;
    tspec(:,2) = tmpsp;
end

if tspec(1,1) > tspec(end,1)
    tspec = flipud(tspec);
end


%%

figure(1); clf(1); hold all;
set(gcf,'position',[556 357 1363 638]);
subplot(2,1,1)
hold all
set(gca,'box','on')
plot(tspec(:,1),tspec(:,2),'linewidth',1.5);
plot(neonlamp(:,1),neonlamp(:,2),'linewidth',1)
legend('Test Spectrum','Standard')
xlim([500 800])
xlabel('Wavelength (nm)')

[pk loc] = findpeaks(tspec(:,2),tspec(:,1),'minpeakprominence',.01);


subplot(2,1,2)
hold all
plot(neonlamp(:,1),neonlamp(:,2),'linewidth',1.5);
set(gca,'box','on')
xlim([500 800])
xlabel('Wavelength (nm)')

set(gcf,'color','w')

match_peak_number = input('How many peaks would you like to use? ');

[testx testy stanx stany] = deal(zeros(match_peak_number, 1));

for m = 1:match_peak_number
    fprintf('Select peaks from Test Spectrum\n')
    clear x y
    subplot(2,1,1)
    [x, y] = ginput(1);
    xt = loc(dsearchn(loc,x));
    x = dsearchn(tspec(:,1),xt);
    
    plot(xt,y,'r*','handlevisibility','off')
    testx(m) = x;
    testy(m) = y;
    
    clear x y
    subplot(2,1,2)
    fprintf('Select corresponding peak in Standard spectrum\n')
    [x, y] = ginput(1);  
    xt = dsearchn(neonpeaks,x);
    x = neonpeaks(xt);
    plot(x,y,'r*','handlevisibility','off')
    stanx(m) = x;
    stany(m) = y;
end

pfit = polyfit(testx,stanx,2); %Fit with second order polynomial
pval = polyval(pfit,1:size(tspec,1));

subplot(2,1,2)
plot(pval,tspec(:,2),'linewidth',1);
legend('Standard Spectrum','Fit Calibration Spectrum')

figname = [testspec '_calibration-result'];
saveas(gcf,[sfold sla figname '.png'])
saveas(gcf,[sfold sla figname '.fig'])

%%
calwl.prewl = tspec(:,1);
calwl.wl = pval;

save([sfold sla testspec '.mat'],'calwl','pfit');
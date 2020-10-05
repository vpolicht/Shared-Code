% Look for spectral drift during an experiment (look at the scan by scan
% data) and then reaverage over a select part of the scans.
% Script currently requires read_bin.m and sgauss.m to run.
% V Policht, 10/2020

% This scipt requires a small amount of user input:
    % When prompted, you must select a "bin0" to define t1 = 0 fs
    % When prompted, you must select whether to window t1 with a Guassian
        % or tukey window function
    % When prompted, you must select a subplot where you will select which
        % scans to average over in creating a new .mat data file.
        
clear;
close all;
clc;

c = 299.792458; % Speed of light in nm/fs

%%

% Use these to define the folder where the data is and also what to name
% figures
% Where is the data
foldername = ['/Volumes/Utini/Data/V2O3/7-24-2020/'];
% Data set name
dataname = '2D_V2O3-spot1-100K-470flu-100ps-1';
 
% Calibration set name
calibname = '2DVIS-4_CALIB';
 
% Make an initial guess at the t1 axis
tlim = [-30 220]; % Motor t1 in fs

% Pump limits (nm)
pulm = [490 620];
% Probe limits (nm)
prlm = [490 620];

% Starting bin0
bin0 = 100;

% If you want to select data to reaverage: 0 - no, 1 - yes
rmman = 1;

%% Load and modify Calibration

% Load data
calib = load([foldername calibname '.dat']);

phint = calib(4:end,2);
datcal = calib(4:end,3:end);
wlcal = calib(1,3:end);
pumpspec = calib(2,3:end);
mstepcal = calib(4:end,1);

[srx sry] = size(datcal);
%%
figure(1); clf(1); hold all;
plot(wlcal,pumpspec);

figure(2); clf(2); hold all;
plot(mstepcal,phint)

figure(3); clf(3); hold all;
contour(wlcal,mstepcal,datcal,4,'color','k')

% Remove noisy ends of the t1 scan
t1lb = round(srx*.02);
t1ub = round(srx*.985);

mstlm = mstepcal(t1lb:t1ub);
phintlm = phint(t1lb:t1ub);
datcal = datcal(t1lb:t1ub,:);

% Reduce w3 axis based on threshold
thrsh = .007;
k = max(pumpspec)*thrsh;

w3lb = find(pumpspec > k,1);
w3ub = sry - find(fliplr(pumpspec) > k,1);

pumpspec = pumpspec(w3lb:w3ub);
wlcal = wlcal(w3lb:w3ub);
datcal = datcal(:,w3lb:w3ub);

datcal(isnan(datcal)) = 0;

%% Subtract the pump spectrum
datcal = datcal - mean(datcal,1);

figure(1);
plot(wlcal,pumpspec);

figure(2); 
plot(mstepcal,phint)

figure(3);
contour(wlcal,mstlm,datcal,4)

[srx sry] = size(datcal);

%% Perform w1 calibration

nfft = 2^14;

w3cal = (c*1e15)./wlcal; % in Hz

% Pseudofrequency axis
dx = mean(abs(diff(mstepcal)));
f_fft = ((1/dx)/2).*linspace(0,1,nfft/2+1);

ftdat = fft(datcal,nfft,1);
ftdatcal = abs(ftdat(1:nfft/2+1,:));
[a maxcorr] = max(ftdatcal);

tmp = abs(diff(diff(diff(maxcorr))));
mm = find(tmp<5);

b = maxcorr(mm);

% Fitting axis
ll = abs(mstlm(1) - mstlm(end));
pp = 1/(ll/srx*nfft);
psfr = maxcorr*pp.*1e15;

figure(1); clf(1); hold all;
contour(w3cal,f_fft,ftdatcal)
plot(w3cal(mm),f_fft(b),':')

p = polyfit(f_fft(b),w3cal(mm),1);

w1 = p(1)*f_fft + p(2);
wl1 = (c*1e15)./w1;


figure(2); clf(2); hold all;
contour(w3cal,w1,ftdatcal)
line(w3cal,w3cal,'color','r')
ylim([1e14 .8e15])

%% Select bin0 based on phase

good = 1;
unwrapd = 1;
srx = length(mstepcal);

padsize = 2^(nextpow2(srx)+1) - srx;
phintpad = padarray(phint,padsize,'post');

nfft = length(phintpad);
f_fft = ((1/dx)/2).*linspace(0,1,nfft/2+1);

w1 = p(1)*f_fft + p(2);
wl1 = (c*1e15)./w1;
w1lm = dsearchn(wl1',pulm')';

df = mean(abs(diff(w1)));
t1tmp = ((1e15/df)/2).*linspace(0,1,nfft/2+1);

while good == 1

    phintsft = circshift(phintpad,-bin0);

    fftphint = fft(phintsft,nfft);
    phipha = angle(fftphint);
    phiamp = abs(fftphint);

    close all;
    figure(1); clf(1); 
    subplot(2,1,1);
    hold all;
    plot(phintsft)
    title(['bin_0 = ' num2str(bin0)])
    set(gca,'box','on')
    xlabel('Motor Step (zero-padded)')

    subplot(2,1,2);
    yyaxis left
    plot(wl1,phiamp(1:nfft/2+1))
    ylabel('FFT amplitude')
    yyaxis right
    plot(wl1,phipha(1:nfft/2+1))
    ylabel('FFT phase')
    xlim([450 800]); xlabel('\lambda (nm)')

    set(gcf,'color','w');
    
    good = input('Change bin0? 0 - no, 1 - yes\n');
    
    if unwrapd
        phipha = unwrap(phipha);
    end
    
    phipha = flipud(phipha);
    expha = exp(1i.*phipha);
    
    if good 
        bin0 = input('New bin0 = \n');
    end
end

% Find maximum of circleshifted phase
x0 = find(abs(phintsft) >= max(abs(phintsft)));
%% Pick a filter type for along t1

datcal = calib(4:end,3:end);
datcal = datcal(bin0:end,w3lb:w3ub);
datcal(isnan(datcal)) = 0;
datcal = datcal - mean(datcal,1);

SRX = length(bin0:length(mstepcal));
nmr = @(a) a./max(a);

figure(2); clf(2); hold all;
plot(1:SRX,nmr(sum(datcal,2)))

filtype = input('Choose t1 filter: 0 - Gaussian, 1 - Tukey\n');

filt = zeros(2*SRX,1);
tmpax = 1:2*SRX;
fact = .8;

if ~filtype
    filt = gaussmf(tmpax,[fact*SRX/2 SRX+x0-1]);
else
    filt = tukeywin(2*SRX,fact);
end

filt = filt(SRX+1:end);
filt = nmr(filt);

plot(1:SRX,filt)
plot(1:SRX,nmr(sum(datcal.*filt,2)))

dcalfilt = datcal.*filt;
%%

filtpd = zeros(2*srx,1);
tmpax = 1:2*srx;
fact = .8;

if ~filtype
    filtpd = gaussmf(tmpax,[fact*srx/2 srx+x0-1]);
else
    filtpd = tukeywin(2*srx,fact);
end

filtpd = filtpd(srx+1:end);
filtpd = nmr(filt);

ftphi = fft(-1.*phint(bin0:end).*filtpd,nfft);
ftphi = abs(ftphi(1:nfft/2+1));

figure(1); clf(1); hold all;
plot(wl1,ftphi./max(ftphi))
plot(wlcal,pumpspec./max(pumpspec))
xlim(pulm)

ftdcal = fft(dcalfilt,nfft);
ftdcal = abs(ftdcal(1:nfft/2+1,:));

figure(2); clf(2); hold all
tmp = sum(ftdcal,2);
plot(wl1,tmp./max(tmp))
plot(wlcal,pumpspec./max(pumpspec))
plot(wl1,ftphi./max(ftphi))
legend('Calibrated FT','Pump Spectrum','PD FT')
xlim([400 1000]); xlabel('Wavelength (nm)')
set(gca,'box','on')
set(gcf,'color','w')
title(calibname)

saveas(gcf,[foldername 'calibrated_' calibname '.fig'])
saveas(gcf,[foldername 'calibrated_' calibname '.png'])

%% Load in the scan by scan 2D data

% *Define folder path
datfolder = [foldername 'temp/'];


files = dir(datfolder);
files(1:3) = [];

% Look for all files matching the 2D data name
count = 1;
for m = 1:length(files)
    
    if contains(files(m).name,dataname) && contains(files(m).name,'.bin')
        bfile(count) = files(m);
        count = count+1;
    end
    
end

clear files 
%% Define the t2 axis and the number of scans 

for m = 1:length(bfile)
    clear tmp a
    
    tmp = bfile(m).name;
    a = strfind(tmp,'_','ForceCellOutput',true);
    time(m) = str2double(tmp(a{1}(2)+1:a{1}(3)-3));
    avnum(m) = str2double(tmp(a{1}(3)+1:end-4));
   
end

time = unique(time);
avnum = unique(avnum);
avnum = avnum(1:end-2);

nt2 = length(time);
nav = length(avnum);

%%
% Pick a t2 value to load in
tvalue = 100;
tind = dsearchn(time',tvalue);

% Read in one file to define the motor step and probe wavelength
temp = read_bin([datfolder bfile(1).name]);
mstep = temp(4:end,1);
wpr = temp(1,3:end);
srx = length(mstep);
sry = length(wpr);

clear temp dat phint dtm cp hms mn

phint = zeros(nt2,nav,length(mstep)); % Define pump-pump interferogram

for m = tind
    parfor n = 1:nav
        
        filename = [dataname '_' num2str(time(m)) 'fs_' num2str(avnum(n)) '.bin'];
        temp = read_bin([datfolder filename]);
        
        phint(m,n,:) = temp(4:end,2); % Pump-pump interferogram
       
    end
end

%% Look at pump-pump interferogram as a function of average number

for M = tind
    tindc = M;
    
    figure(1); clf(1);
    hold all;
    imagesc(mstep,avnum,squeeze(phint(tindc,:,:)))
    set(gca,'box','on','layer','top')
    set(gcf,'color','w')
    xlabel('Motor Step')
    ylabel('Scan number')
    title(['t_2 = ' num2str(time(M)) ' fs'])
    xlim([1 mstep(end)])
    ylim([1 nav-1])
    
    figure(2); clf(2);
    hold all;
    plot(mstep,mean(squeeze(phint(tindc,:,:)),1))
    set(gca,'box','on','layer','top')
    set(gcf,'color','w')
    xlabel('Motor Step')
    ylabel('PD intensity')
    title(['t_2 = ' num2str(time(M)) ' fs'])

    % Fourier transform over motor steps
    ftint = zeros(nt2,nav,nfft/2+1);

    for m = 1:nt2
        parfor n = 1:nav
            temp = fft(squeeze(phint(m,n,:)),nfft);
            ftint(m,n,:) = abs(temp(1:nfft/2+1));
        end
    end

    %
    figure(3); clf(3); hold all;
    contour(wl1,avnum,squeeze(ftint(tindc,:,:)),20,'linewidth',2)
    xlim(pulm)
    set(gca,'box','on','layer','top')
    set(gcf,'color','w')
    xlabel('Pump Wavelength')
    ylabel('Scan number')
    title(['t_2 = ' num2str(time(M)) ' fs'])

    %
    figure(4); clf(4); hold all;
    plot(wl1,ftphi./max(ftphi),'b','linewidth',1.5)
    plot(wlcal,pumpspec./max(pumpspec),'k','linewidth',1.5)
   
    clear tm2
    for m = 1:nav
        clear temp
        temp = squeeze(ftint(tindc,m,:));
%         plot(cal.wl1,temp./max(temp)) 
        tm2(m,:) = temp./max(temp);
    end
    
    mn = mean(tm2);
    sd = std(tm2);
    errorbar(wl1,mn,sd,'linewidth',1.2)
    
    legend('FT PPI from Calibration','Wavelength resolved Pump',...
        'Mean FFT of PPI during Experiment')
    
    set(gca,'box','on','layer','top')
    set(gcf,'color','w')
    xlabel('Pump Wavelength')
    ylabel('FFT Intensity')    
    title(['t_2 = ' num2str(time(M)) ' fs'])    
    xlim(pulm)
 
end

%% Determine expertimental time in lab frame

for n = 1:nav
    filename = [dataname '_' num2str(time(end)) 'fs_' num2str(avnum(n)) '.bin'];
    temp = dir([datfolder filename]);
    dtm{n} = datestr(temp.datenum,13);
    cp = strfind(dtm{n},':');
    hms(n,:) = [str2num(dtm{n}(1:cp(1)-1)) str2num(dtm{n}(cp(1)+1:cp(2)-1))...
     str2num(dtm{n}(cp(2)+1:end))];
    mn(n) = hms(n,1)*60 + hms(n,2) + hms(n,3)/60;
    
    labtime(n) = hms(n,1) + hms(n,2)/60 + hms(n,3)/(60*60);
end

actime = mean(abs(diff(mn))); % Time between scans of the same t2 in min
midnight = dsearchn(labtime',0);
labtime(midnight:end) = labtime(midnight:end) + 24;

%% Look at the 2D signal as a function of average number for a specific t2

dat = zeros(nav,srx*sry);

parfor n = 1:nav
    filename = [dataname '_' num2str(time(tind)) 'fs_' num2str(avnum(n)) '.bin'];
    temp = read_bin([datfolder  filename]);

    dat(n,:) = reshape(temp(4:end,3:end),srx*sry,1);
end

%%

% Wavelength limits for the probe axis
ylb = dsearchn(wpr',prlm(1));
yub = dsearchn(wpr',prlm(end));

% Wavelength limits for the FFT pump axis
xub = dsearchn(wl1',pulm(1));
xlb = dsearchn(wl1',pulm(end));
w1 = (c*1e15)./wl1;
dw = mean(abs(diff(w1)));

% Define variables
clear ftmp pp mx v1l3 corr ftcorr tmp3
ftmp = zeros(nav,length(xlb:xub),length(ylb:yub));
t1l3 = zeros(nav,srx,length(ylb:yub));
corr = zeros(nav,length(xlb:xub));
ftcorr = zeros(nav,length(xlb:xub),length(ylb:yub));

for m = 1:nav
    clear tmp tmp2 tmp3 tmp4 datasub t1filt B

    tmp = reshape(dat(m,:),srx,sry);
    t1l3(m,:,:) = tmp(:,ylb:yub); 

    tmpift = fft(squeeze(t1l3(m,:,:)),nfft,1);

    B = pulm(1) + (pulm(end)-pulm(1))/2;
    B = (c*1e15)/B;
    t1filt = sgauss(w1,75*dw,B,6);
    t1filt = cat(2,t1filt,fliplr(t1filt(2:end-1)))';

    datasub = ifft(tmpift.*t1filt,nfft,1);   
    t1l3(m,:,:) = datasub(1:srx,:);

    tmp4 = fft(t1l3(m,:,:),nfft,2); % FFT 2D(nav,t1,l3) --> 2D(nav,l1,l3)
    ftcorr(m,:,:) = abs(tmp4(:,xlb:xub,:)); % Keep only some l1
    corr(m,:) = sum(ftcorr(m,:,:),3); % Sum over probe wl

    pp(m,:) = sum(squeeze(ftcorr(m,:,:)),1); % Sum over the pump wl
    mx(m) = max(max(squeeze(ftcorr(m,:,:)))); % Find the maximum for each scan
end

%%
tmphi = mean(squeeze(phint(tind,:,:)),1); % Look at PPI for t2 of interest
dcorr = sum(mean(t1l3,1),3); % Mean of the scans, sum over probe wl

figure(2); clf(2);
hold all;
plot(mstep,tmphi./max(tmphi))
plot(mstep,dcorr./max(dcorr))
set(gca,'box','on','layer','top')
set(gcf,'color','w')
xlabel('Motor Step')
ylabel('PD intensity')
title(['t_2 = ' num2str(time(tind)) ' fs'])
% legend({['PPI for t_2 = ' num2str(time(tind)) ' fs']},{'FFT Data, Mean over n_{av}, sum over \lambda_3'})

sumephi = sum(squeeze(ftint(tind,:,:)),1); % Select t2 and sum over pump wl
mnphi = mean(squeeze(ftint(tind,:,3:10:end)),1); % take mean
stphi = std(squeeze(ftint(tind,:,3:10:end)),1); % take STD

sumcorr = sum(corr,1); % Sum over scans
mncorr = mean(corr(:,1:10:end),1); % Take mean
mncorr = mncorr - min(mncorr); % Offset mean
stcorr = std(corr(:,1:10:end),1); % Take STD

%%

% Define more pump wl limits
ub = dsearchn(wl1',pulm(1));
lb = dsearchn(wl1',pulm(end));
ub = length(wl1);
lb = 1;
ww1 = wl1(lb:ub);
skipgif = 0; % If you don't want to make the gif, enter 1

% Define where to save the figures
figfold = [foldername dataname '/Figures/'];

% Make if doesn't exist
if ~exist(figfold,'dir')
    mkdir(figfold)
end

% Start writing the .gif and .mp4
if ~skipgif
    gifname = [figfold dataname '_scannumber.gif'];
    dt = .1;

    writerObj = VideoWriter([figfold dataname '_scannumber.mp4'],'MPEG-4');
    writerObj.FrameRate = 5;
    open(writerObj);
end

% Iterate over number of scans (nav)
for M = 1:nav
    clear ppi corft pfit pfev cfit cfev

    ppi = squeeze(ftint(tind,M,lb:ub));
    ppi = ppi./max(ppi);

    corft = smoothdata((corr(M,:)-min(corr(M,:)))./max(corr(M,:)-min(corr(M,:))),'movmean',1);

    % Calculate PPI width
    [a b] = max(ppi);
    pedge = [ww1(dsearchn(ppi(1:b),.5)) ww1(b+dsearchn(ppi(b+1:end),.5))];
    pwid(M) = abs(diff(pedge));
    px = [pedge(2)+5;pedge(1)+5];
    py = .5;

    [a b] = max(corft);
    cedge = [ww1(xlb-1+dsearchn(corft(1:b)',.5)) ww1(xlb-1+b+dsearchn(corft(b:end)',.5))];
    cwid(M) = abs(diff(cedge));
    cx = [cedge(2)-40;cedge(1)+5];
    cy = .5;
    
    if ~skipgif
        figure(3); clf(3); hold all;
        plot(wl1,ftphi./max(ftphi),'k','linewidth',1.2)

        ax = gca;
        ax.ColorOrderIndex = 1;

        errorbar(wl1(3:10:end),mnphi./max(mnphi),stphi./max(mnphi),'linewidth',1)
        errorbar(wl1(xlb:10:xub),mncorr./max(mncorr),stcorr./max(mncorr),'linewidth',1)

        ax.ColorOrderIndex = 4;
        plot(ww1,ppi,'linewidth',2)
        text(px(2),py,['\Delta\lambda_1 = ' num2str(round(pwid(M))) ' nm'],...
            'fontweight','bold','fontsize',12,'Color',ax.ColorOrder(4,:));

        ax.ColorOrderIndex = 7;
        plot(ww1(xlb:xub),corft,'linewidth',2)
        text(cx(1),cy,['\Delta\lambda_1 = ' num2str(round(cwid(M))) ' nm'],...
            'fontweight','bold','fontsize',12,'Color',ax.ColorOrder(7,:));


        xlim([pulm(1)-50 pulm(end)+50])
        ylim([-.1 1.1])
        set(gca,'box','on','layer','top','linewidth',1.5)
        set(gcf,'color','w')
        legend('Calibration PPI FT',...
            'Experiment PPI Stats','Experiment \omega_1 stats',...
            'Single Experiment PPI','Single Experiment \omega_1')
        title(['t_2 = ' num2str(time(tind)) '; Scan number ' num2str(M)])
        xlabel('Pump Wavelength')
        ylabel('Int (a.u.)')
        
        % Write this image to the .gif and .mp4
        drawnow
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);

        if M == 1 % For the first frame
           imwrite(imind,cm,gifname,'gif','Loopcount',inf,'DelayTime',dt);
        else % for every other frame
           imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',dt);
        end
        writeVideo(writerObj,frame);
        clear frame

    end

end
close(writerObj); % NEED to close the .mp4 object or it will not save

%%
pout = ~isoutlier(pwid);

pwid = pwid(pout);
lin = [1:nav];
lin = lin(pout);
pfit = fit(lin',pwid','poly1');
pfev = feval(pfit,1:nav);

cfit = fit([1:nav]',cwid','poly1');
cfev = feval(cfit,1:nav);

figure(4); clf(4); hold all;
plot(lin,pwid,'ro')
plot(1:nav,pfev,'r','linewidth',1.2,'handlevisibility','off')
plot(1:nav,cwid,'mx')
plot(1:nav,cfev,'m','linewidth',1.2,'handlevisibility','off')

set(gca,'box','on','layer','top','linewidth',1.5)
set(gcf,'color','w')
xlabel('Scan Number')
ylabel('FFT FWHM (nm)')
title(['t_2 = 100 ps'])
legend('PPI','2D Data')

figname = [figfold dataname '_fftwidths-100ps'];
print([figname '.png'],'-dpng');

%%
figure(1); clf(1); hold all;
plot(mx)
title('Signal Maximum as a function of scan number')

figure(2); clf(2); hold all;
imagesc(wpr(ylb:yub),1:nav,pp)
title('"Pump-Probe" signal as a function of scan number')

figure(3); clf(3); hold all;
plot(1:nav,sum(pp,2))
title('Total "Pump Probe" signal as a function of scan number')

%% Reaverage the data over selected scans

close all;

ppa = sum(squeeze(ftint(tind,:,:)),2);
nvec = 1:nav;
pbckg = smoothdata(ppa,'movmean',floor(nav/10));
ppb = ppa - pbckg;
pout = isoutlier(ppb,'threshold',2);

figure(1); clf(1); hold all;
plot(nvec,ppb)
plot(nvec(pout),ppb(pout),'o')
set(gca,'box','on')
set(gcf,'color','w')
xlabel('N_{scans}')
legend('Background subtracted average PP signal','Outlier scans','location','northoutside')

reav = nvec(~pout);
nvec = nvec(reav);

if rmman
    %%
    wl1lm = dsearchn(wl1',[pulm(end)+100 pulm(1)-100]');
    pwl1 = wl1(wl1lm(1):wl1lm(end));
    
    ppitr = squeeze(ftint(tind,reav,wl1lm(1):wl1lm(2)));
    [ma ml] = max(ppitr,[],2);
    mloc = pwl1(ml);
   
    
    [lb ub wid] = deal(zeros(1,length(nvec)));
    for n = 1:length(nvec)
        clear tmpp
        tmpp = ppitr(n,:);
        tmpp = tmpp./max(tmpp);
        lb(n) = pwl1(ml(n)+1+find(tmpp(ml(n)+1:end) < 0.57 & tmpp(ml(n)+1:end) > 0.43,1));
        ub(n) = pwl1(find(tmpp(1:ml(n)) < 0.57 & tmpp(1:ml(n)) > 0.43,1,'last'));
        wid(n) = (ub(n)) - (lb(n));
    end
    
    
    figure(1); clf(1); hold all;
    contour(pwl1,nvec,ppitr,20,'linewidth',2,'handlevisibility','off')
    
    figure(2); clf(2); hold all;
    contour(pwl1,nvec,ppitr,20,'linewidth',2,'handlevisibility','off')
    colormap('gray')
    plot(mloc,nvec,'b')
    plot(lb,nvec,'r')
    plot(ub,nvec,'r')
    legend('Peak','HM Lower Bound','HM Upper Bound','location','northoutside','Orientation','horizontal')
    set(gca,'box','on','layer','top','xgrid','on','xminorgrid','on')
    set(gcf,'color','w')
    xlabel('Pump Wavelength')
    ylabel('Scan number')
    title(['Pump-Pump Interferogram for t_2 = ' num2str(time(tind)) ' fs'])

    figure(3); clf(3); set(gcf,'color','w')
    subplot(2,1,1); hold all;
    plot(nvec,mloc,'x')
    set(gca,'box','on','ygrid','on','yminorgrid','on')
    title('FFT PPI (Pump spectrum) Peak location')
    
    subplot(2,1,2); hold all;
    plot(nvec,wid,'x')
    set(gca,'box','on','ygrid','on','yminorgrid','on')
    title('FFT PPI (Pump spectrum) FWHM')
  
    num = input('Select subplot to select scans from: 1 - Peak location, 2 - FWHM\n');
    if num == 1
        subplot(2,1,1)
    end
    
    for n = 1:4
       
        fprintf('Draw rectangle (4 points) around the scans you want to keep\n')
        [a(n) b(n)] = ginput(1);
        plot(a(n),b(n),'k.')
        
        if n > 1
            line([a(n) a(n-1)],[b(n) b(n-1)],'color','k')
        end
        
        if n == 4
            line([a(1) a(n)],[b(1) b(n)],'color','k')
        end
    end
        
    slm = [min(round(a)) max(round(a))];
    wlm = [min(round(b)) max(round(b))];
    
    if num == 1
        inds = find(mloc <= wlm(end) & mloc >= wlm(1));
    else
        inds = find(wid <= wlm(end) & wid >= wlm(1));
    end
    
    nvec = nvec(inds);    
    
end

%% Reaverage using new indecies
avdat = {length(time)};
avphi = {length(time)};

for m = 1:length(time)

    ttx = time(m);
    
    data = zeros(srx,sry);
    phint = zeros(srx,1);

    parfor n = 1:length(nvec)

        filename = [dataname '_' num2str(ttx) 'fs_' num2str(nvec(n)) '.bin'];
        temp = read_bin([datfolder sla filename]);

        data = squeeze(temp(4:end,3:end)) + data;
        
        phint = temp(4:end,2) + phint;
        
    end

    avdat{m} = data./length(nvec);
    avphi{m} = phint./length(nvec);
end

%%
matfile = [foldername dataname '.mat'];
matsavefile = [foldername dataname '_av' num2str(length(nvec)) '.mat'];

if exist(matfile,'file')
    load(matfile);
    clear NumScans Data PhInt
    
    NumScans = length(nvec);
    Data = avdat;
    PhInt = avphi;
    
    save(matsavefile,'Data','PhInt','T2','NumScans','wpr','mStep')
else
    NumScans = length(nvec);
    Data = avdat;
    T2 = time;
    PhInt = avphi;
    
    save(matsavefile,'Data','PhInt','T2','NumScans','wpr','mStep');
end

    
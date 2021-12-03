% Look for spectral drift during an experiment 
% (Look at the scan by scan data)

clear;
close all;
clc;

WhichOS;
set(0,'defaulttextInterpreter','latex');
%% Define input parameters

% Use these to define the folder where the data is & any important
% parameters
% VP - I organize my data by project/subproject/date/ 
% Modify as necessary and check that 'foldname' locates your data folder
% properly
proj = 'Perovskites';
subproj = '';
twist = '';
temp = 'RT';
date = '11-10-2021';
tdname = '2D_1-pervTest';

% Which pump calibration file
calibname = '2DVIS-1_CALIB';
 
% Make an initial guess at the t1 axis
tlim = [-30 220]; % Motor t1 in fs
singlet2 = '';

% Pump limits (nm)
pulm = [500 700];
% Probe limits (nm)
prlm = [500 700];
%% Load and modify Calibration

foldname = [basefolder sla proj sla subproj sla date sla];
% Load data
calib = load([foldname calibname '.dat']);

phint = calib(4:end,2);
datcal = calib(4:end,3:end);
wlcal = calib(1,3:end);
pumpspec = calib(2,3:end);
mstepcal = calib(4:end,1);

[srx sry] = size(datcal);
%% Plot Spectra
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

% Reduce w3 axis
w3lb = dsearchn(wlcal',pulm(1));
w3ub = dsearchn(wlcal',pulm(2));

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
line(w3cal,w3cal,'color','r');
ylim([1e14 .8e15])
%% Pick a filter type for along t1

SRX = length(mstlm);%length(bin0:length(mstepcal));
nmr = @(a) a./max(a);

figure(2); clf(2); hold all;
plot(1:SRX,nmr(sum(datcal,2)))

filtype = input('Choose t1 filter: 0 - Gaussian, 1 - Tukey\n');

tmpax = 1:length(datcal);
fact = .7;

[filt] = deal(zeros(2*SRX,1));

if ~filtype
    [a,x0] = max(abs(sum(datcal(1:SRX/2,:),2)));
    filt = (gaussmf(tmpax,[fact*SRX/2 x0]))';
    filt = filt(1:SRX);
    filt = filt - min(filt);
else
    filt = tukeywin(2*SRX,fact);
    filt = filtcirc(SRX+1:end);
end

filt = nmr(filt);

plot(1:SRX,filt)
plot(1:SRX,nmr(sum(datcal.*filt,2)))

datcalf = datcal.*filt;
phintf = phintlm.*filt;
%% Load in some bin0 finding parameters
srx = length(phintf);

nfft = 2^(nextpow2(srx)+1);

f_fft = ((1/dx)/2).*linspace(0,1,nfft/2+1);

w1 = p(1)*f_fft + p(2);
wl1 = (c*1e15)./w1;
w1lm = dsearchn(wl1',pulm')';

df = mean(abs(diff(w1)));
t1tmp = ((1e15/df)/2).*linspace(0,1,nfft/2+1);

[a mx] = max(abs(phintf));

padsize = nfft - srx;
phintpad = padarray(phintf,padsize,'post');
datpad = padarray(datcalf,padsize,'post');
%% New bin0 seeking method
[a,x0] = max((phintpad));

x0wind = 25;
yax = x0-x0wind:x0+x0wind;

wl3strong = wlcal(find(pumpspec>=max(pumpspec)*.4));
wslm = dsearchn(wl1',[wl3strong(end) wl3strong(1)]');
wax = wl1(wslm(1):wslm(2));

ftv = zeros(x0wind*2+1,nfft);
[fsnd off fslope] = deal(zeros(1,length(yax)));
fitz = zeros(length(yax),length(wax));
immin = zeros(length(yax),1);

for m = 1:length(yax)
    M = yax(m);
    ftv(m,:) = fft(circshift(phintpad,-M),nfft);
    
    ff = fit(wax',smoothdata((angle(ftv(m,wslm(1):wslm(2)))),10)','poly1');
%     fsnd(m) = ff.p1;
    fslope(m) = ff.p1;
    foff(m) = ff.p2;
    fitz(m,:) = feval(ff,wax);
    
    immin(m) = min(sum(imag(ftv(m,wslm(1):wslm(2))),2));
end

figure(5); clf(5); hold all;
yyaxis left
plot(yax,phintf(yax))
yyaxis right
hold all;
plot(yax,abs(fslope./max(fslope)),'o')
plot(yax,foff,'rd')
% plot(yax,fsnd./max(fsnd),'g+')
plot(yax,abs(immin./max(immin)),'kx')
ylim([-1 1])

[a b] = min(abs(fslope));
bin0 = yax(b);
xline(bin0);
%% Select bin0 based on phase

good = 1;
unwrapd = 1;

while good == 1

    phintsf = circshift(phintpad,-bin0);
    datpsf = circshift(datpad,-bin0);
    
    fftphint = fft(phintsf./max(phintsf),nfft);
    fftdat = fft(-1.*sum(datpsf,2)./max(sum(datpsf,2)),nfft,1);
    
    phipha = angle(fftphint);
    phiamp = abs(fftphint);

    datpha = angle(fftdat);
    datamp = abs(fftdat);

    figure(1); clf(1); 
    subplot(2,1,1);
    hold all;
    yline(0,'handlevisibility','off');
    plot(phintlm./max(phintlm))
    plot(-1.*sum(datcal,2)./max(sum(datcal,2)))
    xline(bin0,'handlevisibility','off');
    xlim([mx - 2e2 mx + 2e2])
    title(['$bin_0$ = ' num2str(bin0) ', filt w. ' num2str(fact)])
    set(gca,'box','on','XTickLabel','')
    xlabel('Motor Step (a.u.)')

    subplot(2,1,2);
    yyaxis left
    hold all;
    plot(wl1,phiamp(1:nfft/2+1)./max(phiamp(1:nfft/2+1)))
    plot(wl1,datamp(1:nfft/2+1)./max(datamp(1:nfft/2+1)))
    plot(wlcal,pumpspec./max(pumpspec));
    ylabel('FFT amplitude')
    yyaxis right
    hold all;
    plot(wl1,phipha(1:nfft/2+1))
    plot(wl1,sum(datpha(1:nfft/2+1,:),2));
    yline(0,'handlevisibility','off');
    ylabel('FFT phase')
    xlim([450 900]); xlabel('$\lambda$ (nm)')
    set(gcf,'color','w');
    set(gca,'box','on')
    
    good = input('Change bin0? 0 - no, 1 - yes\n');
    
    if unwrapd
        phipha = unwrap(phipha);
    end
    
    phipha = flipud(phipha);
    expha = exp(1i.*phipha);
    
    if good 
        bin0 = input('New bin0 = \n');
        clear phintmp
    end
end

saveas(gcf,[basefolder sla proj sla subproj sla date sla 'bin0_' calibname '.fig'])
saveas(gcf,[basefolder sla proj sla subproj sla date sla 'bin0_' calibname '.png'])

% Find maximum of circleshifted phase
x0 = find(abs(phintsf) >= max(abs(phintsf)));

datamp = datamp(1:nfft/2+1);
phiamp = phiamp(1:nfft/2+1);
%% Plot comparison of spectra amplitudes 

figure(2); clf(2); hold all
plot(1240./wl1,datamp./max(datamp))
plot(1240./wl1,phiamp./max(phiamp))
plot(1240./wlcal,pumpspec./max(pumpspec))
legend('Calibrated FT','PD FT','Pump Spectrum','PD CS FT')
xlim([1240./[1000 400]]); xlabel('$\hbar\omega$ (eV)')
set(gca,'box','on')
set(gcf,'color','w')
title(calibname)

saveas(gcf,[basefolder sla proj sla subproj sla date sla 'calibrated_' calibname '.fig'])
saveas(gcf,[basefolder sla proj sla subproj sla date sla 'calibrated_' calibname '.png'])
%% Load in the scan by scan 2D data 

% *Define folder path
if ~isempty(subproj)
    datfolder = [basefolder sla proj sla subproj sla date sla 'temp'];
else
    datfolder = [basefolder sla proj sla date sla 'temp'];
end

files = dir(datfolder);

% Look for all files matching the 2D data name
count = 1;
count2 = 1;
for m = 1:length(files)
    
    if contains(files(m).name,'(2)')
        delete([datfolder sla files(m).name]);
    elseif contains(files(m).name,tdname) && contains(files(m).name,'.bin') && ~contains(files(m).name,'ppscatter')
        bfile(count) = files(m);
        count = count+1;
    elseif contains(files(m).name,tdname) && contains(files(m).name,'.bin') && contains(files(m).name,'ppscatter')
        sfile(count2) = files(m);
        count2 = count2 + 1;
    end
    
end

clear files 
%% Define the t2 axis and the number of scans 
clear ltime time avnum 

numund = strfind(tdname,'_','ForceCellOutput',true);
vsts = length(numund{1});

for m = 1:length(bfile)
    clear tmp a
    
    tmp = bfile(m).name;
    ltmp = bfile(m).date;
    ltime(m) = datetime(ltmp,'Format','dd-MMM-yyyy HH:mm:ss');
    
    a = strfind(tmp,'_','ForceCellOutput',true);
    time(m) = str2double(tmp(a{1}(vsts+1)+1:a{1}(vsts+2)-3));
    avnum(m) = str2double(tmp(a{1}(vsts+2)+1:end-4));

    if m == 1
        a = read_bin([datfolder sla tmp]);
        cohtime = length(a(4:end,1));
    end
    
end

ltime = sort(ltime);
time = unique(time);
avnum = unique(avnum);
avnum = avnum(1:end-2);

nt2 = length(time);
nav = length(avnum);

expend = max(ltime);
expsts = min(ltime);
explength = expend-expsts;
%% Read in data 
if input('Analyze all t2 - 1 or a single t2 - 0\n')
    tvalues = time;
else
    % Pick a t2 value to load in
    if singlet2
        tvalues = input('select t2 point\n');
    else
        tvalues = 300;
    end
end

phint = zeros(length(tvalues),nav,cohtime);
tic
for m = 1:length(tvalues)
    tvalue = tvalues(m);
    tind = dsearchn(time',tvalue);
    tval = time(tind);
    
    if m == 1
        % Read in one file to define the motor step and probe wavelength
        temp = read_bin([datfolder sla bfile(1).name]);
        mstep = temp(4:end,1);
        wpr = temp(1,3:end);
        srx = length(mstep);
        sry = length(wpr);

        clear temp dat dtm cp hms mn
    end
    
    parfor n = 1:nav

        filename = [tdname '_' num2str(tval) 'fs_' num2str(avnum(n)) '.bin'];
        temp = read_bin([datfolder sla filename]);

        phint(m,n,:) = temp(4:end,2);

    end
end
toc

%% Read in Scatter data, if exists

if exist('sfile','var')
    for m = 1:nav
        if m == 1

            temp = read_bin([datfolder sla sfile(m).name]);
            smstep = temp(4:end,1);
            swpr = temp(1,3:end);
            sdat = zeros(size(temp(4:end,3:end)));
            scat = zeros([nav,size(temp(4:end,3:end))]);
            
            sphint = zeros(size(temp(4:end,2)));
            
            clear temp
        end
        
        filename = [tdname '_0fs_ppscatterwithdiv_' num2str(avnum(m)) '.bin'];
        temp = read_bin([datfolder sla filename]);
        
        scat(m,:,:) = temp(4:end,3:end);
        sdat = squeeze(temp(4:end,3:end)) + sdat;
        sphint = temp(4:end,2) + sphint;
    
    end
    
    sdat = sdat./nav;
    sphint = sphint./nav;
    
end

%% Look at pump-pump interferogram as a function of average number 

[ftint t2avwp] = deal(zeros(length(tvalues),nav,floor(nfft/2+1)));
showfig = 0;

for m = 1:length(tvalues)
    if length(tvalues) == 1
        tindc = 1;
    else
        tindc = dsearchn(time',tvalues(m));
    end
    
    if showfig
        figure(1); clf(1);
        hold all;
        imagesc(mstep,avnum,squeeze(phint(tindc,:,:)))
        set(gca,'box','on','layer','top')
        set(gcf,'color','w')
        xlabel('Motor Step')
        ylabel('Scan number')
        title(['t_2 = ' num2str(time(m)) ' fs'])
        xlim([1 mstep(end)])
        ylim([1 nav-1])

        figure(2); clf(2);
        hold all;
        plot(mstep,mean(squeeze(phint(tindc,:,:)),1))
        set(gca,'box','on','layer','top')
        set(gcf,'color','w')
        xlabel('Motor Step')
        ylabel('PD intensity')
        title(['t_2 = ' num2str(time(m)) ' fs'])
    end
    
    % Fourier transform over motor steps
    parfor n = 1:nav
        tmp = squeeze(phint(m,n,:));
        tmp = tmp(t1lb:t1ub);
        tmpf = tmp.*filt;
        tmfpd = padarray(tmpf,padsize,'post');
        tmfpd = circshift(tmfpd,-bin0);
        temp = fft(tmfpd,nfft);
        ftint(m,n,:) = abs(temp(1:nfft/2+1));
    end

    if showfig
        figure(3); clf(3); hold all;
        contour(wl1,avnum,squeeze(ftint(tindc,:,:)),20,'linewidth',2)
        xlim(pulm)
        set(gca,'box','on','layer','top')
        set(gcf,'color','w')
        xlabel('Pump Wavelength')
        ylabel('Scan number')
        title(['t_2 = ' num2str(time(m)) ' fs'])

        %
        figure(4); clf(4); hold all;
        plot(wl1,ftphi./max(ftphi),'b','linewidth',1.5)
        plot(wlcal,pumpspec./max(pumpspec),'k','linewidth',1.5)
    end
    
    clear tm2
    for M = 1:nav
        clear temp
        temp = squeeze(ftint(tindc,M,:));
%         plot(cal.wl1,temp./max(temp)) 
        tm2(M,:) = temp./max(temp);
    end
    
    t2avwp(m,:,:) = tm2;
    mn = mean(tm2);
    sd = std(tm2);
    
    if showfig
        errorbar(wl1,mn,sd,'linewidth',1.2)
        legend('FT PPI from Calibration','Wavelength resolved Pump',...
            'Mean FFT of PPI during Experiment')

        set(gca,'box','on','layer','top')
        set(gcf,'color','w')
        xlabel('Pump Wavelength')
        ylabel('FFT Intensity')    
        title(['t_2 = ' num2str(time(m)) ' fs'])    
        xlim(pulm)
    end

bfold = [basefolder sla proj sla subproj sla date sla tdname sla];

if ~exist(bfold,'dir')
    mkdir(bfold)
end

% save([bfold 'scan-by-scan_PPI_stats.mat'],'tm2','mn','sd','wl1','phint','mstep');

end
%% Reshape the data to show the full lab-time progression 
expstsvec = day(expsts)+hour(expsts)/24+minute(expsts)/6000;
expendvec = day(expend)+hour(expend)/24+minute(expend)/6000;

if length(time)>=10 && length(tvalues)~=1
    pnum = 10;
else
    pnum = 1;
end

yticks = [0:hours(explength)].*60;%.*length(time);
yticklabels = {0:hours(explength)};
% 
% oneytick =  [hour(ltime(1:nt2:end))+minute(ltime(1:nt2:end))/60];
% oneytick = oneytick(1:size(ftint(2)));

wllim = dsearchn(wl1',[800 500]');

[sfx sfy sfz] = size(ftint);
ftcont = reshape(ftint,sfx*sfy,sfz);

% tscan = ones(sfx,sfy).*linspace(1,nav,nav);
% tscan = tscan + time'./1e4;
% tscan = reshape(tscan,sfx*sfy,1);
tscan = ones(sfx*sfy,1);
tscan = linspace(0,hours(explength)*60,length(tscan));

figure(1); clf(1); hold all;
subplot(1,2,1)
contourf(wl1(wllim(1):wllim(end)),tscan,ftcont(:,wllim(1):wllim(end)))
xlabel('$\lambda_1$ (nm)')
ylabel('Experiment time (hours)')
title('Pump Spectrum as a function of Experiment Time')
grid on
set(gca,'layer','top','box','on','linewidth',2,'YTick',yticks,'YTickLabel',yticklabels)

subplot(1,2,2)
contourf(wl1(wllim(1):wllim(end)),(1:nav),squeeze(ftint(pnum,:,wllim(1):wllim(end))))
xlabel('$\lambda_1$ (nm)')
ylabel('Experiment time (a.u.)')
title(['Pump Spectrum as a function of Scan Number for $t_2$ = ' num2str(tvalues(pnum)) ' fs'])
grid on
set(gca,'layer','top','box','on','linewidth',2)

sgtitle(['Experiment length = ' char(explength)])
saveas(gcf,[foldname 'pumpdrift_' tdname '.fig'])
print([foldname 'pumpdrift_' tdname '.png'],'-dpng')
%% Choose which scans to include in re-averaging 

avsel = input('How to select avnum? User input - 0, Statistical calc - 1\n');

if ~avsel
    figure(2); clf(2);
    contour(squeeze(ftint(pnum,:,wllim(1):wllim(end))))    
    fprintf('Pick which average number to stop at\n')
    [a b] = ginput(1)
    b = round(b);
    nvec = b;
    ftfin = ftint(:,1:b,:);
    [sfx sfy sfz] = size(ftfin);
    ftfcont = reshape(ftfin,sfx*sfy,sfz);
    
    tscan = ones(sfx*sfy,1);
    tscan = linspace(0,hours(explength)*60*(b/nav),length(tscan));

    figure(2); clf(2);
    contourf(wl1(wllim(1):wllim(end)),tscan,ftfcont(:,wllim(1):wllim(end)))

else
    mntot = mean(ftcont,1);
    sdtot = std(ftcont,1);

    figure(2); clf(2); hold all;
    errorbar(wl1,mntot,sdtot)
    xlim([500 800])
end
%% Find outliers in chosen set of scans 

wllm = dsearchn(wl1',[pulm(end) pulm(1)]');
nv = 1:nvec;
outlr = zeros(length(tvalues),length(nv));
showplot = 0;

for m = 1:length(tvalues)
    
    clear tmp
    tmp = squeeze(ftfin(m,:,:));
    tmpsum = sum(tmp,2);
    tmout = isoutlier(tmpsum);
    outlr(m,:) = ~tmout;
    
    if showplot
        figure(1); clf(1);
        contourf(tmp);
        xlim([190 300])

        figure(2); clf(2); hold all;
        plot(nv,tmpsum)   
        plot(nv(tmout),tmpsum(tmout),'o')
    end
    
end
%% Re-average the data 

m = 1; n = 1;
filename = [tdname '_' num2str(time(m)) 'fs_' num2str(avnum(n)) '.bin'];
temp = read_bin([datfolder sla filename]);
clear bigdata

for m = 1:length(time)
    
    data = zeros(size(temp(4:end,3:end)));
    phint = zeros(size(temp(4:end,2)));
    
    for n = 1:length(nv)
        if length(tvalues) > 1
            oval = outlr(m,n);
        else
            oval = outlr(1,n);
        end
    
        if oval
        
            filename = [tdname '_' num2str(time(m)) 'fs_' num2str(avnum(n)) '.bin'];
            temp = read_bin([datfolder sla filename]);
            if exist('sfile','var')
                tmpd = temp(4:end,3:end) - squeeze(scat(n,:,:));
            else
                tmpd = temp(4:end,3:end);
            end
            data = tmpd + data;
            phint = temp(4:end,2) + phint;
            
        end
    end
    
    if length(tvalues) > 1
        sout = sum(outlr(m,:));
    else
        sout = sum(outlr(1,:));
    end
        
    avdat = data./sout;
    avphi = phint./sout;
    dat(m,:,:) = avdat;
    phi(m,:) = avphi;
    T2(m) = temp(1,1);

end

pumpspec = temp(2,3:end);
wlcal = temp(1,3:end);
mstepcal = temp(4:end,1);
%% Save the re-averaged data out
matfile = [foldname tdname '_av' num2str(length(nv)) sla];
matsavefile = [foldname tdname '_av' num2str(length(nv)) sla tdname '_av' num2str(length(nv)) '.mat'];

if ~exist(matfile,'dir')
    mkdir(matfile);
end

if exist(matsavefile,'file')
    load(matsavefile);
    clear NumScans Data PhInt
    
    NumScans = length(nv);
    tddat.Data = dat;
    tddat.PhInt = phi;
    tddat.mStep = mstepcal;
    tddat.wl3 = wlcal;
    tddat.t2 = T2;
    tddat.pump = pumpspec;
    
    if exist('sfile','var')
        tddat.scat = scat;
        tddat.sav = sdat;
    end
    
    save(matsavefile,'tddat','NumScans')
else
    NumScans = length(nv);
    tddat.Data = dat;
    tddat.PhInt = phi;
    tddat.mStep = mstepcal;
    tddat.wl3 = wlcal;
    tddat.t2 = T2;
    
    if exist('sfile','var')
        tddat.scat = scat;
        tddat.sav = sdat;
    end
    
    save(matsavefile,'tddat','NumScans')
end
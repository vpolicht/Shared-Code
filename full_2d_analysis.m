clear 
close all
clc;
WhichOS; clc;
%% Modify variables in this section

% Where is the data

proj = 'Perovskites';
subproj = '';
twist = '';
temp = 'RT';
date = '11-10-2021';

tdname = '2D_1-pervTest';

% Calibration set name
calibfile = '2DVIS-1_CALIB'; % phi1,phi2 = 0,0

% How remove background
ifftwin = 1;
polyfz = 0;

% Force calc reading in .bin data
forcecalc = 0; 
% Force calc performing filtering/circshift/fft/etc
bigforcecalc = 1;

% Pump limits (nm)
pulm = [506 750];
% Probe limits (nm)
prlm = [506 750];%[500 740];

[excitonev excitonlab] = excitonloc(proj,subproj,twist,temp);
zlim = [1.75 2.35];%1240./fliplr(pulm);%[];
nlim = [1.85 2.35];%1240./fliplr(pulm);%[];
%% Define Folder locations, etc.
phiname = {''};
calibfiles = {calibfile};

filename = [basefolder sla proj sla subproj sla date sla];
figfold = [filename sla tdname sla 'Figures' sla];

if ~exist(figfold,'dir')
    mkdir(figfold)
end
%% Load calibration data 
sfiles = length(calibfiles);
clear cal

calibname = calibfile;
calib = load([filename calibname '.dat']);

phint = calib(4:end,2);
datcal = calib(4:end,3:end);
wlcal = calib(1,3:end);
pumpspec = calib(2,3:end);
mstepcal = calib(4:end,1);

[srx sry] = size(datcal);
clear calib

figure(1); clf(1); hold all;
plot(wlcal,pumpspec);

figure(2); clf(2); hold all;
plot(mstepcal,phint)

figure(3); clf(3); hold all;
contour(wlcal,mstepcal,datcal,4,'color','k')

% Remove noisy ends of the t1 scan
t1lb = round(srx*.02);
t1ub = round(srx*.985);

mstepcal = mstepcal(t1lb:t1ub);
phint = phint(t1lb:t1ub);
datcal = datcal(t1lb:t1ub,:);

% Reduce w3 axis
w3lb = dsearchn(wlcal',pulm(1));
w3ub = dsearchn(wlcal',pulm(2));

pumpspec = pumpspec(w3lb:w3ub);
wlcal = wlcal(w3lb:w3ub);
datcal = datcal(:,w3lb:w3ub);

datcal(isnan(datcal)) = 0;
%% Subtract the pump spectrum
datcalf = datcal - mean(datcal,1);

figure(1);
plot(wlcal,pumpspec);

figure(2); 
plot(mstepcal,phint)

figure(3);
contour(wlcal,mstepcal,datcalf,4)

[srx sry] = size(datcalf);
%% Perform w1 calibration

nfft = 2^14;

w3cal = (c*1e15)./wlcal; % in Hz

% Pseudofrequency axis
dx = mean(abs(diff(mstepcal)));
f_fft = ((1/dx)/2).*linspace(0,1,nfft/2+1);

ftdat = fft(datcalf,nfft,1);
ftdatcal = abs(ftdat(1:nfft/2+1,:));
[a maxcorr] = max(ftdatcal);

tmp = abs(diff(diff(diff(maxcorr))));
mm = find(tmp<5);

b = maxcorr(mm);

% Fitting axis
ll = abs(mstepcal(1) - mstepcal(end));
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
%% Filter along t1
fact = .7;

SRX = length(mstepcal);
nmr = @(a) a./max(a);

figure(2); clf(2); hold all;
plot(1:SRX,nmr(sum(datcalf,2)))

filtype = input('Choose t1 filter: 0 - Gaussian, 1 - Tukey\n');
tmpax = 1:length(phint);

datcalf = squeeze(datcalf);
[filt] = deal(zeros(2*SRX,1));

if ~filtype
    [a,x0] = max(abs(sum(datcalf(1:SRX/2,:),2)));
    filt = (gaussmf(tmpax,[fact*SRX/2 x0]))';
    filt = filt(1:SRX);
    filt = filt - min(filt);
else
    filt = tukeywin(2*SRX,fact);
    filt = filtcirc(SRX+1:end);
end

filt = nmr(filt);

plot(1:SRX,filt)
plot(1:SRX,nmr(sum(datcalf.*filt,2)))

datcalf = datcalf.*filt;
phintf = phint.*filt;
%% Pad data & load in bin0 parameters

srx = length(phintf);
nfft = 2^(nextpow2(srx)+1);
padsize = nfft - srx;

dx = mean(abs(diff(mstepcal)));
f_fft = ((1/dx)/2).*linspace(0,1,nfft/2+1);
  
w1 = p(1)*f_fft + p(2);
wl1 = (c*1e15)./w1;
w1lm = dsearchn(wl1',pulm')';

df = mean(abs(diff(w1)));
t1tmp = ((1e15/df)/2).*linspace(0,1,nfft/2+1);

[a mx] = max(abs(phintf));

phintpad = padarray(phintf,padsize,'post');
datpad = padarray(squeeze(datcalf),padsize,'post');
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

for MM = 1:length(yax)
    M = yax(MM);
    ftv(MM,:) = fft(circshift(phintpad,-M),nfft);
    
    ff = fit(wax',smoothdata((angle(ftv(MM,wslm(1):wslm(2)))),10)','poly1');
%     fsnd(MM) = ff.p1;
    fslope(MM) = ff.p1;
    foff(MM) = ff.p2;
    fitz(MM,:) = feval(ff,wax);
    
    immin(MM) = min(sum(imag(ftv(MM,wslm(1):wslm(2))),2));
end

figure(1); clf(1); hold all;
yyaxis left
plot(yax,phintf(yax))
yyaxis right
hold all;
plot(yax,abs(fslope./max(fslope)),'o')
plot(yax,foff,'rd')

[a b] = min(abs(fslope));
xline(yax(b));
ylim([-1 1])

bin0 = yax(b);
%% Select bin0 based on phase

good = 1;
unwrapd = 1;
mnpump = pumpspec;

while good == 1
clear phintsf datpsf fftphint fftdat phipha phiamp datpha datamp

    phintsf = circshift(phintpad,-bin0);
    datpsf = circshift(squeeze(-1.*datpad),-bin0);
    datpss = sum(datpsf,2);

    fftphint = fft(phintsf,nfft);%./max(phintsf),nfft);
    fftdat = fft(datpss,nfft);%./max(datpss(m,:,:)),nfft);

    phipha = angle(fftphint);
    phiamp = abs(fftphint);

    datpha = angle(fftdat);
    datamp = abs(fftdat);
     
    figure(3); clf(3); 
    set(gca,'Position',[680 189 560 789])
    subplot(2,1,1);
    hold all;
    yline(0,'handlevisibility','off');
    plot(phint./max(phint)) 
    xline(bin0,'handlevisibility','off');
    xlim([mx - 2e2 mx + 2e2])
    ylim([-1.5 1.5])
    title(['$bin_0$ = ' num2str(bin0) ', filt w. ' num2str(fact)])
    set(gca,'box','on','XTickLabel','')
    xlabel('Motor Step (a.u.)')

    subplot(2,1,2);
    yyaxis left
    hold all;
    plot(wl1,phiamp(1:nfft/2+1)./max(phiamp(1:nfft/2+1)))
    plot(wl1,datamp(1:nfft/2+1)./max(datamp(1:nfft/2+1)))
    plot(wlcal,pumpspec./max(pumpspec))
    ylim([0 1.1])
    
    yyaxis right
    hold all;
    plot(wl1,phipha(1:nfft/2+1))
    plot(wl1,datpha(1:nfft/2+1))

    yyaxis left 
    ylabel('FFT amplitude')
    
    yyaxis right       
    yline(0,'handlevisibility','off');
    ylabel('FFT phase')
    xlim([450 800]); xlabel('$\lambda$ (nm)')
    ylim([-pi*1.5 pi*1.5]); 
    set(gcf,'color','w');
    set(gca,'box','on','YTick',[-pi 0 pi],'YTickLabel',{'-\pi','0','\pi'},...
        'layer','top','linewidth',1.5)
    yline(-pi,':','handlevisibility','off');
    yline(0,':','handlevisibility','off');
    yline(pi,':','handlevisibility','off');
    
    if unwrapd
        phipha = unwrap(phipha);
    end
    
    phipha = flipud(phipha);
    expha = exp(1i.*phipha);
    
    good = input('Change bin0? 0 - no, 1 - yes\n');
    if good 
        bin0 = input('New bin0 = \n');
        clear phintmp
    end
end

saveas(gcf,[figfold 'bin0_' calibname '.fig'])
saveas(gcf,[figfold 'bin0_' calibname '.png'])

phintsft = phintsf;
phiamp = phiamp(1:nfft/2+1);
datamp = datamp(1:nfft/2+1);
expha = expha(1:nfft/2+1);
clear phintsf datpsf fftphint fftdat phipha datpha datamp

% Find maximum of circleshifted phase
x0 = find(abs(phintsft) >= max(abs(phintsft)));

figure(4); clf(1); hold all;
plot(wlcal,pumpspec./max(pumpspec),'k:','linewidth',1.5)
plot(wl1,phiamp./max(phiamp),'r:','linewidth',1.5)
%% Filter out pump spectrum from datcal
w1 = c*1e15./wl1;
dw = mean(abs(diff(w1)));

if ifftwin

    tfftp = fft(datcal,nfft);

    B = pulm(1) + (pulm(end)-pulm(1))/2;
    B = (c*1e15)/B;
    t1filt = sgauss(w1,75*dw,B,6);
    t1filt = cat(2,t1filt,fliplr(t1filt(2:end-1)))';
    
    datcalosc = ifft(tfftp.*t1filt,nfft,1);
    datcalosc = datcalosc(1:size(datcal,1),:);
elseif polyfz    
    Pp = polyfit(tmpax,datcal,3);
    polbck = polyval(Pp,tmpax);
    datcalosc = datcal - polbck';
else
    datcalosc = datcal - mean(datcal,1);
end

clear data B
%% Load in p/m45 data
datfile = tdname;
fprintf('Load in 2D data...\n')
tic
if ~exist([filename datfile sla 'processed-' datfile '.mat'],'file') || bigforcecalc
%% Load in from bin0 or .mat 
clear tddat

if ~exist([filename datfile sla datfile '.mat'],'file') || forcecalc
    files = dir([filename sla]);

    % Look for all files matching the 2D data name
    count = 1;
    for m = 1:length(files)
        if contains(files(m).name,datfile) && contains(files(m).name,'.bin')
            bfile(count) = files(m);
            count = count+1;
        end

    end

    if exist('bfile','var')

        for m = 1:length(bfile) 
            temp = read_bin([filename sla bfile(m).name]);

            tphi(m,:) = temp(4:end,2);
            tData(m,:,:) = temp(4:end,3:end);
            twpr = temp(1,3:end);
            tmStep = temp(4:end,1);
            tpumpspec(m,:) = temp(2,3:end);
            
            tT2(m) = temp(1,1);

        end

        [va ord] = sort(tT2);
        
        tddat.t2 = tT2(ord);
        tddat.Data = tData(ord,:,:);
        tddat.PhInt = tphi(ord,:);
        tddat.mStep = tmStep;
        tddat.wl3 = twpr;
        tddat.pump = tpumpspec(ord,:);
        
        save([filename datfile sla datfile '.mat'],'tddat','bin0','expha','wl1');

    else
        tmp = load([filename datfile '.mat']);

        for m = 1:length(tmp.T2)
           tddat.Data = tmp.Data{m};
           tddat.PhInt = tmp.PhInt{m};
        end

        tddat.mStep = tmp.mStep;
        tddat.wl3 = tmp.wpr;
        tddat.t2 = tmp.T2;
        srx = length(tmp.mStep);
        sry = length(tmp.wpr);
    end

    plot(tddat.wl3,tddat.pump(1,:)./max(tddat.pump(1,:)))
    clear tmStep tData tt2 twpr temp
else
    load([filename datfile sla datfile '.mat']);

    if size(tddat.t2,2) == 1
        tddat.t2 = tddat.t2';
        tddat.mStep = tddat.mStep';
    end
end
toc
%% Prepare phasing for p/m45 data
fprintf('Processing the 2D data before FFT...\n')
tic
tddat.SRZ = length(tddat.t2);

tind = dsearchn(tddat.t2',50);
if tind == 1 && tddat.SRZ > 1
    tind = tind+1;
end

tmpdat = sum(squeeze(tddat.Data(tind,:,:)),1)';
% tmpdat = abs(tmpdat)-min(abs(tmpdat));

% Force data w3 limits to match the calibration wl axis
wl3lm = [dsearchn(tddat.wl3',wlcal(1)) dsearchn(tddat.wl3',wlcal(end))];
tddat.pwl3 = tddat.wl3(wl3lm(1):wl3lm(end));
tddat.pw3 = 1240./tddat.pwl3; % eV
dpw3 = mean(abs(diff(tddat.pw3)));
tddat.SRY = length(tddat.pwl3);

tddat.wl1 = wl1;
tddat.w1 = (c*1e15)./tddat.wl1;

tddat.mStep = tddat.mStep(t1lb:t1ub);
tddat.PhInt = tddat.PhInt(:,t1lb:t1ub);

try 
    sav = tddat.sav;
catch
end
% 
% if exist('scat','var')
% %     for m = 1:length(tddat.t2)
% %         tddat.Data(m,:,:) = squeeze(tddat.Data(m,:,:)) - sav;    
% %     end
% tddat.sav = tddat.sav(t1lb:t1ub,wl3lm(1):wl3lm(end));
% end

wl1lm = dsearchn(wl1',pulm')';
tddat.pwl1 = wl1(wl1lm(end):wl1lm(1));
tddat.pw1 = 1240./tddat.pwl1; % eV
tddat.dw = mean(abs(diff(w1)));
dpw1 = mean(abs(diff(tddat.pw1)));
%% Remove background from p/m45 data
cshft = bin0/2;

if length(tddat.t2) == 1
    dim = 1;
else
    dim = 2;
end

tddat.datasub = zeros(size(tddat.Data));
SRX = size(tddat.Data,2);

if ifftwin

%     tfftp = fft(circshift(tddat.Data,cshft,2),nfft,dim);
    tfftp = fft(tddat.Data,nfft,dim);
    
    B = pulm(1) + (pulm(end)-pulm(1))/2;
    B = (c*1e15)/B;
    t1filt = sgauss(tddat.w1,75*tddat.dw,B,6);
    t1filt = cat(2,t1filt,fliplr(t1filt(2:end-1)))';

    if exist('sav','var')
%         tfts = fft(circshift(tddat.sav,cshft),nfft,1);
        tfts = fft(tddat.sav,nfft,1);
        tmps = ifft(tfts.*t1filt,nfft,1);
%         tddat.savsub = circshift(tmps(1:SRX,:),-cshft);
        tddat.savsub = tmps(1:SRX,:);

        figure(2); clf(2); hold all;
        plot(tddat.wl1,mean(abs(tfts(1:nfft/2+1,:)),2))
        plot(tddat.wl1,t1filt(1:nfft/2+1))
        plot(tddat.wl1,mean(abs(t1filt(1:nfft/2+1).*tfts(1:nfft/2+1,:)),2))
        xlim([pulm(1)-300 pulm(end)+300])
        
        figure(3); clf(3); hold all;
        plot(mean(tddat.sav,2))
        plot(mean(tddat.savsub,2))
        yline(0)
                    
        figure(1); clf(1); hold all;
        plot(mean(squeeze(tddat.Data(dsearchn(tddat.t2',500),:,100:200)),2))
        plot(mean(squeeze(tddat.datasub(dsearchn(tddat.t2',500),:,100:200)),2))
        yline(0)

        figure(4);clf(4); hold all;
        plot(mean(tddat.savsub,2))
        plot(mean(squeeze(tddat.datasub(dsearchn(tddat.t2',500),:,100:200)),2))
        plot(mean(squeeze(tddat.datasub(dsearchn(tddat.t2',500),:,:))-tddat.savsub,2))

    end
            
    for n = 1:tddat.SRZ
        if dim ~= 1
            tmpp = ifft(squeeze(tfftp(n,:,:)).*t1filt,nfft,1);   
        else
            tmpp = ifft(tfftp.*t1filt,nfft,1);
        end
%         tddat.datasub(n,:,:) = circshift(tmpp(1:SRX,:),-cshft,1);
        tddat.datasub(n,:,:) = tmpp(1:SRX,:);
    end
    
    adl = '-ifft';
elseif polyfz    
    for n = 1:tddat.SRZ
        for m = 1:tddat.SRY
            Pp = polyfit(tmpax,squeeze(tddat.Data(n,:,m)),3);
            polbck = polyval(Pp,tmpax);
            tddat.datasub(n,:,m) = tddat.Data(n,:,m) - polbck';
        end
    end
    adl = '-poly';
else
    for n = 1:tddat.SRZ
        tddat.datasub(n,:,:) = squeeze(tddat.Data(n,:,:)) - mean(squeeze(tddat.Data(n,:,:)),1);
    end
    adl = '-mean';
end

tddat.SRX = length(tddat.mStep);

tddat.datasub = tddat.datasub(:,t1lb:t1ub,wl3lm(1):wl3lm(end));
if exist('sav','var')
    tddat.savsub = tddat.savsub(t1lb:t1ub,wl3lm(1):wl3lm(end))
    
    figure(3);
    plot(t1lb:t1ub,mean(tddat.savsub,2))

end

figure(1);
plot(t1lb:t1ub,mean(squeeze(tddat.datasub(dsearchn(tddat.t2',500),:,100:200)),2))

clear data B
%% Filter along t1, Pad data & circshift

tddat.dataf = zeros(size(tddat.datasub));
[tddat.datapad,tddat.datasf] = deal(zeros(length(tddat.t2),size(tddat.datasub,2)+padsize,size(tddat.datasub,3)));

tddat.phf = zeros(size(tddat.PhInt));
[tddat.phpad,tddat.phsf] = deal(zeros(length(tddat.t2),size(tddat.datasub,2)+padsize));

for n = 1:tddat.SRZ
   tddat.dataf(n,:,:) = squeeze(tddat.datasub(n,:,:)).*filt;
   tddat.datapad(n,:,:) = padarray(squeeze(tddat.dataf(n,:,:)),padsize,'post');
   tddat.datasf(n,:,:) = circshift(squeeze(tddat.datapad(n,:,:)),-bin0);

   tddat.phf(n,:) = tddat.PhInt(n,:).*filt';
   tddat.phpad(n,:) = padarray(tddat.phf(n,:)',padsize,'post');
   tddat.phsf(n,:) = circshift(tddat.phpad(n,:),-bin0);
end

if exist('sav','var')
    tddat.savf = squeeze(tddat.savsub.*filt);
    tddat.savpad = padarray(tddat.savf,padsize,'post');
    tddat.savsf = circshift(tddat.savpad,-bin0);
end
toc
%% FFT Data
fprintf('FFT the Data...\n')
tic;
smdat = 0;
count = 1;

tddat.expha = expha(wl1lm(end):wl1lm(1));
tddat.ftdat = zeros(tddat.SRZ,tddat.SRY,length(tddat.pw1));

for MM = 1:tddat.SRZ
   clear data tmpdat legeset tmp tmpp tmpm tmps

    % Pick a filter type for along t1
    nmr = @(a) a./max(a);

    % Perform FFT
    tddat.SRY = size(squeeze(tddat.datapad(MM,:,:)),2);

    tmpdat = fft(squeeze(tddat.datasf(MM,:,:)),nfft);
    tmpph = fft(tddat.phsf(MM,:),nfft);
   
    tddat.ftph(MM,:) = tmpph(wl1lm(end):wl1lm(1));
    tddat.ftdat(MM,:,:) = -1.*permute(tmpdat(wl1lm(end):wl1lm(1),:).*tddat.expha,[2 1]);

    count = count+1;

end   

if exist('sav','var')
    tmpftsav = fft(tddat.savsf,nfft,1);
    tddat.ftsav = -1.*permute(tmpftsav(wl1lm(end):wl1lm(1),:).*tddat.expha,[2 1]);
end

figure(4); 
set(gca,'ColorOrderIndex',1);
plot(1240./tddat.pw1,abs(real(tddat.ftph(tind,:)))./max(abs(real(tddat.ftph(tind,:)))),'linewidth',1.5)
plot(1240./tddat.pw1,abs(sum(real(squeeze(tddat.ftdat(tind,:,:))),1))./max(abs(sum(real(squeeze(tddat.ftdat(tind,:,:))),1))),'linewidth',1.5)
xlim([480 1000])
legend('Calib Pump','Calib PPI','Data PPI','Data Int over \omega_1');
set(gca,'box','on','layer','top')
xlabel('\lambda_1 (nm)')
ylabel('Norm Amp (a.u.)')
title('Comparison of Calibration & Data')
set(gcf,'color','w')

for m = 1:length(excitonev)
    xline(1240./excitonev(m),':','handlevisibility','off','linewidth',1.5);
end

print([figfold sla 'PPI_v_PumpSpec.png'],'-dpng')
 toc
%% Save out processed data
fprintf('Saving out data...\n')
tic

tmp = tddat;

clear tddat;

tddat.t2 = tmp.t2;
tddat.pwl3 = tmp.pwl3;
tddat.pw3 = tmp.pw3;
tddat.pwl1 = tmp.pwl1;
tddat.pw1 = tmp.pw1;
tddat.Data = tmp.Data;
tddat.datasub = tmp.datasub;
tddat.dataf = tmp.dataf;
tddat.datapad = tmp.datapad;
tddat.datasf = tmp.datasf;
tddat.adl = adl;
tddat.mstep = tmp.mStep;

save([filename datfile sla 'filtered-' datfile '.mat'],'tddat','bin0','-v7.3');

clear tddat;

tddat.t2 = tmp.t2;
tddat.pwl3 = tmp.pwl3;
tddat.pw3 = tmp.pw3;
tddat.pwl1 = tmp.pwl1;
tddat.pw1 = tmp.pw1;
tddat.ftdat = tmp.ftdat;
tddat.adl = adl;
tddat.ftph = tmp.ftph;
save([filename datfile sla 'processed-' datfile '.mat'],'tddat','bin0');

toc
clear tmp
else
    load([filename datfile sla 'processed-' datfile '.mat'])
end
%% Shift t0

[tddat.SRZ tddat.SRY tddat.SRX] = size(tddat.ftdat);
exnum = 3;

try
if length(tddat.t2) >= 3

    e1 = dsearchn(tddat.pw1',excitonev(exnum));
    e3 = dsearchn(tddat.pw3',excitonev(exnum));

    t2tr = squeeze(real(tddat.ftdat(:,e3,e1)));
    ttlim = dsearchn(tddat.t2',100);
    [A tstp] = max(abs(t2tr(1:ttlim)));
    A = A.*sign(t2tr(tstp));
    tstr = tddat.t2(tstp)-80;

    figure(1); clf(1); hold all;
    plot(tddat.t2,t2tr)
    t2tr = t2tr - mean(t2tr(1:dsearchn(tddat.t2',tddat.t2(tstp)-70)));
    plot(tddat.t2,t2tr)
    plot(tddat.t2(2:end),abs((diff(squeeze(tddat.ftdat(:,e3,e1))))),'o')
    set(gca,'box','on','layer','top','linewidth',1.5,'fontname','arial')
    xline(0,'handlevisibility','off','linestyle',':');
    yline(0,'handlevisibility','off');
    xlabel('t_2 (fs)')
    set(gcf,'color','w')

    xline(tddat.t2(tstp),'r');
    xline(tstr,'g');
    outlie = excludedata(tddat.t2,t2tr','domain',[tstr tddat.t2(tstp)]);

    ft = fit(tddat.t2(~outlie)',t2tr(~outlie),'gauss1','robust','LAR',...
        'Upper',[A Inf 40],'Lower',[A -Inf 0]);
    fte = feval(ft,tddat.t2);
    fcoeff = coeffvalues(ft);

    fwhm = 2*sqrt(log(2))*fcoeff(3);
    t2zero = dsearchn(tddat.t2',fcoeff(2)-fwhm/2);
    xline(tddat.t2(t2zero),'b');
    xline(0,'linewidth',3);

    plot(tddat.t2,fte,'linewidth',1.5,'linestyle','--')

    figure(1)
    ttemp = tddat.t2-tddat.t2(t2zero);

    plot(ttemp,t2tr)

    xlim([-250 250])

    H = input('Does the new t2 = 0 look good? 1 - yes, 0 - no\n');
    while H == 0

        fprintf('Select best t2 = 0:\n')
        [l n] = ginput(1);
        if size(ttemp,1) == 1
            ttemp = ttemp';
        end
       
        t2zt = dsearchn(ttemp,l);
        t2zero = dsearchn(tddat.t2',ttemp(t2zt));
        ttemp = tddat.t2 - tddat.t2(t2zero);

        plot(tddat.t2-tddat.t2(t2zero),(squeeze(tddat.ftdat(:,e3,e1))))
        H = input('Does the new t2 = 0 look good? 1 - yes, 0 - no\n');
    end
    clear b e1 e3 ttemp

    tddat.t2 = tddat.t2 - tddat.t2(t2zero);   
end
catch
end
%% Calculate the Pump-Subtracted data

for m = 1:tddat.SRZ
    tddat.ftPS(m,:,:) = squeeze(real(tddat.ftdat(m,:,:)))./(abs(tddat.ftph(m,:))./max(abs(tddat.ftph(m,:))));
end
%% Plot spectrum

smdat = 5;
zmplt = 0;
count = 1;
PS = 1;

nlmx = dsearchn(tddat.pw1',nlim');
nlmy = flipud(dsearchn(tddat.pw3',nlim'));

tddat.dw1 = mean(abs(diff(tddat.pw1)));
tddat.dw3 = mean(abs(diff(tddat.pw3)));

if smdat
    Adl = [tddat.adl '-sm'];
else 
    Adl = tddat.adl;
end

if ~PS
    [cmax b] = max(max(max(real(tddat.ftdat(:,nlmy(1):nlmy(end),nlmx(1):nlmx(end))))));
    [cmin b] = min(min(min(real(tddat.ftdat(:,nlmy(1):nlmy(end),nlmx(1):nlmx(end))))));
    AAdl = [Adl];
else
    [cmax b] = max(max(max(tddat.ftPS(:,nlmy(1):nlmy(end),nlmx(1):nlmx(end)))));
    [cmin b] = min(min(min(tddat.ftPS(:,nlmy(1):nlmy(end),nlmx(1):nlmx(end)))));
    AAdl = [Adl '-PS'];
end

clear color
colors = b2r(cmin,cmax);

for MM = 1:tddat.SRZ 
    clear tdmap
    
    if ~PS
        tdmap = squeeze(real(tddat.ftdat(MM,:,:)));
    else
        tdmap = squeeze(tddat.ftPS(MM,:,:));
    end
    
    if smdat 
        smv = smdat*tddat.dw3;
        smx = round(smv/tddat.dw1);
        
        tdmap = smoothdata(tdmap,1,'movmean',smdat);
        tdmap = smoothdata(tdmap,2,'movmean',smx);
    end
            
    lmin = min(min(tdmap(nlmy(1):nlmy(end),nlmx(1):nlmx(end))));
    lmax = max(max(tdmap(nlmy(1):nlmy(end),nlmx(1):nlmx(end))));
    if abs(lmax) > abs(lmin)
        intv = lmax/10;
        contlp = intv:intv:lmax;
        contln = -1.*contlp;
    else
        intv = lmin/10;
        contln = intv:intv:lmin;
        contlp = -1.*contln;
    end
    
    figure(1); clf(1); hold all;
    contourf(tddat.pw1,tddat.pw3,tdmap,linspace(cmin,cmax,30),'linecolor','none');
    colormap(colors); crange = caxis; 

    contour(tddat.pw1,tddat.pw3,tdmap,contln,'linecolor','k','linestyle','--');
    contour(tddat.pw1,tddat.pw3,tdmap,contlp,'linecolor','k');
    caxis([cmin,cmax])
    
    line(tddat.pw1,tddat.pw1,'color','k')
    daspect([1 1 1])
    set(gca,'box','on','Xlim',[tddat.pw1(1) tddat.pw1(end)],'YLim',[tddat.pw3(end) tddat.pw3(1)],...
        'layer','top','linewidth',1.5)
    xlabel('$\hbar\omega_1 (eV)$','interpreter','latex')
    ylabel('$\hbar\omega_3 (eV)$','interpreter','latex')
    colorbar
    title(['$t_2$ = ' num2str(tddat.t2(MM)) ' fs'])
    
     for X = 1:length(excitonev)
        xline(excitonev(X),'k:','linewidth',1);
        yline(excitonev(X),'k:','linewidth',1);
     end
     
     print([figfold sla 'real_abs-t2-' num2str(tddat.t2(MM)) AAdl '.png'],'-dpng')

    if zmplt
        xlim([zlim(1) zlim(2)])
        ylim([zlim(1) zlim(2)])
        print([figfold sla 'zoom-real_abs-t2-' num2str(tddat.t2(MM)) AAdl '.png'],'-dpng')
    end
    count = count + 1;
end
%% Plot t2 traces at the exciton positions

smdat = 5;
smst = 5; 

avx = floor(0.005/tddat.dw1);
avy = floor(0.005/tddat.dw3);
if ~avx
    avx = 1;
end

if sum(diff(unique(diff(tddat.t2))) > 100)
    logplot = 1;
else
    logplot = 0;
end

for m = 1:length(excitonev)
    figure(m); clf(m); hold all;
    ax1 = gca;
    xind = dsearchn(tddat.pw1',excitonev(m));
    
    if xind == 1
        xind = avx+1;
    elseif xind == tddat.SRX
        xind = tddat.SRX - avx-1;
    end
        
    for n = 1:length(excitonev)
        
        yind = dsearchn(tddat.pw3',excitonev(n));
        
        if yind == 1
            yind = avy+1;
        elseif yind == tddat.SRY
            yind = tddat.SRY - avy - 1;
        end
        
        if ~PS
            t2tr = mean(mean(squeeze(real(tddat.ftdat(:,yind-avy:yind+avy,xind-avx:xind+avx))),2),3);
        else
            t2tr = mean(mean(squeeze(real(tddat.ftPS(:,yind-avy:yind+avy,xind-avx:xind+avx))),2),3);
        end
        
        if logplot
    
            if n == 1
                ax1 = subplot(121);
                hold all;
                ax2 = subplot(122);
                hold all;
            end
            
            axsp = 'log';
            spc = min(diff(tddat.t2(tddat.t2 > 0 & tddat.t2 < 1000)));
            tmp = find(diff(tddat.t2) == spc);
            lnl = tmp(end) + 1;
            
            axes(ax1)
            aa = plot(tddat.t2,t2tr,':','handlevisibility','off')
            set(ax1,'ColorOrderIndex',n)
            t2sm = [t2tr(1:dsearchn(tddat.t2',0)+smst)',smoothdata(t2tr(dsearchn(tddat.t2',0)+smst+1:end),'movmean',smdat)'];
            ab = plot(tddat.t2,t2sm)
      
            axes(ax2)
            bb = plot(tddat.t2,t2tr,':','handlevisibility','off')
            set(ax2,'ColorOrderIndex',n)
            t2sm = [t2tr(1:dsearchn(tddat.t2',0)+smst)',smoothdata(t2tr(dsearchn(tddat.t2',0)+smst+1:end),'movmean',smdat)'];
            bc = plot(tddat.t2,t2sm)
            
            set(ax1,'units','normalized','position',[0.1 0.1 0.2 0.8]);
            set(ax1,'xlim',[tddat.t2(1) tddat.t2(lnl)])
            set(ax2,'units','normalized','position',[0.3 0.1 0.6 0.8]);
            set(ax2,'xscale','log','xlim',[tddat.t2(lnl) tddat.t2(end)],'yticklabel','');
        else           
            plot(tddat.t2,t2tr,':','handlevisibility','off')
            set(ax1,'ColorOrderIndex',n)
            t2sm = [t2tr(1:dsearchn(tddat.t2',0)+smst)',smoothdata(t2tr(dsearchn(tddat.t2',0)+smst+1:end),'movmean',smdat)'];
            plot(tddat.t2,t2sm)
        end
        
        ymn(n) = floor(1.1*min(t2tr)*100)/100;
        ymx(n) = round(1.1*max(t2tr)*100)/100;

        leg{n} = ['$\hbar\omega_3$ = ' num2str(round(1e2.*tddat.pw3(yind))/1e2) ' eV'];
    end
        
    xline(0,'handlevisibility','off');
    yline(0,'handlevisibility','off');
    title(['$\hbar\omega_1$ = ' num2str(round(1e2.*tddat.pw1(xind))/1e2) ' eV'],'interpreter','latex')
    legend(leg,'Interpreter','latex')
    set(gcf,'color','w')
    xlabel('$t_2$ (fs)','interpreter','latex')
    
    if logplot
        tcks = linspace(min(ymn),max(ymx),11);

        axes(ax1)
        ylabel('2DES Amplitude (a.u.)')
        set([ax1 ax2],'box','on','layer','top','linewidth',1)
        set([ax1 ax2],'ylim',[tcks(1) tcks(end)],'ytick',tcks,'box','off');
        set([ax1 ax2],'box','on','layer','top','linewidth',1.5,...
            'fontname','arial','xminortick','on');%'xgrid','on','ygrid','on','gridcolor','k'
        xline(0); yline(0); axes(ax2); yline(0);
        axes(ax1)
    else
        set(gca,'box','on','layer','top','linewidth',1)
        ylabel('2DES Amplitude (a.u.)')
    end
    
    print([figfold sla 'hw1-' excitonlab{m} '-t2trace' AAdl '.png'],'-dpng')
end        
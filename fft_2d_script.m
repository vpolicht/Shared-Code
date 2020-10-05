% Load in the data S(t1,t2,w3) and perform the FFT to obtain S(w1,t2,w3)
% Saves out the data in the form of S(t2,w1,w3)
% Veronica Policht, 10/2020

% This scipt requires a small amount of user input:
    % When prompted, you must select a "bin0" to define t1 = 0 fs
    % When prompted, you must select whether to window t1 with a Guassian
        % or tukey window function

 clear;
 close all;
 clc;
 
 c = 299.792458; % Speed of light in nm/fs
 %% Modify variables in this section
 
 % Where is the data
 filename = ['/Volumes/Utini/Data/V2O3/7-24-2020/'];
 % Data set name
 dataname = '2D_WMS2ua-spot14-1000fs-RT-1_av46';
 
% Calibration set name
calibname = '2DVIS-4_CALIB';

% Pump limits (nm)
pulm = [490 600];
% Probe limits (nm)
prlm = [490 590];

% Toggle to plot 2D maps: 0 - no, 1 - yes
plotcont = 1;

% If wish to smooth the data: 0 - no, 1 - yes
smdat = 0;

% Toggle to plot pre-FFT maps: 0 - no, 1 - yes
noplot = 1;

% Starting bin0
bin0 = 100;

%% Load and modify Calibration

% Load data
calib = load([filename calibname '.dat']);

phint = calib(4:end,2);
datcal = calib(4:end,3:end);
wlcal = calib(1,3:end);
pumpspec = calib(2,3:end);
mstepcal = calib(4:end,1);

[srx sry] = size(datcal);

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

if plotcont
    figfold = [filename sla dataname sla 'Figures' sla];
    
    if ~exist(figfold,'dir')
        mkdir(figfold)
    end
end


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

a = isoutlier(maxcorr);
st = 1:length(w3cal);
st(a) = [];
b = round(interp1(st,maxcorr,1:length(w3cal)));

% Fitting axis
ll = abs(mstlm(1) - mstlm(end));
pp = 1/(ll/srx*nfft);
psfr = maxcorr*pp.*1e15;

figure(1); clf(1); hold all;
contour(w3cal,f_fft,ftdatcal)
plot(w3cal,f_fft(b),'color','k')

p = polyfit(f_fft(b),w3cal,1);

w1 = p(1)*f_fft + p(2);
wl1 = (c*1e15)./w1;

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
%% Load in 2D data [S(t1,t2,w3)]

tddat = load([filename dataname '.mat']);

t2 = tddat.T2;
SRZ = length(t2);

wl3 = tddat.wpr;
wl3lm = dsearchn(wl3',prlm')';
pwl3 = wl3(wl3lm(1):wl3lm(end));
pw3 = 1240./pwl3; % eV
SRY = length(pwl3);

wl1lm = dsearchn(wl1',pulm')';
pwl1 = wl1(wl1lm(end):wl1lm(1));
pw1 = 1240./pwl1; % eV
% srx = length(pwl1);

mstep = tddat.mStep(bin0:end);
SRX = length(mstep);

phapad = -1.*expha.*ones(length(expha),SRY);
phapad = phapad(wl1lm(end):wl1lm(1),:);
phapad = permute(phapad,[2 1]);

count = 1;

t1ax = ((1e15/df)/2).*linspace(0,1,SRX);
%% Calculate the FFT wrt t1 for each t2 point

sw1t2w3 = zeros(SRZ,SRY,length(pwl1));

for MM = 1:SRZ
    clear data
    
    data = tddat.Data{MM}(bin0:end,wl3lm(1):wl3lm(end));

    if noplot
        figure(2); clf(2); subplot(1,2,1); hold all
        contour(pw3,t1ax,(data))
        xlabel('Detection Frequency (eV)')
        ylabel('t_1 (fs)')
        set(gca,'box','on','layer','top')
        colorbar('location','northoutside')
    end
    
    %% Subtract out power fluctuations during t1
    tmpax = 0:SRX-1;

    datasub = zeros(size(data));

    for m = 1:SRY
        P = polyfit(tmpax,data(:,m),3);
        polbck = polyval(P,tmpax);
        datasub(:,m) = data(:,m) - polbck';
    end
    
    if noplot
        figure(2);
        subplot(1,2,2); hold all
        contour(pw3,t1ax,abs(datasub),[.2 .4 .6 .8 1].*max(max(abs(datasub))))
        xlabel('Detection Frequency (eV)')
        ylabel('t_1 (fs)')
        set(gca,'box','on','layer','top')
        colorbar('location','northoutside')
        sgtitle(['S(t_1,\omega_3) at t_2 = ' num2str(t2(MM)) ' fs'])
    end
    
    %% Pick a filter type for along t1
    nmr = @(a) a./max(a);

    if count == 1

        figure(3); clf(3); hold all;
        plot(1:SRX,nmr(sum(datasub,2)))

        filtype = input('Choose t1 filter: 0 - Gaussian, 1 - Tukey\n');

        filt = zeros(2*SRX,1);
        tmpax = 1:2*SRX;
        fact = .8;

        if ~filtype
            filt = gaussmf(tmpax,[fact*SRX/2 SRX+x0-1]);
        else
            filt = tukeywin(2*SRX,fact);
        end

        % filt = circshift(filt,-bin0);
        filt = filt(SRX+1:end);
        filt = nmr(filt);

        plot(1:SRX,filt)

    end

%% Apply filter along t1

    bckg = repmat(mean(datasub,1),size(datasub,1),1);
    datasb = datasub - bckg;

    datafilt = zeros(size(datasb));

    for m = 1:SRY
        datafilt(:,m) = datasb(:,m).*filt;
    end
    
    if noplot
        figure(2); subplot(1,2,2);
        contour(pw3,t1ax,abs(datafilt),[.2 .4 .6 .8 1].*max(max(abs(datafilt))),'color','k')
    end
    
%% Perform FFT

    datapad = padarray(datafilt,padsize,'post');
    SRY = size(datapad,2);

    fftdat = fft(datapad,nfft);
    
    sw1t2w3(MM,:,:) = permute(fftdat(wl1lm(end):wl1lm(1),:),[2 1]);
    
    if smdat 
        sw1t2w3(MM,:,:) = smoothdata(sw1t2w3(MM,:,:),3,'movmean',5);
        sw1t2w3(MM,:,:) = smoothdata(sw1t2w3(MM,:,:),2,'movmean',5);
    end
    
    sw1t2w3(MM,:,:) = squeeze(sw1t2w3(MM,:,:)).*phapad;
    
    if plotcont
        figure(1); clf(1); hold all;
        contourf(pw1,pw3,squeeze(real(sw1t2w3(MM,:,:))),10,'linecolor','k')
        set(gca,'box','on','layer','top')
        line(pw3,pw3,'color','k')
        daspect([1 1 1])
        xlabel('Excitation Frequency (eV)'); ylabel('Detection Frequency (eV)')
        title([dataname ', t_2 = ' num2str(t2(MM)) ' fs'])
        colorbar 
        set(gcf,'color','w')
        
        print([figfold num2str(t2(MM)) '_' dataname '.png'],'-dpng')
    end

    count = count + 1;
end

%% Save out fft-ed data

save([filename dataname '_fft.mat'],'sw1t2w3','pw1','pw3','t2')
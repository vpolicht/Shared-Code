userpath = '/Users/veronicapolicht/Documents/MATLAB';

clear all;
close all;
clc;
WhichOS;

%%

proj = 'eumelanin';
subproj = '';
twist = '';
date = '9-24-2020';
temp = 'RT';%'83 K';    
tdname = '2D_eumelanin-spot2-10ps-2';

calibname = '2DVIS_CALIB';

othermod = '';
binnum = '';
tlim = '';
avernum = '';
pumpsub = 0; % If the pump spectrum is subtracted enter 1
risetrack = 0;

xlm = [1.775 2.275];% in eV
ylm = [1.725 2.175]; % in eV
noilim = [1.8 2.3];%[1.88 2.35];

save2dmap = 1;
makegif = 0;
savepdf = 0;
probet2 = 100;
bnnum = 0;
zeropt = 2.04;%2.05; % Frequency to look for t2 = 0 
forcecalc = 0;
bgrplt = 0;
nolog = 1;

calibspec = ''; % If relevant load .mat calibrated probe wavelength array

%% Load the Excitaiton axis calibration

if ~isempty(subproj)
    datfolder = [basefolder sla proj sla subproj sla date];
else
    datfolder = [basefolder sla proj sla date];
end

figtitles = [proj ' - ' subproj ' T = ' temp];
 
if contains(proj,'HS')
    if contains(subproj,'WMS2')
        if contains(twist,'ua')
            m = 1;
        else
            m = 0;
        end
        
        if contains(temp,'RT')
            excitonev = excitonevset(1+2*m,:);
        else
            excitonev = excitonevset(2+2*m,:);
        end
        excitonlab = excitonlabset{1};
    end
elseif contains(proj,'ML')
    if contains(subproj,'MoS2')
        if contains(temp,'RT')
            excitonev = excitonevset(5,1:4);
        else
            excitonev = excitonevset(5,3:4);
        end
        excitonlab = excitonlabset{2};
    elseif contains(subproj,'WS2')
        if contains(temp,'RT')
            excitonev = excitonevset(6,1:4);
        else
            excitonev = excitonevset(6,3:4);
        end
        excitonlab = excitonlabset{3};
    end
elseif contains(proj,'V2O3')
    excitonlab = '';%excitonlabset{4};
    excitonev = 1;%excitonevset(7,:);
else
    excitonlab = {''};
    excitonev = [1];
end

if contains(proj,'HS')
    adl = ['-' twist];
else
    adl = '';
end

if exist([basefolder sla 'Linear Spectra' sla subproj adl '.txt'],'file')
tm = load([basefolder sla 'Linear Spectra' sla subproj adl '.txt']);
else
tm = zeros(1,2);
end
lin.spec = tm(:,2);
lin.wl = tm(:,1);
clear tm

if exist([datfolder sla 'temp'],'dir') && exist([datfolder sla 'temp' sla tdname '_' num2str(probet2) 'fs_' num2str(bnnum) '.bin'],'file')
tm = read_bin([datfolder sla 'temp' sla tdname '_' num2str(probet2) 'fs_' num2str(bnnum) '.bin']);
pr.l1 = tm(1,3:end);
pr.probe = tm(2,3:end);
clear tm
else
pr.l1 = '';
pr.probe = '';
end

if ~isempty(calibname)
    pat = [datfolder sla calibname '.dat'];
else
    pat = [datfolder sla '2DVIS_CALIB.dat'];
end

cal = dlmread(pat);

calib.dat = cal(4:end,3:end);
calib.phase = cal(3:end,2);
calib.l1 = cal(1,3:end);
calib.pump = cal(2,3:end);
calib.mstep = cal(3:end,1);
clear cal

if isempty(avernum)
savefold = [datfolder sla tdname];
else
savefold = [datfolder sla tdname '_aver' avernum];
end

if ~isempty(binnum)
savefold = [savefold '_bin' binnum];
end

if ~isempty(othermod)
savefold = [savefold '-' othermod];
end

if ~exist(savefold,'dir')
    mkdir(savefold)
end

figfold = [savefold sla 'Figures'];
if ~exist(figfold,'dir')
    mkdir(figfold)
end

load(['gradients.mat']);

%%
adl = '';

if isempty(avernum)
    datfile = [datfolder sla 'complete-' tdname '.mat'];
else
    datfile = [datfolder sla 'complete-' tdname '.mat'];
    adl = ['_aver' avernum];
end

if ~isempty(binnum)
    datfile = [datfolder sla 'complete-' tdname '_bin' binnum '.mat'];
    adl = [adl '_bin' binnum];
end

if ~isempty(othermod)
    datfile = [datfolder sla 'complete-' tdname '-' othermod '.mat'];
    adl = [adl '-' othermod];
end

if pumpsub
    inpt = strfind(datfile,'lete-');
    datfile = [datfile(1:inpt+4) 'PS-' datfile(inpt+5:end)];
    datfilem = [savefold sla 'complete-mod-PS-' tdname adl '.mat'];

    if isempty(adl)
    adl = strcat(adl,'_PS-');
    else
    adl = strcat(adl,'PS-');
    end
else
    datfilem = [savefold sla 'complete-mod-' tdname adl '.mat'];
end

if ~isempty(calibspec)
    load([datfolder sla calibspec '.mat']);
    l3calib = calwl;
    clear calib
    intp = 1;
else
    intp = 0;
end

%%
if ~exist(datfilem,'file') || forcecalc
    tddat = load(datfile);
    td.t2 = tddat.T2In;
    td.l3 = tddat.wprCutaver;
    
    td.l1 = tddat.wlpump;
    td.e1 = 1240./td.l1; 
    td.e3 = 1240./td.l3;
   
    td.data = zeros(size(tddat.B,1),size(tddat.B,2),length(td.t2));
    for m = 1:length(td.t2)
        td.data(:,:,m) = tddat.Map2DSave{m};
    end

    if ~isempty(zeropt)
        e1 = dsearchn(td.e1',zeropt);
        e3 = dsearchn(td.e3,zeropt);
    else
        [a e1] = max(max(max(td.data,3),1),2);
        [a e3] = max(max(td.data(:,e1,:),1),1);
    end

    if length(td.t2) > 2
        figure(1); clf(1); hold all;

        t2tr = squeeze(td.data(e1,e3,:));
        ttlim = dsearchn(td.t2,100);
        [A tstp] = max(abs(t2tr(1:ttlim)));
        A = A.*sign(t2tr(tstp));
        tstr = td.t2(tstp)-80;
        plot(td.t2,t2tr)
        t2tr = t2tr - mean(t2tr(1:dsearchn(td.t2,td.t2(tstp)-70)));
        plot(td.t2,t2tr)
        plot(td.t2(2:end),abs((diff(squeeze(td.data(e1,e3,:))))),'o')
        set(gca,'box','on','layer','top','linewidth',1.5,'fontname','arial')
        xline(0,'handlevisibility','off','linestyle',':')
        yline(0,'handlevisibility','off')
        xlabel('t_2 (fs)')
        set(gcf,'color','w')

        xline(td.t2(tstp),'r');
        xline(tstr,'g');
        outlie = excludedata(td.t2,t2tr,'domain',[tstr td.t2(tstp)]);

        ft = fit(td.t2(~outlie),t2tr(~outlie),'gauss1','robust','LAR',...
            'Upper',[A Inf 40],'Lower',[A -Inf 0]);
        fte = feval(ft,td.t2);
        fcoeff = coeffvalues(ft);

        fwhm = 2*sqrt(log(2))*fcoeff(3);
        t2zero = dsearchn(td.t2,fcoeff(2)-fwhm/2);
        xline(td.t2(t2zero),'b');
        xline(0,'linewidth',3);

        plot(td.t2,fte,'linewidth',1.5,'linestyle','--')

        figure(1)
        ttemp = td.t2-td.t2(t2zero);

        plot(ttemp,t2tr)

        xlim([-250 250])

        H = input('Does the new t2 = 0 look good? 1 - yes, 0 - no\n');
        while H == 0

            fprintf('Select best t2 = 0:\n')
            [l n] = ginput(1);
            t2zt = dsearchn(ttemp,l);
            t2zero = dsearchn(td.t2,ttemp(t2zt));
            ttemp = td.t2 - td.t2(t2zero);

            plot(td.t2-td.t2(t2zero),(squeeze(td.data(e1,e3,:))))
            H = input('Does the new t2 = 0 look good? 1 - yes, 0 - no\n');
        end
        clear b e1 e3

        td.t2 = td.t2 - td.t2(t2zero);
        td.data = permute(td.data,[2 1 3]);
    end
    
    size(td.data)
    clear tddat

    save([datfilem],'td');

    if ~isempty(tlim)
        tstp = dsearchn(td.t2,tlim);
        td.t2 = td.t2(1:tstp);
        td.data = td.data(:,:,1:tstp);
        adl = strcat(adl,[num2str(tlim) '-']);
    end
    
    t2 = td.t2;
    e1 = td.e1;
    e3 = td.e3;
    tddata = td.data;
    l1 = td.l1;
    l3 = td.l3;
    
    save([datfilem],'t2','e1','e3','tddata','l1','l3');
    
    clear t2 e1 e3 tddata l1 l3 
else
    dat = load([datfilem]);
    
    td.t2 = dat.t2;
    td.e1 = dat.e1;
    td.e3 = dat.e3;
    td.data = dat.tddata;
    td.l1 = dat.l1;
    td.l3 = dat.l3;
%     td.data = permute(td.data,[2 1 3]);

    
    clear dat
end

%% Select 2D plotting bounds
if isempty(xlm) || isempty(ylm)
    [xind yind] = deal(zeros(2,1));
    for m = 1:2;
        if m == 1; fprintf('Select upper-bounds for \omega_1 & \omega_3:\n');
        else fprintf('Select lower-bounds for \omega_1 & \omega_3:\n');
        end
        [a b] = ginput(1);
        xind(m) = dsearchn(td.e1',a);
        yind(m) = dsearchn(td.e3,a);
    end
    xind = (xind');
    yind = (yind');
    xlm = fliplr(td.e1(xind));
    ylm = fliplr(td.e3(yind));
else
    for m = 1:2
        xind(m) = dsearchn(td.e1',xlm(m));
        yind(m) = dsearchn(td.e3,ylm(m));
    end
    xind = fliplr(xind);
    yind = fliplr(yind);
end

if isempty(noilim)
    nxin = xind;
    nyin = yind;
else
    for m = 1:2
        nxin(m) = dsearchn(td.e1',noilim(m));
        nyin(m) = dsearchn(td.e3,noilim(m));
    end
    nxin = fliplr(nxin);
    nyin = fliplr(nyin);
end

Nline = 10;
levb = .1;
 
m = dsearchn(td.t2,150);
cmin = min(min(squeeze(td.data(:,:,m))));
cmax = max(max(squeeze(td.data(:,:,m))));
contf = 30;
contln = linspace(cmin,cmin*levb,Nline);
contlp = linspace(levb*cmax,cmax,Nline);
clear colors
colors = b2r(cmin,cmax);

figure(1); clf(1); hold all;
contourf(td.e1,td.e3,squeeze(td.data(:,:,m)),30,'linecolor','none');
colormap(colors); crange = caxis; 
contour(td.e1,td.e3,squeeze(td.data(:,:,m)),contln,'linecolor','k','linestyle','--');
contour(td.e1,td.e3,squeeze(td.data(:,:,m)),contlp,'linecolor','k','linestyle','-');
caxis(crange)
set(gca,'box','on','layer','top')%,'xgrid','on','ygrid','on')
ylabel('Detection frequency (eV)'); xlabel('Excitation frequency (eV)');
set(gcf,'color','w')
for M = 1:sum(~contains(excitonlab,'BGR'))
    xline(excitonev(M),'linestyle','--','handlevisibility','off');
    yline(excitonev(M),'linestyle','--','handlevisibility','off');
end
line(td.e1,td.e1,'color',[.9 .9 .9],'linewidth',1.5)
xlim([td.e3(end) td.e3(1)])
ylim([td.e1(end) td.e1(1)])
daspect([1 1 1])

text(1.7,2.35,['t_2 = ' num2str(td.t2(m)) ' fs'],'fontsize',12)
title(figtitles)
colorbar

%% Plot pulse spectra & reflectance spectrum

lmax = max(lin.spec(dsearchn(lin.wl,xlm(2)):dsearchn(lin.wl,xlm(1))));

if lmax == 0
    lmax = 1;
end

figure(1); clf(1); hold all;
a = area(lin.wl,smoothdata(lin.spec,'movmean',10),'linewidth',1.3);
a.FaceAlpha = .8;
a.FaceColor = [0.9290, 0.6940, 0.1250];
a.EdgeColor = [0.9290, 0.6940, 0.1250];
plot(1240./calib.l1,lmax.*calib.pump./max(calib.pump),'k','linewidth',1.5);
plot(1240./pr.l1,lmax.*pr.probe./max(pr.probe),'color',[0.6350, 0.0780, 0.1840],'linewidth',1.5);

if ~isempty(excitonlab)
for m = 1:sum(~contains(excitonlab,'BGR'))
    xline(excitonev(m),'linestyle','--','handlevisibility','off');
    text(excitonev(m)+.01,.8,excitonlab{m},'fontsize',18);
end
end

legend('Linear Spectrum','Pump spectrum from Calibration File','Probe Spectrum','location','best')
title(['Pump Spectrum, ' date])
ylabel([{'Normalized Int. (a.u.) - Pump'},{'OD - Spec'}])
xlabel('\omega (eV)')
set(gcf,'color','w')
set(gca,'box','on','layer','top','linewidth',1.5,'fontname','arial')
xlim(xlm)
% xlim([1.5 2.5])
print([figfold sla 'pump-spectrum.png'],'-dpng')
if savepdf
    print([figfold sla 'pump-spectrum.pdf'],'-dpdf')
end
saveas(gcf,[figfold sla 'pump-spectrum.fig'])

pump.e = 1240./calib.l1;
pump.spectrum = calib.pump;
probe.e = 1240./pr.l1;
probe.spectrum = pr.probe;
abso.e = lin.wl;
abso.spectrum = lin.spec;

save([savefold sla tdname adl '_linear-spectra.mat'],'pump','probe','abso')

%%
if savepdf; print([figfold sla 'pump-spectrum.pdf'],'-r300','-dpdf'); end

%% Plot 2D maps
% t2times = [0 10 50 100 200 1000];

if save2dmap
    for M = 1:length(td.t2) 
%         m = dsearchn(td.t2,t2times(M)');
        m = M
        cmin = min(min(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),m))));
        cmax = max(max(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),m))));
        contf = 30;
        contln = linspace(cmin,cmin*levb,Nline);
        contlp = linspace(levb*cmax,cmax,Nline);
        clear colors
        colors = b2r(cmin,cmax);

        figure(1); clf(1); hold all;
        contourf(td.e1(xind(1):xind(end)),td.e3(yind(1):yind(end)),squeeze(td.data(yind(1):yind(end),xind(1):xind(end),m)),30,'linecolor','none');
        colormap(colors); crange = caxis; 
        contour(td.e1(xind(1):xind(end)),td.e3(yind(1):yind(end)),squeeze(td.data(yind(1):yind(end),xind(1):xind(end),m)),contln,'linecolor','k','linestyle','--');
        contour(td.e1(xind(1):xind(end)),td.e3(yind(1):yind(end)),squeeze(td.data(yind(1):yind(end),xind(1):xind(end),m)),contlp,'linecolor','k','linestyle','-');
        caxis([cmin cmax])
        set(gca,'box','on','layer','top','linewidth',1.5,'gridcolor','k','fontname','arial')% 'xgrid','on','ygrid','on','yminorgrid','on',
        ylabel('Detection energy (eV)'); xlabel('Excitation energy (eV)');
        line(td.e1,td.e1,'color','k','linewidth',1.1);
        ylim([td.e3(fliplr(yind))])
        xlim([td.e1(fliplr(xind))])
        daspect([1 1 1])
        for n = 1:sum(~contains(excitonlab,'BGR'))
        xline(excitonev(n),'linestyle','--','color','k','linewidth',1);
        yline(excitonev(n),'linestyle','--','color','k','linewidth',1);
        end

        text(td.e1(xind(end))+.1,td.e3(yind(1))-.1,['t_2 = ' num2str(td.t2(m)) ' fs'],'fontsize',12)
        title(figtitles)
        colorbar
        set(gcf,'color','w')

        if save2dmap
        print([figfold sla tdname adl '_t2-' num2str(td.t2(m)) 'fs.png'],'-dpng')
    %     saveas(gcf,[figfold sla tdname adl '_t2-' num2str(td.t2(m)) 'fs.fig'])

        set(gcf, 'DefaultFigureRenderer', 'Painters');
        set(gcf,'Renderer','Painters');

        if savepdf
            print([figfold sla tdname adl '_t2-' num2str(td.t2(m)) 'fs.pdf'],'-r300','-dpdf')
        end

        end
    end
end
%% Make a gif of the t2 dependence
if makegif
tlm = dsearchn(td.t2,0);
cmin = min(min(min(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),tlm:end)))));
cmax = max(max(max(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),tlm:end)))));

gifname = [figfold sla tdname adl '_t2.gif'];

dt = .5;

writerObj = VideoWriter([figfold sla tdname adl '_t2.mp4'],'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);
tst = dsearchn(td.t2,-50);

for m = tst:length(td.t2)

    clear color
    colors = b2r(cmin,cmax);
    lmin = min(min(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),m))));
    lmax = max(max(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),m))));
    contln = linspace(lmin,lmin*levb,Nline);
    contlp = linspace(levb*lmax,lmax,Nline);
    
    
    figure(1); clf(1); hold all;
    contourf(td.e1(xind(1):xind(end)),td.e3(yind(1):yind(end)),squeeze(td.data(yind(1):yind(end),xind(1):xind(end),m)),30,'linecolor','none');
    colormap(colors); crange = caxis; 
    contour(td.e1(xind(1):xind(end)),td.e3(yind(1):yind(end)),squeeze(td.data(yind(1):yind(end),xind(1):xind(end),m)),contln,'linecolor','k','linestyle','--');
    contour(td.e1(xind(1):xind(end)),td.e3(yind(1):yind(end)),squeeze(td.data(yind(1):yind(end),xind(1):xind(end),m)),contlp,'linecolor','k','linestyle','-');
    caxis([cmin,cmax])
    set(gca,'box','on','layer','top','linewidth',1.5,'gridcolor','k','fontname','arial')%'xgrid','on','ygrid','on',
    ylabel('Detection energy (eV)'); xlabel('Excitation energy (eV)');
    line(td.e1,td.e1,'color','k','linewidth',1.1);
    ylim([td.e3(fliplr(yind))])
    xlim([td.e1(fliplr(xind))])
    daspect([1 1 1])
    for n = 1:sum(~contains(excitonlab,'BGR'))
    xline(excitonev(n),'linestyle','--','color','k','linewidth',1);
    yline(excitonev(n),'linestyle','--','color','k','linewidth',1);
    end
    
    text(td.e1(xind(end))+.1,td.e3(yind(1))-.1,['t_2 = ' num2str(td.t2(m)) ' fs'],'fontsize',12)
    title(figtitles)
    colorbar
    set(gcf,'color','w')
    
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if m == tst
       imwrite(imind,cm,gifname,'gif','Loopcount',inf,'DelayTime',dt);
    else
       imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',dt);
    end
    writeVideo(writerObj,frame);
    clear frame
end
close(writerObj);
end

%%
avsize = 5;

tr.t2 = td.t2;
tr.exset = excitonlab;
tr.exval = excitonev;
tr.dat = zeros(length(excitonlab(~contains(excitonlab,'BGR'))),...
    length(tr.exval),length(tr.t2));

for m = 1:length(excitonlab(~contains(excitonlab,'BGR')))
    for h = 1:length(tr.exval)
        xx = dsearchn(td.e1',tr.exval(m));
        yy = dsearchn(td.e3,tr.exval(h));
        
        if xx <= avsize
            xx = avsize+1;
        elseif xx == length(td.e1)
            xx = length(tr.t2)-1-avsize;
        end
        
        if yy <= avsize
            yy = avsize+1;
        elseif yy == length(td.e3)
            yy = length(tr.t2)-1-avsize;
        end
        
        tr.dat(m,h,:) = squeeze(mean(mean(td.data(yy-avsize:yy+avsize,xx-avsize:xx+avsize,:),1),2));
        tr.point(m,h,:) = [td.e1(xx) td.e3(yy)];
    end
end

save([savefold sla tdname '_t2-traces-easy-import' adl '.mat'],'tr')

%% Select 2D points and Plot against t2
avsize = 5;

% excitonlab = {'MoS2A','MoS2B','WS2A'};

% adlmset = {'-diag','-MS2A','-MS2B','-WS2A','-WS2B','-MS2B-rise','-xp-rise'};
adlmset = {'-diag'};
if ~bgrplt
    for g = 1:sum(~contains(excitonlab,'BGR'))
        adlmset{1+g} = ['-' excitonlab{g}];
    end
    exn = sum(~contains(excitonlab,'BGR'));   

else
    for g = 1:length(excitonlab)
        adlmset{1+g} = ['-' excitonlab{g}];
    end
    exn = length(excitonlab);   
end

if risetrack
    adlmset{size(adlmset,2)+1} = ['-' subproj '-rise'];
    adlmset{size(adlmset,2)+1} = ['-xp-rise'];
end

lenex = exn;
exval = cat(1,exn+1,linspace(1,exn,exn)',0,0)';
t2lms = cat(1,td.t2(end)*ones(1,dsearchn(exval',0)-1)',200,200)';
nmset = cat(1,td.t2(end)*zeros(1,dsearchn(exval',0)-1)',1,1)';

clear mst2tr

t2plt = td.t2;

if min(unique(diff(td.t2))) < max(unique(diff(td.t2)))*100 && ~nolog
    
    axsp = 'log';
    spc = min(diff(td.t2(td.t2 > 0 & td.t2 < 1000)));
    tmp = find(diff(td.t2) == spc);
    lnl = tmp(end) + 1;
    
else
    axsp = 'linear';
end

for h = 1:length(adlmset)
    
    clear t2tr
    
    if pumpsub == 1
        adlm = ['_PS' adlmset{h}];
    else
        adlm = adlmset{h};
    end

    fprintf(adlm)
    adl = adlm;
    close all;

    m = dsearchn(td.t2,250);
    cmin = min(min(min(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),m)))));
    cmax = max(max(max(squeeze(td.data(nyin(1):nyin(end),nxin(1):nxin(end),m)))));
    contf = 30;
    contln = linspace(cmin,cmin*levb,Nline);
    contlp = linspace(levb*cmax,cmax,Nline);
    clear colors
    colors = b2r(cmin,cmax);

    figure(1); clf(1); hold all;
    contourf(td.e1,td.e3,squeeze(td.data(:,:,m)),30,'linecolor','none');
    colormap(colors); crange = caxis; 
    contour(td.e1,td.e3,squeeze(td.data(:,:,m)),contln,'linecolor','k','linestyle','--');
    contour(td.e1,td.e3,squeeze(td.data(:,:,m)),contlp,'linecolor','k','linestyle','-');
    caxis([cmin,cmax])
    set(gca,'box','on','layer','top','linewidth',1.5,'fontname','arial')
    xlabel('Excitation energy (eV)'); ylabel('Detection energy (eV)');
    set(gcf,'color','w')
    line(td.e1,td.e1,'color',[.9 .9 .9],'linewidth',1.5)
    ylim([td.e3(fliplr(yind))])
    xlim([td.e1(fliplr(xind))])
    daspect([1 1 1])

        for n = 1:sum(~contains(excitonlab,'BGR'))
            xline(excitonev(n),'linestyle','--','color','k','linewidth',1);
            yline(excitonev(n),'linestyle','--','color','k','linewidth',1);
        end

    text(xlm(1)+.05,ylm(1)-.05,['t_2 = ' num2str(td.t2(m)) ' fs'],'fontsize',12)
    title(figtitles)
    colorbar

    colb = col.blbar;
    clear t2ind; 
    t2ind = {};
    count0 = 1;

    if exval(h) > exn

        exnx = linspace(1,lenex,lenex);
        exny = exnx;
        xnd = dsearchn(td.e1',excitonev(exnx)');
        ynd = dsearchn(td.e3,excitonev(exny)');
        for s = 1:length(xnd)
            t2ind{count0} = [xnd(s) ynd(s)]; 

            ab = plot(excitonev(exnx(s)),excitonev(exny(s)),'o','markersize',10); 
            ab.MarkerFaceColor = ab.Color;
            count0 = count0+1; 
        end    

    elseif exval(h) > 0 && exval(h) <= sum(~contains(excitonlab,'BGR'))  || bgrplt

        exnx = exval(h)*ones(1,lenex);
        exny = linspace(1,lenex,lenex);
        xnd = dsearchn(td.e1',excitonev(exnx)');
        ynd = dsearchn(td.e3,excitonev(exny)');
        for s = 1:length(xnd)
            t2ind{count0} = [xnd(s) ynd(s)]; 

            ab = plot(excitonev(exnx(s)),excitonev(exny(s)),'o','markersize',10); 
            ab.MarkerFaceColor = ab.Color;
            count0 = count0+1; 
        end

    elseif contains(adlm,'rise')

        fprintf('Select (w1,w3) positions horizontal limits t_2:\n');
        [a b] = ginput(2);
        xnd = dsearchn(td.e1',a);
        ynd = dsearchn(td.e3,b);
        trnum = floor(abs(diff(xnd))./(avsize+2));
            for m = 1:trnum
            xdiv = xnd(1)-(avsize+2)*(m-1);
            t2ind{m} = [xdiv,ynd(1)];
            ab = plot(td.e1(xdiv),b(1),'o','markersize',10); 
            ab.MarkerFaceColor = ab.Color;
            ab.MarkerFaceColor = colb(size(colb,1)-1*(m-1),:);
            end
        clear xnd ynd    

    else % Freeform pick (x,y) indices
        m = ''; count = 1;
        fprintf('Select (w1,w3) positions to plot as a function of t_2:\n');
        while isempty(m)
            [a b] = ginput(1);
            xnd = dsearchn(td.e1',a);
            ynd = dsearchn(td.e3,b);
            t2ind{count} = [xnd,ynd];
            ab = plot(a,b,'o','markersize',10); 
            ab.MarkerFaceColor = ab.Color;
            ab.MarkerFaceColor = colb(size(colb,1)-1*(count-1),:);
            m = input('Cancel?\n');
            clear xnd ynd
            count = count+1;
        end
    end

    print([figfold sla tdname adl '_t2-points.png'],'-dpng')

    set(0, 'DefaultFigureRenderer', 'painters');
    set(gcf,'Renderer','Painters');


    if savepdf
        print([figfold sla tdname adl '_t2-points.pdf'],'-r300','-dpdf')
    end

    
    adl = adlm;

    mx = nmset(h);
    sm = input('Enter smoothing value or "0"\n');
    smmeth = 'movmean';

    if isempty(sm)
        sm = 1;
    end
    
    tlim = t2lms(h);
    load('gradients.mat');
    clear legset
    count = 1;
    colb = col.blbar;

    figure(2); clf(2); hold all;

for m = 1:length(t2ind)
%     xx = t2ind{length(t2ind)-m+1}(1);
%     yy = t2ind{length(t2ind)-m+1}(2);
    xx = t2ind{m}(1);
    yy = t2ind{m}(2);
    
    [ypt xpt] = deal(1);
    if yy <= avsize
        yy = avsize+1;
        ypt = 0;
    elseif yy == length(td.e3)
        yy = length(td.e3)-1-avsize;
        ypt = 0;
    end
    
    if xx <= avsize
        xx = avsize+1;
        xpt = 0;
    elseif xx == length(td.e1)
        xx = length(td.e1)-1-avsize;
        xpt = 0;
    end
    
    if xpt && ypt
        
        t2tr = squeeze(mean(mean(td.data(yy-avsize:yy+avsize,xx-avsize:xx+avsize,:),1),2));
        
        if mx
            [a b] = max(abs(t2tr));
            sgn = sign(t2tr(b));
            t2tr = sgn.*t2tr./t2tr(b);
        end

        set(gca,'ColorOrderIndex',m);
        
        if contains(axsp,'log')
            
            if exval(h) == 0
                
                bc = plot(t2plt,t2tr,'linewidth',1,'color',colb(size(colb,1)-1*(m-1),:));
            
            else
                
                ymn(m) = floor(1.1*min(t2tr)*100)/100;
                ymx(m) = round(1.2*max(t2tr)*100)/100;

                if m == 1
                ax1 = subplot(121);
                hold all;
                ax2 = subplot(122);
                hold all;
                end
                
                axes(ax1)
                bc = plot(t2plt(1:lnl),t2tr(1:lnl),'o')%,'linewidth',1);
                axes(ax2)
                bc2 = plot(t2plt(lnl:end),t2tr(lnl:end),'o')%'linewidth',1);
                
                set(ax1,'units','normalized','position',[0.1 0.1 0.2 0.8]);
                set(ax1,'xlim',[t2plt(1) t2plt(lnl)])
                set(ax2,'units','normalized','position',[0.3 0.1 0.7 0.8]);
                set(ax2,'xscale','log','xlim',[t2plt(lnl) t2plt(end)],'yticklabel','');
            end

            legset{count} = ['(\omega_1, \omega_3) = (' num2str(round(td.e1(xx)*100)./100) ', ' num2str(round(td.e3(yy)*100)./100) ')'];

            if sm
                t0 = dsearchn(td.t2,0); tend = dsearchn(td.t2,200);
                [a,B] = max(abs(diff(diff(t2tr(t0:tend)))));
                B = B+1; 
                clear smd;    
                if B+t0-sm < 1
                    smst = B+t0-1;
                    ctst = 3;
                else
                    smst = B+t0-sm;
                    ctst = 2+sm;
                end
                smd = smoothdata(t2tr(smst:end),smmeth,sm);

                t2tr = cat(1,t2tr(1:B+t0),smd(ctst:end));
        %         t2tr = smoothdata(t2tr,'gaussian',sm);

                if mx
                    t2tr = sgn.*t2tr./t2tr(b);
                end

                bc.LineStyle = ':';
                bc.LineWidth = 1.5;
                bc2.LineStyle = ':';
                bc2 .LineWidth = 1.5;
                
                            
                if exval(h) == 0    
                    ab = plot(t2plt,t2tr,'linewidth',1,'color',colb(size(colb,1)-1*(m-1),:));
                else
                    axes(ax1)
                    ab = plot(t2plt(1:lnl),t2tr(1:lnl),'linewidth',1);
                    axes(ax2)
                    ab2 = plot(t2plt(lnl:end),t2tr(lnl:end),'linewidth',1)
                    axes(ax1)
                end

                ab.Color = bc.Color;
                ab2.Color = bc.Color;
                
                legset{count+1} = ['Smoothed (\omega_1, \omega_3) = (' num2str(round(td.e1(xx)*100)./100) ', ' num2str(round(td.e3(yy)*100)./100) ')'];
                count = count+1;
            end
            
        else
            
            ymn(m) = floor(.8*min(t2tr)*100)/100;
            ymx(m) = round(1.2*max(t2tr)*100)/100;
                
            if exval(h) == 0
                bc = plot(t2plt,t2tr,'linewidth',1,'color',colb(size(colb,1)-1*(m-1),:));
            else
                bc = plot(t2plt,t2tr,'linewidth',1);
            end

            legset{count} = ['(\omega_1, \omega_3) = (' num2str(round(td.e1(xx)*100)./100) ', ' num2str(round(td.e3(yy)*100)./100) ')'];

            if sm
                t0 = dsearchn(td.t2,0); tend = dsearchn(td.t2,200);
                [a,B] = max(abs(diff(diff(t2tr(t0:tend)))));
                B = B+1; 
                clear smd;    
                if B+t0-sm < 1
                    smst = B+t0-1;
                    ctst = 3;
                else
                    smst = B+t0-sm;
                    ctst = 2+sm;
                end
                smd = smoothdata(t2tr(smst:end),smmeth,sm);

                t2tr = cat(1,t2tr(1:B+t0),smd(ctst:end));
        %         t2tr = smoothdata(t2tr,'gaussian',sm);

                if mx
                    t2tr = sgn.*t2tr./t2tr(b);
                end

                bc.LineStyle = ':';
                bc.LineWidth = 1;

                if exval(h) == 0    
                    ab = plot(t2plt,t2tr,'linewidth',1,'color',colb(size(colb,1)-1*(m-1),:));
                else
                    ab = plot(t2plt,t2tr,'linewidth',1);
                end

                ab.Color = bc.Color;
                legset{count+1} = ['Smoothed (\omega_1, \omega_3) = (' num2str(round(td.e1(xx)*100)./100) ', ' num2str(round(td.e3(yy)*100)./100) ')'];
                count = count+1;
            end
            if isempty(pumpsub)
                xlim([t2plt(1) t2plt(dsearchn(t2plt,tlim))])
            else
                xlim([t2plt(1) t2plt(dsearchn(t2plt,tlim))])
            end
        end
        
        mst2tr{h,m,1} = adlm;
        mst2tr{h,m,2} = [td.e1(xx) td.e3(yy)];
        mst2tr{h,m,3} = td.t2;
        mst2tr{h,m,4} = squeeze(mean(mean(td.data(yy-avsize:yy+avsize,xx-avsize:xx+avsize,:),1),2));
        mst2tr{h,m,5} = t2tr;

        count = count+1;
    end
    
    tcks = linspace(min(ymn),max(ymx),11);
    
    avy = abs(td.e3(yy+avsize)-td.e3(yy-avsize));
    avx = abs(td.e1(xx+avsize)-td.e1(xx-avsize));
   

    
end

    if contains(axsp,'log')
        set([ax1 ax2],'ylim',[tcks(1) tcks(end)],'ytick',tcks,'box','off');
        set([ax1 ax2],'box','on','layer','top','linewidth',1.5,...
            'fontname','arial','xminortick','on');%'xgrid','on','ygrid','on','gridcolor','k'
        xline(0); yline(0); axes(ax2); yline(0);
        axes(ax1)
    else
        set(gca,'box','on','layer','top','linewidth',1.5,...
            'fontname','arial','xminortick','on');
        xline(0); yline(0);
    end
    
    if mx
        if sgn > 0
        ylim([-0.3 1.1])
        else
        ylim([-1.1 0.3])
        end
    adl = strcat(adl,'-norm');
    end
    if sm
    adl = strcat(adl,'-smoothed');
    end

    fg = gcf;
    set(fg,'Units','Centimeters');
    pos = fg.Position;
%     set(fg,'Position',[4 4 20 12]);
%     xlim([-150 950])
    xlabel('t_2 (fs)'); ylabel('2D Amplitude (a.u.)');
    set(gcf,'color','w')
    
    if contains(axsp,'log')
        axes(ax2)
    end
    
    legend(legset,'location','northeastoutside')
    title([{'t_2 traces of selected peaks'},{['Averaged over ' num2str(round(mean([avy,avx]).*100)/100) ' eV']}])

    print([figfold sla tdname adl '_t2-traces.png'],'-dpng')

%     saveas(gcf,[figfold sla tdname adl '_t2-traces.fig'])

    set(0, 'DefaultFigureRenderer', 'painters');
    set(gcf,'Renderer','Painters');

    if savepdf
        print([figfold sla tdname adl '_t2-traces.pdf'],'-r300','-dpdf'); 
    end

end

save([savefold sla tdname '_t2-traces.mat'],'mst2tr','avsize')
%%
clear t2ind

count = 1;
for m = 1:length(excitonev)
    e3ind = dsearchn(td.e3,excitonev(m));
    for n = 1:length(excitonev)
    e1ind = dsearchn(td.e1',excitonev(n));   
    t2ind{count} = [e1ind,e3ind];
    ptname{count} = strcat(excitonlab{n},excitonlab{m});
    count = count+1;
    clear e1ind
    end
    clear e3ind
end

nm = @(x) x./max(x);

[smpk pkdat] = deal(zeros(length(ptname),length(td.t2)));
figure(3); clf(3); hold all;
count = 1;
for m = 1:length(ptname)
    clear B t0 tend smd
    
    e1 = t2ind{m}(1);
    e3 = t2ind{m}(2);
    
    if e1 == 1
        e1 = avsize+1;
    elseif e1 == length(td.e1)
        e1 = length(td.e1)-avsize+1;
    end
    
    if e3 == 1
        e3 = avsize+1;
    elseif e3 == length(td.e3)
        e3 = length(td.e3)-avsize+1;
    end

    pkdat(m,:) = squeeze(mean(mean((td.data(e3-avsize:e3+avsize,e1-avsize:e1+avsize,:)),1),2));
    
    t0 = dsearchn(td.t2,0); tend = dsearchn(td.t2,200);
    [a,B] = max(abs(diff(diff(pkdat(m,t0:tend)))));
    B = B+1; 
    if B+t0-sm < 1
        smst = B+t0-1;
        ctst = 3;
    else
        smst = B+t0-sm;
        ctst = 2+sm;
    end
    smd = smoothdata(pkdat(m,smst:end),smmeth,5);

    smpk(m,:) = cat(2,pkdat(m,1:B+t0),smd(ctst:end));
        
%     smpk(m,:) = smoothdata(pkdat(m,:),smmeth);
    a = plot(td.t2,nm(pkdat(m,:)),'linewidth',2,'linestyle',':');
    b = plot(td.t2,nm(smpk(m,:)),'linewidth',1.5);
    b.Color = a.Color;
    legname{count} = ['\nu_1 = ',num2str(round(100*td.e1(e1))/100),'eV, \nu_3 = ',num2str(round(100*td.e3(e3))/100),' eV'];
    legname{count+1} = ['Smoothed - \nu_1 = ',num2str(round(100*td.e1(e1))/100),' eV, \nu_3 = ',num2str(round(100*td.e3(e3))/100),' eV'];
    count = count + 2;
end
fg = gcf;
set(fg,'Units','Centimeters');
pos = fg.Position;
set(fg,'Position',[4 4 20 12]);
xlim([td.t2(1) td.t2(dsearchn(td.t2,tlim))])
avy = abs(td.e3(e3+avsize)-td.e3(e3-avsize));
avx = abs(td.e1(e1+avsize)-td.e1(e1-avsize));
legend(legname,'location','northeastoutside')
title([{'t_2-dependent diagonal & cross peaks'},{['Averaged over ' num2str(round(mean([avy,avx]).*100)/100) ' eV']}])
xlabel('t_2 (fs)')
ylabel('2D Amplitude (a.u.)')
set(gca,'box','on','layer','top','linewidth',1.5,'fontname','arial')
set(gcf,'color','w')

print([figfold sla adl 'pk-decays.png'],'-dpng')
% saveas(gcf,[figfold sla adl 'pk-decays.fig'])
if savepdf; print([figfold sla adl 'pk-decays.pdf'],'-r300','-dpdf'); end

%%

[smpk pkdat] = deal(zeros(length(ptname),length(td.t2)));
figure(3); clf(3);
subplot(1,2,1); hold all;
for m = 1:length(ptname)
    e1 = t2ind{m}(1);
    e3 = t2ind{m}(2);
    pkdat(m,:) = squeeze(mean(mean((td.data(e3-avsize:e3+avsize,e1-avsize:e1+avsize,:)),1),2));
    smpk(m,:) = smoothdata(pkdat(m,:),smmeth);
    a = plot(td.t2,pkdat(m,:),'linewidth',1);
%     plot(td.t2,smpk(m,:),'HandleVisibility','off','color',)
   legname{m} = strcat(ptname{m},' - (',num2str(td.e1(e1)),', ',num2str(td.e3(e3)),')') 
end
legend(legname,'location','best')
title([{'t_2-dependent diagonal & cross peaks'},{'Averaged over 0.1 eV'}])
xlabel('t_2 (fs)')
ylabel('2D Amplitude (a.u.)')
set(gca,'box','on','layer','top','linewidth',1.5,'fontname','arial')

subplot(1,2,2); hold all;
plot(td.t2,pkdat(3,:)./smpk(1,:),'color',[.8 .8 .8],'linewidth',1,'linestyle','--')
plot(td.t2,pkdat(2,:)./pkdat(4,:),'color',[.5 .5 .5],'linewidth',1)
plot(td.t2,smoothdata(pkdat(3,:)./smpk(1,:),smmeth),'linewidth',1.5)
plot(td.t2,smoothdata(pkdat(2,:)./pkdat(4,:),smmeth),'linewidth',1.5)
xlim([0 td.t2(end)]); ylim([0 1.5]);
legend('AB/AA','BA/BB','AB/AA - smoothed','BA/BB - smoothed','location','best')
xlabel('t_2 (fs)')
ylabel('Ratio (a.u.)')
title(['t_2-dependent ratio of diagonal & cross peaks'])

fg = gcf;
set(fg,'Units','Centimeters');
pos = fg.Position;
set(fg,'Position',[4 4 30 18]);

set(gca,'box','on')%,'xgrid','on','ygrid','on','yminorgrid','on','yminortick','on')

set(gcf,'color','w')
sgtitle([proj ', ' subproj ' - ' temp ', ' date]) 
print([figfold sla adl 'pk-ratio-comp.png'],'-dpng')

if savepdf; print([figfold sla adl 'pk-ratio-comp.pdf'],'-r300','-dpdf'); end

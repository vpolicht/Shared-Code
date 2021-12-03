% This script takes the global fit calculated using run_global_fit.mat %
% and subtracts background kinetics from phased data, leaving the residual
% which is then fourier transformed over t2 to yeild the 3D frequency
% solid. Frobenius spectra is calculated and peak frequencies are found. %
% This script is a modified version of FrobeniusSpectrumRepshasing_v3.mat
% by S Senlik. V Policht 09/07/2016

% Second to run

close all; clear all; clc; WhichOS;

% User input required for defining: 
% 1 Experimental set; 2 windows. %

%% 1 Experimental set

proj = 'V2O3';
subproj = '';
twist = '';
date = '7-3-2020';
temp = '100 K';%'83 K';    
tdname = '2D_V2O3-spot1-100k-1ps-1';
fit = '2020-7-3_17-14-abs';

othermod = '';
binnum = '';
avernum = '';

t2trunc = '';
t2st = '80';
t2end = '';

forcecalc = 1;

excitonwn = [2.32];%, 1.95];%[1.88 2.03];

t2trpt = '';

xlm = [2.1 2.53];% in eV
ylm = [2.23  2.53]; % in eV

% If want to reduce size of 2D matrix, enter number other than 1
spacenum = 1;

% Recalculate the fft
forcecalc = 1;

%% Begin code

FitFolder = fit;

if ~isempty(t2trunc)
    if isempty(t2st)
    t2txt = [num2str(t2trunc) '-'];
    else
    t2txt = [num2str(t2st) '-'];
    end
else
t2txt = [''];
end

datfolder = [basefolder sla proj sla subproj sla date];
loadname = [datfolder sla tdname sla 'complete-mod-' tdname '.mat'];
prfx = 'abs_';
oloadname = [datfolder sla tdname sla 'compressed_complete-mod-' tdname '.mat'];
fitresfld = [datfolder sla tdname sla 'FittingResults'];
fitfld = [fitresfld sla FitFolder];
resdat = [fitresfld sla FitFolder sla  t2txt prfx 'bckgfree.mat'];

if exist(resdat,'file') && ~forcecalc
    fprintf('Loading background-free data...')
    remat = load(resdat);
    
    e1 = remat.e1;
    e3 = remat.e3;
    T = remat.T;
    Zff = remat.data;
    
    srx = length(e3);
    sry = length(e1);
    
    nfft = 2^10;
    dT = mean(abs(diff(T)));
    f_fft_fs = ((1/dT)/2)*linspace(0,1,nfft/2 +1);
    k_fft_cm = f_fft_fs./(sol*1e-7);
    
    fprintf('loaded.\n')
    
else

    fprintf('Loading global fit data...\n')
    load([fitresfld sla FitFolder sla 'best_model.mat']);
    fprintf('Loading data...')
    load(loadname);
    fprintf('Data loaded.\n')

    if isempty(spacenum)
        spacenum = 1;
    end

    w1_indlb = dsearchn(e1',xlm(2)) ;
    w1_indup = dsearchn(e1',xlm(1));
    w3_indub = dsearchn(e3,ylm(1));
    w3_indlb = dsearchn(e3,ylm(2));

    x_axis_w_roi = e1;
    y_axis_w_roi = e3;

    e1 = e1(w1_indlb:w1_indup);
    e3 = e3(w3_indlb:w3_indub);

    T2vec = t2;
    clear t2 
    
    if isempty(t2trunc)
        t2trunc = 0;
    end

    lonum = t2trunc - 5; upnum = t2trunc + 5;
    
    if ~isempty(t2st)
        lonum = t2st - 5; upnum = t2st + 5;
    end
    
    if ~isempty(t2end)
        lonume = t2end - 5; upnume = t2end + 5;
    end

    if isempty(t2end)
        Trange = find(T2vec>lonum&T2vec<upnum):length(T2vec);
    else
        Trange = find(T2vec>lonum&T2vec<upnum):find(T2vec>lonume&T2vec<upnume);
    end

    T = T2vec(Trange);  

    all_data_mat_r = tddata; 
    clear tddata

    if size(size(all_data_mat_r),2) < 3
        all_data_mat_r = reshape(all_data_mat_r,length(y_axis_w_roi),length(x_axis_w_roi),length(T2vec));
    end 

    all_data_mat_r = squeeze(all_data_mat_r(w3_indlb:w3_indub,w1_indlb:w1_indup,Trange));

    [srx,sry,srz] = size(all_data_mat_r);
    
    if length(unique(diff(T))) > 1
        smp = min(unique(diff(T)));
        Tintp = linspace(T(1),T(end),(T(end)-T(1))/smp + 1);
        intp = 1;
    else
        intp = 0;
    end
    
    P = best_model.modelfun(best_model.Fit.beta);

    X = vectorize(all_data_mat_r);

    [b a] = max(max((real(X))));
    [c b] = max((real(X(:,a))));

    fprintf('Begin subtraction of kinetics...\n')
    Pf = multi_cexp2(best_model.Fit.beta(2*best_model.Nosc+1:end),0,best_model.StatRates,best_model.ConstBool,T);

    % Define linear least-squares solution to Pf*af = X'
    af = Pf\X';
    Zf = X-(Pf*af)'; %the residual without oscillators

    ffolder = [fitfld sla];

    pf = multi_cexp2(best_model.Fit.beta,0,best_model.StatRates,best_model.ConstBool,T);
    Af = pf\X';
    
    filter_f = @(L,a) [1; subsref(tukeywin(2*L,a),struct('type','()','subs',{{(L+2):(2*L),1}}))];
    filter_mat = repmat(filter_f(srz,0.3)',[srx*sry 1]);

    Zff = filter_mat.*Zf;

    data = Zff;
    
    if intp
        tdata = zeros(size(data,1),length(Tintp));
        
        for m = 1:size(data,1)
        tdata(m,:) = interp1(T,data(m,:),Tintp,'linear','extrap');
        end
        
        T = Tintp;
        [data Zff] = deal(tdata);
    end

    figure; hold all;
    plot(T,data(1000,:));
    
    fprintf('Saving residual data...')
    save([fitfld sla  t2txt prfx 'bckgfree.mat'],'e1','e3','T','data','-v7.3');
    fprintf('Saved.\n')
    
    clear T data e1 e3
    
    remat = load(resdat);
    
    e1 = remat.e1;
    e3 = remat.e3;
    T = remat.T;
    
    srx = length(e3);
    sry = length(e1);
    
    nfft = 2^10;
    dT = mean(abs(diff(T)));
    f_fft_fs = ((1/dT)/2)*linspace(0,1,nfft/2 + 1);
    k_fft_cm = f_fft_fs./(sol*1e-7);
end

%% Calculate Frobenius spectrum

if exist([fitfld sla  t2txt prfx 'bckgfree_full_fft.mat'],'file') && ~forcecalc
    load([fitfld sla  t2txt prfx 'bckgfree_full_fft.mat']);
    load([fitfld sla t2txt prfx 'plotting_set_peaks.mat'])
    
else
    fprintf('Perform fourier transform over t2:\n')

    s1 = size(Zff,1)/32;

    for ii=1:32
    Zffx = Zff((ii-1)*s1+1:s1*ii,:);
    fZfx = fft(Zffx,nfft,2);
    fZf_sum = sum(abs(fZfx).^2,1);
    fZf_total(ii,:) = fZf_sum;
    end
    clear Zffx fZfx fZf_sum

    grandsum = (sum((fZf_total),1))';
    Sfro = sqrt(grandsum);
    clear grandsum fZf_total 

    frobenius_spec = Sfro(1:nfft/2+1); 
    normspec = (frobenius_spec)./sum(frobenius_spec);

    fac = .30;
    [bkg bwn] = findpeaks(-1.*normspec);
    bkg = -1.*bkg;
    bkg = interp1(k_fft_cm(bwn),bkg,k_fft_cm);
    tmp = normspec - smooth(bkg,.2,'rloess');
    figure; hold all
    plot(k_fft_cm,normspec)
    plot(k_fft_cm,tmp)
    plot(k_fft_cm,smooth(bkg,.2,'rloess'))

    minh = fac.*(max(tmp));

    [amps inds] = findpeaks(tmp,'MinPeakHeight',minh,'SortStr','descend');
    modes = k_fft_cm(inds)';
    plot(modes,amps,'kd')

    % Plot frobenius spectrum

    frobspec = figure;
    plot(k_fft_cm,normspec);
    xlim([0 k_fft_cm(end)]);
    xlabel('\omega_{2}')
    ylabel('Intensity (a.u.)');

    fpks = figure;
    hold all
    plot(k_fft_cm,normspec)
    plot(modes,normspec(inds),'kd')
    plot(get(gca,'XLim'),[minh minh])

    for m = 1:length(inds)
        text(k_fft_cm(inds(m))-10,normspec(inds(m))+.05*max(normspec),num2str((round(k_fft_cm(inds(m)))))); 
    end

    xlabel('\omega_2 (cm^{1})')
    ylabel('Intensity (a.u.)');

    saveas(frobspec,[fitfld sla t2txt prfx 'frobenius-spec.fig'])
    saveas(frobspec,[fitfld sla t2txt prfx 'frobenius-spec.png'])

    saveas(fpks,[fitfld sla  t2txt prfx 'frobenius-peaks.fig'])
    saveas(fpks,[fitfld sla  t2txt prfx 'frobenius-peaks.png'])

    save([fitfld sla t2txt prfx 'plotting_set_peaks.mat'],'inds','modes');
    dlmwrite([fitfld sla t2txt prfx 'frob.txt'],normspec);

    save([fitfld sla  t2txt prfx 'frobspec_and_peaks.mat'],'k_fft_cm','normspec','frobenius_spec','modes','inds')
    clear frobenius_spec Sfro

    size(Zff)

    %% Will be calculating and saving out the 3D matrix in chunks

    vex = sry*srx;
    num = 3;

    % Number of chunks
    if length(find(factor(vex) == 2)) >= num
        Nbin = 2^num;
    else
        a = factor(vex);
        Nbin = a(find(factor(vex) > 3,1));
    end

    Nlength = vex/Nbin;

    % Define matrices to be saved into; have save abs and angle separatly
    % because cannot save a complex number into a real matrix and can't define
    % a complex matrix for matfile...

    save([fitfld sla  t2txt prfx 'bckgfree_full_fft.mat'],'k_fft_cm','e1','e3','-v7.3');
    % Allow for writeable MatFile
    r = matfile([fitfld sla  t2txt prfx 'bckgfree_full_fft.mat'],'Writable',true);

    %% Perform FFT in chunks along w1-w3 vector
    clear X Y

    tic
    fprintf(1,'Calculate full 3D solid using %g bins...',Nbin)

    for s = 1:Nbin

        arange = (s*Nlength-Nlength)+1:s*Nlength;
        data = remat.data(arange,:);

        parfor m = 1:Nlength

            X(m,:) = fft(data(m,:),nfft,2);

        end

        clear data;

        r.full_fft(arange,1:nfft/2+1) = X(:,1:nfft/2+1);

        fprintf('bin %g saved in %g sec...\n',s,round(toc))
        clear X

    end

    fprintf('3D matrix calculated and saved in %g minutes.\n',round(toc/60))
end

%%

load([fitfld sla  t2txt prfx 'bckgfree_full_fft.mat'])

srx = length(e3); sry = length(e1);

for m = 1:2%:length(modes)
    figure(1); clf(1);
    hold all;
    
    tmp = reshape(full_fft(:,inds(m)),srx,sry);
    
    contourf(e1,e3,abs(tmp),30,'linecolor','none')
    contour(e1,e3,abs(tmp),10,'linecolor','k')
    plot(e1,e1,'k','linewidth',1.5)
%     plot(e1,e1 - wn2eV(modes(m)),'k','linewidth',1.5,'linestyle','--')
%     plot(e1,e1 + wn2eV(modes(m)),'k','linewidth',1.5,'linestyle','--')
%     
    for n = 1:length(excitonwn)
        xline(excitonwn(n),'k','linestyle','--','linewidth',2)
        yline(excitonwn(n),'k','linestyle','--','linewidth',2)
    end
    
    
    xlim([xlm])
    ylim([ylm])
    daspect([1 1 1])
    
    set(gca,'layer','top','box','on','linewidth',2)
    xlabel('e_1 (eV)')
    ylabel('e_3 (eV)')
    
    set(gcf,'color','w')
    title(['\omega_2 = ' num2str(round(modes(m))) ' cm^{-1}'])
    
    fignm = [fitfld sla 'cohmap-' num2str(round(modes(m))) 'wn'];
    print([fignm '.png'],'-dpng')
    saveas(gcf,[fignm '.fig'])
    
    clear tmp
end

%%

for m = 1:2%:length(modes)
    figure(1); clf(1);
    hold all;
    
    tmp = reshape(full_fft(:,inds(m)),srx,sry);
    
    contourf(e1,e3,angle(tmp),30,'linecolor','none')
    contour(e1,e3,angle(tmp),10,'linecolor','k')
    plot(e1,e1,'k','linewidth',1.5)
%     plot(e1,e1 - wn2eV(modes(m)),'k','linewidth',1.5,'linestyle','--')
%     plot(e1,e1 + wn2eV(modes(m)),'k','linewidth',1.5,'linestyle','--')
%     
    for n = 1:length(excitonwn)
        xline(excitonwn(n),'k','linestyle','--','linewidth',2)
        yline(excitonwn(n),'k','linestyle','--','linewidth',2)
    end
    
    
    xlim([xlm])
    ylim([ylm])
    daspect([1 1 1])
    
    set(gca,'layer','top','box','on','linewidth',2)
    xlabel('e_1 (eV)')
    ylabel('e_3 (eV)')
    
    set(gcf,'color','w')
    title(['Phase \omega_2 = ' num2str(round(modes(m))) ' cm^{-1}'])
    
    fignm = [fitfld sla 'phmap-' num2str(round(modes(m))) 'wn'];
    print([fignm '.png'],'-dpng')
    saveas(gcf,[fignm '.fig'])
    
    clear tmp
end

% Use this script to calculate a global fit to the 2D data set %
% This script is a modified version of fit_r_nr_Bchla_15_06_08.m %
% from S Senlik in 2015. This script loads the rephasing (with the %
% option of loading the nonrephasing) data, windows the data either %
% with user-defined or predefined and loaded parameters, compresses %
% the data via projection to orthogonal polynomials, and runs a %
% global analysis fit using complex exponentials. V Policht 09/07/2016 %

% First to run in sequence.

close all; clear all; clc; WhichOS;

% User input required for defining: 
% 1 Experimental set; 2 windows; 3 fitting parameters. %

%% 1 Experimental set

proj = 'PGFROG';
subproj = '';
twist = '';
date = '3-5-2020';
temp = 'RT';%'83 K';    
tdname = '2D_V2O3-spot1-100k-1ps-1';

othermod = '';
binnum = '';
avernum = '';

t2trunc = '';
t2end = '';

forcecalc = 1;

excitonwn = 2.32;%[1.88 2.03];

t2trpt = '';

xlm = [2.1 2.53];% in eV
ylm = [2.23  2.53]; % in eV

%% 3 Fitting parameters

% Need to change these, are specific to the sample studied %~~~~~~~~~~~~~~~

FitID.RateBounds.rub = [1/10 1/500];
FitID.RateBounds.rlb = [1/300 1/10000];
% Define one rate that the fitting must find
FitID.MaxConstRates = [1/100000];
FitID.MaxNOsc = 4;
FitID.MinNOsc = 1;
FitID.MaxNdynRates = 3;
FitID.MinNdynRates = 1;

%% Should not to edit code below ------------------------------------------

datfolder = [basefolder sla proj sla subproj sla date];
adl = '';

if isempty(avernum)
datfile = [datfolder sla 'complete-' tdname '.mat'];
else
datfile = [datfolder sla 'complete-' tdname '.mat'];
adl = ['_aver' avernum '-'];
end

if ~isempty(binnum)
    datfile = [datfolder sla 'complete-' tdname '_bin' binnum '.mat'];
    adl = [adl '_bin' binnum];
end

if ~isempty(othermod)
    datfile = [datfolder sla 'complete-' tdname '-' othermod '.mat'];
    adl = [adl '-' othermod];
end

loadname = [datfolder sla tdname sla 'complete-mod-' tdname '.mat'];

prfx = 'abs_';
oloadname = [datfolder sla tdname sla 'compressed_complete-mod-' tdname '.mat'];
 
fitresfld = [datfolder sla tdname sla 'FittingResults'];

if ~exist(fitresfld,'dir')
    mkdir(fitresfld);
end

if exist([oloadname],'file') && ~forcecalc
    
    fprintf('Loading compressed data...')
    load([oloadname])
    fprintf('Loaded.\n')
    
elseif ~exist(oloadname,'file') || forcecalc
    fprintf('Load rephasing data...')
    load(loadname);
    fprintf('loaded.\n')
  
    w1_indlb = dsearchn(e1',xlm(2)) ;
    w1_indup = dsearchn(e1',xlm(1));
    w3_indub = dsearchn(e3,ylm(1));
    w3_indlb = dsearchn(e3,ylm(2));

    x_axis_w_roi = e1;
    y_axis_w_roi = e3;
    T2vec = t2;
    clear t2 
    
    if isempty(t2trunc)
        t2trunc = 80;
    end

    lonum = t2trunc - 7; upnum = t2trunc + 7;
    if ~isempty(t2end)
        lonume = t2end - 7; upnume = t2end + 7;
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

    [SRX, SRY, SRZ] = size(all_data_mat_r);
    mom = factor(SRX*SRY);

    if length(find(mom == 2)) <= 1 && length(find(mom < 20 & mom ~= 2)) < 1;
        input('Try redefining 2D limits; will not factor well...\n')
    end

    %%
    fprintf('Compress rephasing data...')

    Nbf = 60;
    Bx = dop( size(all_data_mat_r,2), Nbf );
    By = dop( size(all_data_mat_r,1), Nbf );

    Rcompressed = zeros(Nbf,Nbf,SRZ);
    for k=1:SRZ
      Rcompressed(:,:,k) = By' * squeeze(all_data_mat_r(:,:,k)) * Bx; 
    end
    R = zeros(Nbf*Nbf,SRZ);
    for k=1:Nbf
      R(((k-1)*Nbf+1):(k*Nbf),:) = squeeze(Rcompressed(:,k,:));
    end

    [Ur, Sr, Vr] = svd(R,0);
    Yr = (Sr*Vr')';

    fprintf(1,'Compressed.\n')

    fprintf(1,'Compressed rephasing data is being saved...')
    save([oloadname],'Rcompressed','Yr','Ur','Sr','Vr','Bx','By','Nbf','T','T2vec','Trange');
    fprintf('Saved.\n');

end

Y = Yr;

%%
fprintf(1,'Prepare Nonlinear Fitting...\n')

% Do not need to change these, are defined by scanning parameters %~~~~~~~~

% Bound the frequency of oscillations that can be found:
FitID.fub = 1/(2*mean(diff(T))); % nyquist frequency
% FitID.fub = 0.0060; % Low frequency only (<200 wn)
% FitID.flb = 1/(T(end)/(0.75*2*pi)); %require at least 3/4 period
FitID.flb = 1/(T(end)/(2*pi)); %require at least 1 full period
% Bound the rates of the Exponentials that can be found:
FitID.rub = 1/10;
FitID.rlb = 1/(T(end));
FitID.gen_rand_beta = @(lb,ub) ((ub-lb).*rand(length(ub),1)+lb);
% Oscillation specifics:
FitID.OscBounds.fub = FitID.fub;
FitID.OscBounds.flb = FitID.flb;
% Damping rates of oscillations
FitID.OscBounds.rub = FitID.rub;
FitID.OscBounds.rlb = FitID.rlb;

%All the possible models
[ModelCell, Diagnostics] = GenManySoLModels(FitID.MaxNOsc, FitID.MinNOsc, FitID.MaxNdynRates, FitID.MinNdynRates, FitID.MaxConstRates, FitID.OscBounds, FitID.RateBounds, T);

%constant fitting parameters
FitID.local_tol = 1e-6;
FitID.wall_clock_time = 10; %seconds
FitID.global_method = NLOPT_G_MLSL_LDS;
FitID.Ntests = 1;
tic

%% Fit Loop
fprintf(1,'Nonlinear Fitting has begun...\n')

parfor k = 1:length(ModelCell)
    
    beta_start = FitID.gen_rand_beta(ModelCell{k}.beta_lb,ModelCell{k}.beta_ub);
    
    Nosc = ModelCell{k}.Nosc;
    %beta_start(2:2:(Nosc*2)) = 1./T2wn(candidate_freqs(randi([1,length(candidate_freqs)],[1,Nosc]))); %#ok<PFBNS>
    
    ada = ModelCell{k}.modelfun
    
    [ModelCell{k}.Fit,ModelCell{k}.Stats,ModelCell{k}.exitcode] =  ...
        varpro_multiway_gradient3(Y,beta_start,ada, ...
        ModelCell{k}.beta_lb,ModelCell{k}.beta_ub,FitID.global_method,...
        FitID.wall_clock_time,FitID.local_tol,NLOPT_LN_SBPLX);
                        if (ModelCell{k}.exitcode <= 5);ModelCell{k}.exitcode;end
                        ModelCell{k}.exitcode
end

a = clock;
s = [num2str(a(1)) '-' num2str(a(2)) '-' num2str(a(3)) '_' num2str(a(4)) '-' num2str(a(5))];

toc
%% Refine Model
tic
aicvec = zeros(1,length(ModelCell)); for k=1:length(ModelCell); aicvec(k) = ModelCell{k}.Stats.AICc; end;
[minval, minind] = min(aicvec);
best_model = ModelCell{minind};
uber_model = best_model;
ada = best_model.modelfun;
[uber_model.Fit,uber_model.Stats,uber_model.exitcode] =  ...
        varpro_multiway_gradient3(Y,best_model.Fit.beta,ada, ...
                            best_model.beta_lb,best_model.beta_ub,FitID.global_method,FitID.wall_clock_time,1e-10,NLOPT_LD_LBFGS);
best_model = uber_model;
uber_model.Stats.relerr

toc
%% print some results to screen
fprintf(1,'------- Results ------- \n');

fprintf(1,'Oscillator Freqs\n');
fprintf(1,'%3.3g\n',T2wn(1./uber_model.Fit.beta(2:2:(2*uber_model.Nosc))));
fprintf(1,'Kinetic Lifetimes\n');
fprintf(1,'%3.3g\n',1./uber_model.Fit.beta((2*uber_model.Nosc+1):end));
toc

fitfld = [fitresfld sla s '-abs'];
mkdir([fitresfld sla s '-abs']);

save(fullfile(fitfld,['FitID' '.mat']),'FitID','-v7.3');
save(fullfile(fitfld,['ModelCell' '.mat']),'ModelCell','-v7.3');
save(fullfile(fitfld,['best_model' '.mat']),'best_model','-v7.3');
save([fitfld sla 'Trange.mat'],'T','T2vec','Trange');

%% 

load(loadname);

w1_indlb = dsearchn(e1',xlm(2)) ;
w1_indup = dsearchn(e1',xlm(1));
w3_indub = dsearchn(e3,ylm(1));
w3_indlb = dsearchn(e3,ylm(2));

x_axis_w_roi = e1;
y_axis_w_roi = e3;

re1 = e1(w1_indlb:w1_indup);
re3 = e3(w3_indlb:w3_indub);

T2vec = t2;
clear t2 

if isempty(t2trunc)
    t2trunc = 0;
end

lonum = t2trunc - 7; upnum = t2trunc + 7;
if ~isempty(t2end)
    lonume = t2end - 7; upnume = t2end + 7;
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

%%
pf = multi_cexp2(best_model.Fit.beta,0,best_model.StatRates,best_model.ConstBool,T);
Af = pf\X';

rateoscmat = reshape(Af',srx,sry,size(Af,1));

tic
% NumberOfRates = size(rateoscmat,3);
% DynamicRates = best_model.Fit.beta((2*best_model.Nosc+1):end); %check for the case when Nosc different from zero
StaticRate = best_model.StatRates;

rateosc = zeros(size(Af,1),1);
rateosc(1:2:2*best_model.Nosc-1) = 1./best_model.Fit.beta(1:2:2*best_model.Nosc-1);
rateosc(2:2:2*best_model.Nosc) = T2wn(1./(best_model.Fit.beta(2:2:2*best_model.Nosc)));
rateosc(2*best_model.Nosc+1:length(best_model.Fit.beta)) = 1./best_model.Fit.beta(2*best_model.Nosc+1:length(best_model.Fit.beta));

if isempty(StaticRate)
    rateosc(size(Af,1)-1) = 1;
    rateosc(size(Af,1)) = 1;
else
    rateosc(size(Af,1)) = 1./best_model.StatRates;
end

v = 30;
ind2d = dsearchn(T,100);
excolor = [.85 .85 .85];
vcol = [0.8500    0.3250    0.0980];
th = 1;

for ll = 1:size(rateoscmat,3)
rm = figure(1);
clf(1);

hold all
if isEven(ll) & ll < 2*best_model.Nosc+1;
    
contourf(re1,re3,abs(rateoscmat(:,:,ll)),v,'LineStyle','none');
crange = caxis;
contour(re1,re3,abs(rateoscmat(:,:,ll)),v/3,'color','k');

line(re1,line_diag_eq(re1,rateosc(ll)),'Color','k','linestyle','--')
line(re1,line_diag_eq(re1,-rateosc(ll)),'Color','k','linestyle','--')
line(re1,line_diag_eq(re1,-2*rateosc(ll)),'Color','k','linestyle','--')
line(re1,line_diag_eq(re1,+2*rateosc(ll)),'Color','k','linestyle','--')

else
    
contourf(re1,re3,(rateoscmat(:,:,ll)),v,'LineStyle','none');
crange = caxis;
contour(re1,re3,(rateoscmat(:,:,ll)),v/3,'color','k');

end

contour(re1,re3,reshape(X(:,ind2d),srx,sry,1),5,'color','w')
caxis(crange)

if isOdd(ll) & ll < 2*best_model.Nosc+1;
    title(['Oscillator \omega_2 = ' num2str(round(rateosc(ll+1))) ' cm^{-1} lifetime \tau = ' num2str(round(rateosc(ll))) ' fs'])
    figsavname = [num2str(ll) '-' num2str(round(rateosc(ll))) 'fs_osclifetime-w2' num2str(round(rateosc(ll+1))) 'wn'];
elseif isEven(ll) & ll < 2*best_model.Nosc+1;
    title(['Oscillator \omega_2 = ' num2str(round(rateosc(ll))) ' cm^{-1}'])
    figsavname = [num2str(ll) '-w2' num2str(rateosc(ll)) 'wn'];
elseif ll > 2*best_model.Nosc & ll < size(Af,1)
    title(['Dynamic rate lifetime \tau = ' num2str(round(rateosc(ll))) ' fs'])
    figsavname = [num2str(ll) '-' num2str(round(rateosc(ll))) 'fs_ratemap'];
else
    title(['Static rate lifetime \tau = ' num2str(round(rateosc(ll))) ' fs'])
    figsavname = ['static-' num2str(round(rateosc(ll))) 'fs_ratemap'];
end

line(re1,line_diag_eq(re1,0),'Color','k')

if length(excitonwn) == 1;
else
for g = 1:length(excitonwn);
line([excitonwn(g) excitonwn(g)],[min(re3) min(re3)+50],'linewidth',th+.3,'color',excolor);
line([excitonwn(g) excitonwn(g)],[max(re3)-50 max(re3)],'linewidth',th+.3,'color',excolor);
line([min(re1) min(re1)+50],[excitonwn(g) excitonwn(g)],'linewidth',th+.3,'color',excolor);
line([max(re1)-50 max(re1)],[excitonwn(g) excitonwn(g)],'linewidth',th+.3,'color',excolor);
end
end

colorbar
xlabel('\omega_1 (cm^{-1})')
ylabel('\omega_3 (cm^{-1})')
xlim([re1(end) re1(1)])
ylim([re3(end) re3(1)])
daspect([1 1 1])

ratemap = fullfile(ffolder,figsavname);
saveas(rm,[ratemap '.png'])
saveas(rm,[ratemap '.fig'])
end
toc
clear Af pf
%%
t2trpt = [excitonwn(1).*ones(size(excitonwn));...
    excitonwn(1).*ones(size(excitonwn));excitonwn;...
    fliplr(excitonwn)]'; 
offdiag = 0;

[b a] = max(max((real(X))));
[c b] = max((real(X(:,a))));

temptr = (Pf*af)';
srz = length(Trange);
if isempty(t2trpt)
    tr1 = squeeze(real(X(b,:)));
    tr2 = real(temptr(b,:));
    t2title = ['maximal signal point'];
    t2trpt = [b];
end

if size(t2trpt,1) == 1
    t2trpl = cat(1,t2trpt,t2trpt);
    offdiag = 0;
else
    offdiag = 1;
    t2trpl = t2trpt;
end

for m = 1:size(t2trpl,2)
    if ~isempty(t2trpl)
    wpt = t2trpl(:,m);
    xpt = dsearchn(re1',t2trpl(1,m))
    ypt = dsearchn(re3,t2trpl(2,m));
    tmtr = reshape(X,srx,sry,srz); 
    tr1 = squeeze(real(tmtr(ypt,xpt,:)));
    tmtr = reshape(temptr,srx,sry,srz);
    tr2 = squeeze(real(tmtr(ypt,xpt,:)));
    clear tmtr
    if ~offdiag
    trtitle = ['(\omega_1 = \omega_3) = ' num2str(t2trpl(1,m)) ' cm^{-1}'];
    else
    trtitle = ['\omega_1 = ' num2str(t2trpl(1,m)) ' cm^{-1}, \omega_3 = ' num2str(t2trpl(2,m)) ' cm^{-1}'];    
    end
    end

if ~isempty(t2trpl) && ~offdiag
    t2trn = [num2str(t2trpl(1,m)) 'wn-'];
elseif ~isempty(t2trpl) && offdiag
    t2trn = ['re1-' num2str(t2trpl(1,m)) '_re3-' num2str(t2trpl(2,m)) '-'];
else
    t2trn = '';
end
lwid = 1.5;

if ispc
t2trace = figure(3); 
clf(3)
ax1 = gca;
axdif = ax1.Position(4)*(1-1/1.6);
ax1.Position(4) = ax1.Position(4)/1.6;
plot(ax1,T,tr1,'linewidth',lwid)
hold all
plot(ax1,T,tr2,'linewidth',lwid);
xlabel('t_2 (fs)','Fontsize',10)
trm = max(abs(tr1));
% title(['t_2 trace of max r mat'])
leg1 = legend(['t_2 trace; S_{max} = ' num2str(trm)],'fit to trace'); set(leg1,'box','off','Location','best')
xlim([T2vec(1) T2vec(end)])
ax2 = axes('Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)+.1.*ax1.Position(2) ax1.Position(3) axdif]);
plot(ax2,T,tr1-tr2,'linewidth',lwid);
hold all
plot(get(gca,'XLim'),[0 0],'color','k')
xlim([T2vec(1) T2vec(end)])
ym = max(abs(tr1-tr2));
leg2 = legend(['Residual of fit; S_{max} = ' num2str(ym)]); set(leg2,'box','off','Location','best')
set(ax2,'XTickLabels','')
ylim([-ym*(1.1) ym*(1.1)])
title(['t_2 trace of rephasing ' trtitle])

set(gcf,'color','white')
t2tc = fullfile(ffolder,[t2trn 'real-t2-fit-trace']);
saveas(t2trace,[t2tc '.png'])
saveas(t2trace,[t2tc '.fig'])
% return
t2trace = figure(3); 
clf(3)
if ispc
ax1 = gca;
else
[ax1 ax2] = deal(get(axes));
end
axes(ax1)
axdif = ax1.Position(4)*(1-1/1.6);
ax1.Position(4) = ax1.Position(4)/1.6;
plot(ax1,T,tr1,'linewidth',lwid)
xlabel('T (fs)','Fontsize',10)
trm = max(abs(tr1));
% title(['t_2 trace of max r mat'])
leg1 = legend(['T trace']); set(leg1,'box','off','Location','best')
xlim([T2vec(1) T2vec(end)])
ax2 = axes('Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)+.1.*ax1.Position(2) ax1.Position(3) axdif]);
plot(ax2,T,tr1-tr2,'linewidth',lwid);
hold all
plot(get(gca,'XLim'),[0 0],'color','k')
xlim([T2vec(1) T2vec(end)])
ym = max(abs(tr1-tr2));
leg2 = legend(['Residual of fit']); set(leg2,'box','off','Location','best')
set(ax2,'XTickLabels','')
ylim([-ym*(1.1) ym*(1.1)])
% title(['t_2 trace of rephasing ' trtitle])
title(['T trace of rephasing ' trtitle])

set(gcf,'color','white')
t2tc = fullfile(ffolder,[t2trn 'real-t2-trace']);
saveas(t2trace,[t2tc '.png'])
saveas(t2trace,[t2tc '.fig'])

if m == 4
    [trcs ress] = deal(zeros(2,size(tr1,1)));
    trcs(1,:) = tr1;
    ress(1,:) = tr1-tr2;
elseif m == 5
    trcs(2,:) = tr1;
    ress(2,:) = tr1-tr2;
end

else
%%plot t2 trace on the data analysis computer (older matlab version)
close all
t2trace = figure(3); 
clf(3)
ax1 = gca;
td = get(ax1,'Position'); axdiff = td(4)*(1-1/1.6);
td(4) = td(4)/1.6;
set(ax1,'Position',td);
plot(ax1,T,tr1,'linewidth',lwid)
hold all
plot(ax1,T,tr2,'linewidth',lwid);
xlabel('t_2 (fs)','Fontsize',10)
trm = max(abs(tr1));
% title(['t_2 trace of max r mat'])
leg1 = legend(['t_2 trace; S_{max} = ' num2str(trm)],'fit to trace'); set(leg1,'box','off','Location','best')
xlim([T2vec(1) T2vec(end)])
td = get(ax1,'Position');
ax2 = axes;
set(ax2,'Position',[td(1) td(2)+td(4)+.1.*td(2) td(3) axdiff]);
plot(ax2,T,tr1-tr2,'linewidth',lwid);
hold all
plot(get(gca,'XLim'),[0 0],'color','k')
xlim([T2vec(1) T2vec(end)])
ym = max(abs(tr1-tr2));
leg2 = legend(['Residual of fit; S_{max} = ' num2str(ym)]); set(leg2,'box','off','Location','best')
set(ax2,'XTickLabel','')
ylim([-ym*(1.1) ym*(1.1)])
title(['t_2 trace of rephasing ' trtitle])

set(gcf,'color','white')
t2tc = fullfile(ffolder,[t2trn 'real-t2-fit-trace']);
saveas(t2trace,[t2tc '.png'])
saveas(t2trace,[t2tc '.fig'])
    
end
end
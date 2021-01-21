%%

excitonevset = [[1.881, 2.022, 1.971, 2.412, 1.798, 1.939, 1];...%WMS2-0 RT[1.9, 2.065, 2, 2.4, 2.152, 2.255];...%[
        [1.84 2.028 1.95 2.43 2.178 2.284 1];...% WMS2-0 80K
        [1.88 2.035 1.9518 2.31 1.796 1.953 1.97];...%WMS2-ua RT [1.8816 2.0097 1.9683 2.31 1.7987 1.939];...%[1.8731 2.0065 1.9497 2.31 1.7987 1.9447];...%[1.8817, 2.0277, 1.9541, 2.31, 1.7987, 1.9447];...
        [1.92 2.07 2.00 2.31 1.796 1.953 1.97];...%WMS2-ua 84K %1.91 2.05 2.00
        [1.884, 2.03, 1.97, 2.5, 2.5, 2.5, 1];...%MoS2 ML RT
        [1.92, 2.077, 1.865, 2.007, 2.152, 2.255, 1];...%MoS2 ML 80K
        [2.0048 2.39 2.85, 1.95, 2.152, 2.255, 2.0048];...%WS2 ML RT
        [2.05458 2.39, 1.95, 2.85, 2.152, 2.255, 2.0048];...%WS2 ML 84K
        [2.318 linspace(2.2545,2.4554,6)]];%[1.9917, 2.397, 2.8, 2.43, 2.152, 2.255]];
excitonlabset = {{'MoS2A','MoS2B','WS2A','WS2B','BGRA','BGRB','AltWS2'};...
        {'MoS2A','MoS2B','BGRA','BGRB'};...
        {'WS2A','WS2B','AltWS2','BGRA','BGRB'};...
        {'A','B','C','D','E','F'}};
    

if exist('proj','var')
    if contains(proj,'HS')
        if contains(subproj,'WMS2')
            if contains(twist,'ua')
                m = 1;
            else
                m = 0;
            end 

            if contains(temp,'RT')
                excitonev = excitonevset(1+2*m,1:3);
            else
                excitonev = excitonevset(2+2*m,1:3);
            end
            excitonlab = excitonlabset{1}(1:3);
        end
    elseif contains(proj,'ML')
        if contains(subproj,'MoS2')
            if contains(temp,'RT')
                excitonev = excitonevset(5,1:2);
            else
                excitonev = excitonevset(6,1:2);
            end
            excitonlab = excitonlabset{2};
        elseif contains(subproj,'WS2')
            if contains(temp,'RT')
                excitonev = excitonevset(7,1:3);
            else
                excitonev = excitonevset(8,1:3);
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
end
%%
 
if ispc
    basefolder = 'E:\Data';
    presfolder = 'C:\Users\User\Dropbox\Polimi\Presentations';
    sla = '\';
    xtrafolder = '';
else
    basefolder = '/Volumes/Utini/Data';
    presfolder = '/Users/User/Dropbox/Polimi/Presentations';
    sla = '/';
    xtrafolder = '/Users/veronicapolicht/MATLAB/';
end

%%

addpath('/usr/local/lib/matlab/')
addpath(fullfile(xtrafolder,'current-codes','2D Analysis','NLFittingNew'));
addpath(fullfile(xtrafolder,'current-codes','2D Analysis','DOPBoxV1_7','DOPBoxV1-7','DOPbox'));
addpath(fullfile(xtrafolder,'current-codes','2D Analysis','freezeColors_v23_cbfreeze','freezeColors'));
addpath(fullfile(xtrafolder,'current-codes','2D Analysis','plot2svg'));
addpath(fullfile(xtrafolder,'current-codes','2D Analysis','2DPhaseUnwrapping','2D_SRNCP_phase_unwrapper_without_a_mask'));

%%

c = 299.792458;
sol = 299.792458; % speed of light in nm/fs
T2wn = @(tau) 1./((c.*tau)./(1E7));
w2wn = @(w) 1E7.*w./(2*pi*c);
wn2f = @(wn) (c.*wn)./(1E7);
wn2w = @(wn) (2*pi*c.*wn)./(1E7);
wn2T = @(wn) 1./((c.*wn)./(1E7));
wn2eV = @(wn) 1.24e-4 * wn;
eV2wn = @(eV) eV/(1.24e-4);
T2wn = @(tau) 1./((c.*tau)./(1E7));
w2l = @(w) 2*pi*c./w;
l2w = @(l) (2*pi*c)./l;
l2wn = @(l) 1e7./l;
l2eV = @(l) 1240./l;
eV2l = @(eV) 1240./eV;
vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);
line_diag_eq = @(x,b) x+b;
rightlim = @(a) [min(a),max(a)];
mat2vec = @(sx,x,y) sx*(y-1) + x;
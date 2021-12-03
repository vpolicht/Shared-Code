%% This script is called at the beginning of nearly every other script:
% This serves to define the loation of your data and defines the slash for
% files (which is different on PCs vs Linux/Mac
%% Define your Data location Here:::

% basefolder - where to find the data
% presfolder - where to save figures if for a specific presentation
% sla - slash operator
% xtrafolder - where to find other scripts to reference

if ispc
    basefolder = 'E:\';
    presfolder = 'C:\Users\User\Dropbox\Polimi\Presentations';
    sla = '\';
    xtrafolder = 'C:\Users\User\Documents\GitHub';
else
    basefolder = '/Volumes/Utini';
    presfolder = '/Users/User/Dropbox/Polimi/Presentations';
    sla = '/';
    xtrafolder = '/Users/veronicapolicht/MATLAB/';
end

%% If using the 2D Coherence analysis, will need to add the following paths::

addpath(fullfile(xtrafolder,'2D Coherence Analysis','2D Analysis','NLFittingNew'));
addpath(fullfile(xtrafolder,'2D Coherence Analysis','2D Analysis','DOPBoxV1_7','DOPBoxV1-7','DOPbox'));
addpath(fullfile(xtrafolder,'2D Coherence Analysis','2D Analysis','freezeColors_v23_cbfreeze','freezeColors'));
addpath(fullfile(xtrafolder,'2D Coherence Analysis','2D Analysis','plot2svg'));
addpath(fullfile(xtrafolder,'2D Coherence Analysis','2D Analysis','2DPhaseUnwrapping','2D_SRNCP_phase_unwrapper_without_a_mask'));

%% Define some useful variables and functions

set(0,'defaulttextInterpreter','latex');

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
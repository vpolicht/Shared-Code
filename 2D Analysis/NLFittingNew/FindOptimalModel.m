close all; clear all; clc
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
c = 299.792458; %nm/fs
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
wn2T = @(wn) 1./((c*wn)./(1E7)); %wave number to period
T2wn = @(tau) 1./((c*tau)./(1E7));
w2l = @(w) 2*pi*c./w;
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFittingNew'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','DOPBoxV1_7','DOPBoxV1-7','DOPbox'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFitting','cm_and_cb_utilities','cm_and_cb_utilities'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','freezeColors_v23_cbfreeze','freezeColors'));
results_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','NewFits');
name_string = 'Optimal_Model_Rephasing_and_NonRephasing';

fit_result_cell = {'fit_results_2014-1-22-1-11',
                    'fit_results_2014-1-21-20-54'};
for k=1:length(fit_result_cell)
    load(fullfile(results_path,fit_result_cell{k}));
    aicvec = zeros(1,length(ModelCell));
    for j=1:length(ModelCell); aicvec(j) = ModelCell{j}.Stats.AICc; end
    [aicminval, aicminind] = min(aicvec);
    modelstruct(k).optmodel = ModelCell{aicminind};
    modelstruct(k).AICc = aicminval;
end
v = zeros(1,length(fit_result_cell));
for k=1:length(fit_result_cell); v(k) = modelstruct(k).AICc; end;
[~, best_model_ind] = min(v);
best_model = modelstruct(k).optmodel;
save(fullfile(results_path,name_string),'best_model');
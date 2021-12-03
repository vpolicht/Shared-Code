function [ Model ] = CreateSoLModelWithIRF( Nosc, NdynRates, StatRates, ConstBool, OscBounds, RateBounds, IRFBounds, t )
%This function creates Sum of Lorenzian models from inputs supplied.
% Nosc: The number of oscillators in the model
% NdynRates: The number of variable (dynamic) rates in the model (zero
%   frequency lorenzians)
% StatRates: A vector of static rates to include in the model.
% ConstBool: Include a constant offset?  0 if no.
% OscBounds: A structure containing 4 fields: fub, flb, rub, rlb.  These correspond
%   to the frequency (f) and rate (r) upper bound (ub) and lower bound (lb).
%   Each field is a vector, and length of all vectors must be the same.  
%   If you put 3 ub/lb pairs, this means that CreateSoLModel will create at
%   most 3 ub/lb unique pairs.  If Nosc < length(ub), then not all of the 
%   ub/lb pairs will be used.  If Nosc == length(ub), then one copy of all ub/lb pairs will be
%   used.  If Nosc > length(ub), then some copies of ub/lb pairs will be
%   repeated, it cycles through them cyclically.  Ex: if 3 pairs 1,2,3 are
%   given, and Nosc = 5, then you'll get 1,2,3,1,2.
% RateBounds: Similar to Osc Bounds, except that these control the upper
%   and lower bounds for the dynamic rates. There are only 2 fields: rub, rlb.
%   The controlling number is NdynRates, rather than Nosc.
% t: The vector of time points the model is to be evaluated at

% Returns:
% Model: A struct with 3 fields: 
%   modelfun: This is an anonymous function handle with a single argument: 
%   the non-linear parameter vector beta.
%   beta_lb: The lower bound of beta
%   beta_ub: the upper bound of beta

%%
assert(length(OscBounds.fub) == length(OscBounds.flb));
assert(length(OscBounds.rub) == length(OscBounds.rlb));
assert(length(OscBounds.fub) == length(OscBounds.rub));
assert(length(IRFBounds.lb)==2);
assert(length(IRFBounds.ub)==2);
Model.beta_ub = zeros(1,2*Nosc + NdynRates);
Model.beta_lb = zeros(1,2*Nosc + NdynRates);
osc_ind = repmat(1:length(OscBounds.rub),[1 ceil(Nosc/length(OscBounds.rub))]);
osc_ind = osc_ind(1:Nosc);
fubv = OscBounds.fub(osc_ind);
flbv = OscBounds.flb(osc_ind);
rubv = OscBounds.rub(osc_ind);
rlbv = OscBounds.rlb(osc_ind);
Model.beta_ub(1:2:(2*Nosc)) = rubv;
Model.beta_ub(2:2:(2*Nosc)) = fubv;
Model.beta_lb(1:2:(2*Nosc)) = rlbv;
Model.beta_lb(2:2:(2*Nosc)) = flbv;

assert(length(RateBounds.rub) == length(RateBounds.rlb));
rate_ind = repmat(1:length(RateBounds.rub),[1 ceil(NdynRates/length(RateBounds.rub))]);
rate_ind = rate_ind(1:NdynRates);
rubv = RateBounds.rub(rate_ind);
rlbv = RateBounds.rlb(rate_ind);
start_ind = (2*Nosc+1);
Model.beta_ub(start_ind:(start_ind+NdynRates-1)) = rubv;
Model.beta_lb(start_ind:(start_ind+NdynRates-1)) = rlbv;

Model.beta_ub = [IRFBounds.ub Model.beta_ub];
Model.beta_lb = [IRFBounds.lb Model.beta_lb];


Model.beta_ub = Model.beta_ub';
Model.beta_lb = Model.beta_lb';

Model.modelfun = @(x) multi_cexp_irf(x,Nosc,StatRates,ConstBool,t);
Model.Nosc = Nosc;
Model.NdynRates = NdynRates;
Model.ConstBool = ConstBool;
Model.StatRates = StatRates;
end


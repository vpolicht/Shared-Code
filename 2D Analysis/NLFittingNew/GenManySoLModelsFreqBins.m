function [ ModelCell, Diagnostics ] = GenManySoLModelsFreqBins(MaxNOsc, MinNOsc, MaxNdynRates, MinNdynRates, MaxConstRates, OscBounds, RateBounds, t)
%This function creates all possible model scenarios within the supplied bounds
%and stores them as a cell array of model structures.

%Inputs
%MaxNOsc: The maximum number of oscillators.  Will sweep models from 0
%   oscillators to the max.
%MaxNdynrates: The maximum number of dynamic rates.  Will sweep models from
%   0 dynamic rates to the max
%MaxConstRates: This is a vector of possible constant rates to choose from.
%   Will sweep from 0 of these rates to all of them.
%OscBounds: This is a structure holding the upper and lower bound
%   information for all oscillators in the model.  See CreateSoLModel for
%   description. If more than one range pair is supplied, this function
%   sweeps over all available range pairs.
%RateBounds: This is a structure holding the upper and lower bound
%   information for all dynamic rates in the model. See CreateSoLModel for
%   description.
%t: A vector of time points all models are evaluated at.

% At the minimum, there must be 1 dynamic rate.  Else the model will not be
% non-linear.

%This function can create A LOT of models in a hurry, so give it a whirl
%first before plugging it into your fitting software.
assert(length(OscBounds.fub) == length(OscBounds.flb));
assert(length(OscBounds.rub) == length(OscBounds.rlb));
assert(length(OscBounds.fub) == length(OscBounds.rub));
assert(length(RateBounds.rub) == length(RateBounds.rlb));
assert(length(OscBounds.fub) == 1);



LR = length(RateBounds.rub);
LS = length(MaxConstRates);

OscCombs = cell(1);
NoscCell = cell(1);
c = 0;
ls1 = @(s,e,n) s:((e-s)/(n)):(e-((e-s)/(n)));
ls2 = @(s,e,n) (s+((e-s)/(n))):((e-s)/(n)):e;
for k=MinNOsc:MaxNOsc
    c = c+1;
    NoscCell{c} = k;
end
if MaxNOsc<1
    NoscCell{1} = 0;
end
Diagnostics.OscCombs = OscCombs;

RateCombs = cell(1);
NrateCell = cell(1);
c = 0;
v = 1:LR;
for k=MinNdynRates:MaxNdynRates
    x = combsrep(v,k);
    if k==1; x = x'; end;
    for j=1:size(x,1);
        c = c+1;
        RateCombs{c} = x(j,:);
        NrateCell{c} = k;
    end
end
Diagnostics.RateCombs = RateCombs;

StatRateCombs = cell(1);
StatRateCombs{1} = [];
c = 1;
v = 1:LS;
for k=1:LS
    x = combsrep(v,k);
    if k==1; x = x'; end;
    for j=1:size(x,1);
        c = c+1;
        StatRateCombs{c} = x(j,:);
    end
end
Diagnostics.StatRateCombs = StatRateCombs;

ModelCell = cell(length(OscCombs)*length(RateCombs)*length(StatRateCombs)*2,1);
c = 0;
for i=1:length(NoscCell)
    if MaxNOsc>0;
        if NoscCell{i}==1;
            tOscBounds.rub = OscBounds.rub;
            tOscBounds.rlb = OscBounds.rlb;
            tOscBounds.fub = OscBounds.fub;
            tOscBounds.flb = OscBounds.flb;
        else
            tOscBounds.rub = OscBounds.rub(ones(NoscCell{i},1));
            tOscBounds.rlb = OscBounds.rlb(ones(NoscCell{i},1));
            tOscBounds.fub = 1./ls2(1/OscBounds.flb,1/OscBounds.fub,NoscCell{i});
            tOscBounds.flb = 1./ls1(1/OscBounds.flb,1/OscBounds.fub,NoscCell{i});
        end
        
    else
        tOscBounds.rub = 0;
        tOscBounds.rlb = 0;
        tOscBounds.fub = 0;
        tOscBounds.flb = 0;
    end
    for j=1:length(RateCombs)
        tRateBounds.rub = RateBounds.rub(RateCombs{j});
        tRateBounds.rlb = RateBounds.rlb(RateCombs{j});
        for k=1:length(StatRateCombs)
            for b=0:1
                c = c+1;
                ModelCell{c} = CreateSoLModel(NoscCell{i},NrateCell{j},MaxConstRates(StatRateCombs{k}),b,tOscBounds,tRateBounds,t);
            end
        end
    end
end



end


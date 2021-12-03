function [ ModelCell, Diagnostics ] = GenManySoLStretchedModels(MaxNOsc, MinNOsc, MaxNdynRates, MinNdynRates, MaxConstRates, OscBounds, RateBounds, t)
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
assert(length(OscBounds.sub) == length(OscBounds.slb));
assert(length(OscBounds.sub) == length(OscBounds.rub));
assert(length(RateBounds.rub) == length(RateBounds.rlb));
assert(length(RateBounds.sub) == length(RateBounds.slb));
if ~(length(RateBounds.sub) == length(RateBounds.rub))
    RateBounds.sub = repmat(RateBounds.sub,[1 length(RateBounds.rub)]);
    RateBounds.slb = repmat(RateBounds.slb,[1 length(RateBounds.rub)]);
end


LO = length(OscBounds.rub);
LR = length(RateBounds.rub);
LS = length(MaxConstRates);

OscCombs = cell(1);
NoscCell = cell(1);
c = 0;
v = 1:LO;
for k=MinNOsc:MaxNOsc
    x = combsrep(v,k);
    if k==1; x = x'; end;
    for j=1:size(x,1);
        c = c+1;
        OscCombs{c} = x(j,:);
        NoscCell{c} = k;
    end
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
for i=1:length(OscCombs)
    if MaxNOsc>0;
        tOscBounds.rub = OscBounds.rub(OscCombs{i});
        tOscBounds.rlb = OscBounds.rlb(OscCombs{i});
        tOscBounds.fub = OscBounds.fub(OscCombs{i});
        tOscBounds.flb = OscBounds.flb(OscCombs{i});
        tOscBounds.sub = OscBounds.sub(OscCombs{i});
        tOscBounds.slb = OscBounds.slb(OscCombs{i});
    else
        tOscBounds.rub = 0;
        tOscBounds.rlb = 0;
        tOscBounds.fub = 0;
        tOscBounds.flb = 0;
        tOscBounds.sub = 0;
        tOscBounds.slb = 0;
    end
    for j=1:length(RateCombs)
        tRateBounds.rub = RateBounds.rub(RateCombs{j});
        tRateBounds.rlb = RateBounds.rlb(RateCombs{j});
        tRateBounds.sub = RateBounds.sub(RateCombs{j});
        tRateBounds.slb = RateBounds.slb(RateCombs{j});
        for k=1:length(StatRateCombs)
            for b=0:1
                c = c+1;
                ModelCell{c} = CreateSoLStretchedModel(NoscCell{i},NrateCell{j},MaxConstRates(StatRateCombs{k}),b,tOscBounds,tRateBounds,t);
            end
        end
    end
end



end


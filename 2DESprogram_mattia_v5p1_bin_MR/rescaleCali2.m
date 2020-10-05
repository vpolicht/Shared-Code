
function [OptFreq,PseudoFreq ] = rescaleCali2( wlCal,MaxCorr,PP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
index=find(wlCal>=520);
index2=find(wlCal>=650);
OptFreq=300000./wlCal(index(1):index2(1));
PseudoFreq=MaxCorr*PP;
PseudoFreq=PseudoFreq(index(1):index2(1));
end

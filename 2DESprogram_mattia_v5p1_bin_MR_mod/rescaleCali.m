function [OptFreq,PseudoFreq ] = rescaleCali( wlCal,MaxCorr,PP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
index = find(wlCal >= 315);
index2 = find(wlCal >= 370);
OptFreq = 299792./wlCal(index(1):index2(1));
PseudoFreq = MaxCorr*PP;
PseudoFreq = PseudoFreq(index(1):index2(1));
end


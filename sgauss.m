function [sgau] = sgauss(x,sig,x0,pow)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

num = (x - x0).^2;
den = 2*(sig.^2);

sgau = exp(-(num./den).^pow);

end


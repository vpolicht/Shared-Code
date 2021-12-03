function [W] = MakeFourierOmegaAxis(T,L)
dT = mean(abs(diff(T)));
total_t = dT*L;
d_nu1 = 1/total_t;
W = 2*pi*((0:d_nu1:(L-1)*d_nu1));
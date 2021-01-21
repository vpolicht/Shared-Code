function [ Phi, dPhi, Ind] = multi_cexp_stretched( beta, Nosc, const_rates, const_bool, t )
%This function returns the real part of a complex exponential model
%convolved with a gaussian instrument response.  Functions here were
%computed symbolically in matlab, then traslated into anonymous functions
%via matlabFunction()
%
% INPUTS:
%
% beta: A vector of the non-linear variables.  These variables must be
% ordered in a specific way to accomodate all the options afforded by this
% function:
%   From elements 1 to 2: First is "time zero", the second is the FWHM of
%   the instrument response in units of time.
%
%   From elements 3 to 2*Nosc+2: Every odd element corresponds to a decay
%   rate, in units of 1/time.  Every even element corresponds to a frequency
%   in units of radians/time.
%   
%   From elements 2*Nosc+3 to the end: Each of these are rates
%   (zero-frequency complex exponentials).
%
% Nosc: The number of non-zero frequency oscillators to include in the
% model.
% 
% const_rates: Vectors describing exponential decays constant wrt non-linear
% parameters.  They are given as a vector of rates (1/time).
%
% const_bool: 1 if you'ld like to include a vector of 1s in the model (a
% constant offset) or 0 if you don't want to.
%
% t: a vector of time points we're evalating the model at.
%
% OUTPUTS:
% 
% Phi: The model, with all terms in it included as columns.  The number of
% columns will be length(beta)-2+size(const_rates,2)+const_bool
% dPhi: The first derivative of the model wrt to the non-linear variables.
%
% This is a matrix of columns, where only the non-zero derivatives of the
% model are returned.
%
% Ind: a matrix with 2 rows.  The first row tells you which column of Phi
% is having its derivative taken, and the second row tells you (in order of
% beta) what variable the derivative is being taken with respect to.

%% Function definitions
% Here are all the anonymous functions used to define the cexponentials and
% their derivatives.  Notation: 
%   pexp(t,k,m,d) = pure exponential decay (no oscillation) with rate k, time
%       zero m, instrument response fwhm d.
%   cexp(t,k,m,d,f) = exponential decay with cosine oscillation.  All
%       independent variables as previous, with f = frequency of oscillation
%   sexp(t,k,m,d,f) = exponential decay with sinusoidal oscillation.
%      Independent variables as for cexp.
%   d[f]d[var]: The general form for a derivative of some function [f] wrt
%   some [var].  For example dpexpdu = derivative of pure decay wrt the
%   time zero.


%***Comment on inclusion of both sine and cosine functions separately***
% The idea here is that we are fitting a real-value function, but allowing
% for an arbitrary phase shift of the oscillations.  Written
% trigonometrically, this means that the real part of a complex exponential
% with a real rate (the decay rate), an imaginary rate (frequency) and
% constant complex coefficient (the phase).  IE: 
% Re[exp[1i*phase]*exp[rate*t]*exp[1i*freq*t]] =
% Re[(a+1i*b)*exp[rate*t]*(cos[freq*t]+1i*sin[freq*t])] = 
% a*exp[rate*t]*cos[freq*t]-b*exp[rate*t]*sin[freq*t]

% Hope this all makes sense now.

% Note that in the conversion from symbollic to double precision math,
% matlab expands constants out to double precision.  So while it looks like
% I've hard coded numbers in here, it's really just things like
% log(sqrt(2)), 2*pi, and so forth.

%These functions can be generated with the script: DeriveStretchFunctional

cexpsin = @(f,r,s,t)exp(-r.*t.^s).*sin(f.*t);
cexpcos = @(f,r,s,t)exp(-r.*t.^s).*cos(f.*t);
cexp = @(r,s,t)exp(-r.*t.^s);
dcexpsinds = @(f,r,s,t)-r.*t.^s.*exp(-r.*t.^s).*sin(f.*t).*log(t);
dcexpsindr = @(f,r,s,t)-t.^s.*exp(-r.*t.^s).*sin(f.*t);
dcexpsindf = @(f,r,s,t)t.*exp(-r.*t.^s).*cos(f.*t);
dcexpcosds = @(f,r,s,t)-r.*t.^s.*exp(-r.*t.^s).*cos(f.*t).*log(t);
dcexpcosdr = @(f,r,s,t)-t.^s.*exp(-r.*t.^s).*cos(f.*t);
dcexpcosdf = @(f,r,s,t)-t.*exp(-r.*t.^s).*sin(f.*t);
dcexpds = @(r,s,t)-r.*t.^s.*exp(-r.*t.^s).*log(t);
dcexpdr = @(r,s,t)-t.^s.*exp(-r.*t.^s);

reg_exp = @(a,t) exp(-t*a); %a regular exponential.  Used for the fixed exponentials.




%% Matrix Filling

rate_params = beta(3*Nosc+1:end);

N = length(beta); %number of total non-linear variables
Nr = length(rate_params)/2; %number of pure stretched exponentials
Nc = length(const_rates); %number of constant rates
L = length(t); %number of time points
const_bool = const_bool>0;

%initialize with zeros
if const_bool;
    Phi = zeros(L,Nosc+Nr+Nc+1);
else
    Phi = zeros(L,Nosc+Nr+Nc);
end

%Fill the oscillators. Cosine in odd columns, sine in even columns.  Sorry,
%it's just tradition.  I know cosine is an even function.
c = 0;
for k=1:3:(3*Nosc)
    c = c+1;
    Phi(:,c) = cexpcos(beta(k),beta(k+1),beta(k+2),t);
    c = c+1;
    Phi(:,c) = cexpsin(beta(k),beta(k+1),beta(k+2),t);
end
for k=(Nosc*3+1):2:N
    c = c+1;
    Phi(:,c) = cexp(beta(k),beta(k+1),t);
end
for k=1:Nc;
    c = c+1;
    Phi(:,c) = reg_exp(const_rates(k),t);
end
if const_bool;
    c = c+1;
    Phi(:,c) = ones(L,1);
end

    if nargout>1
        NdPhi = 6*Nosc + 2*Nr; %2 functions per oscillation, and 3 derivatives per function + 2 derivatives per rate.

        dPhi = zeros(L,NdPhi);
        Ind = zeros(2,NdPhi);
        for k=1:Nosc
            Ind(1,6*(k-1)+1) = 2*(k-1)+1; %derivative of kth cosine oscillator
            Ind(2,6*(k-1)+1) = 3*(k-1)+1; %wrt beta(3*(k-1)+1), the frequency
            dPhi(:,6*(k-1)+1) = dcexpcosdf(beta(3*(k-1)+1),beta(3*(k-1)+2),beta(3*(k-1)+3),t);
            
            Ind(1,6*(k-1)+2) = 2*(k-1)+1; %derivative of kth cosine oscillator
            Ind(2,6*(k-1)+2) = 3*(k-1)+2; %wrt beta(3*(k-1)+2), the rate
            dPhi(:,6*(k-1)+2) = dcexpcosdr(beta(3*(k-1)+1),beta(3*(k-1)+2),beta(3*(k-1)+3),t);
            
            Ind(1,6*(k-1)+3) = 2*(k-1)+1; %derivative of kth cosine oscillator
            Ind(2,6*(k-1)+3) = 3*(k-1)+3; %wrt beta(3*(k-1)+3), the stretch
            dPhi(:,6*(k-1)+3) = dcexpcosds(beta(3*(k-1)+1),beta(3*(k-1)+2),beta(3*(k-1)+3),t);
            
            Ind(1,6*(k-1)+4) = 2*(k-1)+2; %derivative of kth sine oscillator
            Ind(2,6*(k-1)+4) = 3*(k-1)+1; %wrt beta(3*(k-1)+1), the frequency
            dPhi(:,6*(k-1)+4) = dcexpsindf(beta(3*(k-1)+1),beta(3*(k-1)+2),beta(3*(k-1)+3),t);
            
            Ind(1,6*(k-1)+5) = 2*(k-1)+2; %derivative of kth sine oscillator
            Ind(2,6*(k-1)+5) = 3*(k-1)+2; %wrt beta(3*(k-1)+2), the rate
            dPhi(:,6*(k-1)+5) = dcexpsindr(beta(3*(k-1)+1),beta(3*(k-1)+2),beta(3*(k-1)+3),t);
            
            Ind(1,6*(k-1)+6) = 2*(k-1)+2; %derivative of kth sine oscillator
            Ind(2,6*(k-1)+6) = 3*(k-1)+3; %wrt beta(3*(k-1)+3), the stretch
            dPhi(:,6*(k-1)+6) = dcexpsinds(beta(3*(k-1)+1),beta(3*(k-1)+2),beta(3*(k-1)+3),t);
        end
        if NdPhi>6*Nosc;
            for k=1:Nr;               
                Ind(1,6*Nosc+2*(k-1)+1) = 2*Nosc+k; %derv of 2*Nosc + kth function (a pure exponential decay)
                Ind(2,6*Nosc+2*(k-1)+1) = 3*Nosc+2*(k-1)+1; %wrt beta(3*Nosc+k), its rate
                dPhi(:,6*Nosc+2*(k-1)+1) = dcexpdr(beta(3*Nosc+k),beta(3*Nosc+k+1),t);
                
                Ind(1,6*Nosc+2*(k-1)+2) = 2*Nosc+k; %derv of 2*Nosc + kth function (a pure exponential decay)
                Ind(2,6*Nosc+2*(k-1)+2) = 3*Nosc+2*(k-1)+2; %wrt beta(2*Nosc+k), its rate
                dPhi(:,6*Nosc+2*(k-1)+2) = dcexpds(beta(3*Nosc+k),beta(3*Nosc+k+1),t);
            end
        end
    end
end


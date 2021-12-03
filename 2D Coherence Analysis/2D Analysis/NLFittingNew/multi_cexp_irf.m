function [ Phi, dPhi, Ind] = multi_cexp_irf( beta, Nosc, const_rates, const_bool, t )
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

pexp = @(t,k,m,d)exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.0./2.0);
cexp = @(t,k,m,d,f)exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.0./2.0);
sexp = @(t,k,m,d,f)exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.0./2.0);

dsexpdk = @(t,k,m,d,f)t.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(1.0./2.0)-exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(m+d.^2.*k.*3.606737602222409e-1).*(1.0./2.0)-sqrt(2.0).*1.0./sqrt(pi).*d.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*3.002806021966125e-1;
dcexpdk = @(t,k,m,d,f)t.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(1.0./2.0)-exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(m+d.^2.*k.*3.606737602222409e-1).*(1.0./2.0)-sqrt(2.0).*1.0./sqrt(pi).*d.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*3.002806021966125e-1;
dpexpdk = @(t,k,m,d)t.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(1.0./2.0)-exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(m+d.^2.*k.*3.606737602222409e-1).*(1.0./2.0)-sqrt(2.0).*1.0./sqrt(pi).*d.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*3.002806021966125e-1;

dsexpdm = @(t,k,m,d,f)k.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.0./2.0)-(sqrt(2.0).*1.0./sqrt(pi).*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*8.325546111576977e-1)./d;
dcexpdm = @(t,k,m,d,f)k.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.0./2.0)-(sqrt(2.0).*1.0./sqrt(pi).*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*8.325546111576977e-1)./d;
dpexpdm = @(t,k,m,d)k.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.0./2.0)-(sqrt(2.0).*1.0./sqrt(pi).*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*8.325546111576977e-1)./d;

dsexpdd = @(t,k,m,d,f)d.*k.^2.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.803368801111205e-1)-1.0./sqrt(pi).*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*(sqrt(2.0).*k.*6.005612043932249e-1-sqrt(2.0).*1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1);
dcexpdd = @(t,k,m,d,f)-1.0./sqrt(pi).*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*(sqrt(2.0).*k.*6.005612043932249e-1-sqrt(2.0).*1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)-d.*k.^2.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*1.803368801111205e-1;
dpexpdd = @(t,k,m,d)d.*k.^2.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0).*(-1.803368801111205e-1)-1.0./sqrt(pi).*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).^2.*(-1.38629436111989)).*exp(-k.*t).*(sqrt(2.0).*k.*6.005612043932249e-1-sqrt(2.0).*1.0./d.^2.*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1);

dsexpdf = @(t,k,m,d,f)-pi.*t.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*cos(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0);
dcexpdf = @(t,k,m,d,f)pi.*t.*exp(k.*(m+d.^2.*k.*1.803368801111205e-1)).*exp(-k.*t).*sin(pi.*f.*t.*2.0).*(erf((sqrt(2.0).*(m-t+d.^2.*k.*3.606737602222409e-1).*8.325546111576977e-1)./d)-1.0);

reg_exp = @(a,t) exp(-t*a); %a regular exponential.  Used for the fixed exponentials.




%% Matrix Filling

rate_params = beta(2*Nosc+3:end);

N = length(beta); %number of total non-linear variables
Nr = length(rate_params); %number of pure decays
Nc = length(const_rates); %number of constant rates
L = length(t); %number of time points
const_bool = const_bool>0;

%initialize with zeros
if const_bool;
    Phi = zeros(L,N-2+Nc+1);
else
    Phi = zeros(L,N-2+Nc);
end

%Fill the oscillators. Cosine in odd columns, sine in even columns.  Sorry,
%it's just tradition.  I know cosine is an even function.
for k=1:2:(2*Nosc)
    Phi(:,k) = cexp(t,beta(k+2),beta(1),beta(2),beta(k+3));
    Phi(:,k+1) = sexp(t,beta(k+2),beta(1),beta(2),beta(k+3));
end
for k=(Nosc*2+3):N
    Phi(:,k-2) = pexp(t,beta(k),beta(1),beta(2));
end
c = 0;
for k=(N-2+1):(N-2+Nc);
    c = c+1;
    Phi(:,k) = reg_exp(t,const_rates(c));
end
if const_bool;
    Phi(:,N-2+Nc+1) = ones(L,1);
end

    if nargout>1
        NdPhi = 8*Nosc + 3*Nr; %2 functions per oscillation, and 4 derivatives per function + 3 derivatives per rate.

        dPhi = zeros(L,NdPhi);
        Ind = zeros(2,NdPhi);
        for k=1:2:(2*Nosc)
            Ind(1,4*(k-1)+1) = k; %derivative of kth cosine oscillator
            Ind(2,4*(k-1)+1) = 1; %wrt beta(1), ie the time zero.
            dPhi(:,4*(k-1)+1) = dcexpdm(t,beta(k+2),beta(1),beta(2),beta(k+3));
            
            Ind(1,4*(k-1)+2) = k; %derivative of kth cosine oscillator
            Ind(2,4*(k-1)+2) = 2; %wrt beta(2), ie the response fwhm
            dPhi(:,4*(k-1)+2) = dcexpdd(t,beta(k+2),beta(1),beta(2),beta(k+3));
            
            Ind(1,4*(k-1)+3) = k; %derivative of kth cosine oscillator
            Ind(2,4*(k-1)+3) = k+2; %wrt beta(k+2), ie its rate
            dPhi(:,4*(k-1)+3) = dcexpdk(t,beta(k+2),beta(1),beta(2),beta(k+3));
            
            Ind(1,4*(k-1)+4) = k; %derivative of kth cosine oscillator
            Ind(2,4*(k-1)+4) = k+3; %wrt beta(k+3), ie its frequency
            dPhi(:,4*(k-1)+4) = dcexpdf(t,beta(k+2),beta(1),beta(2),beta(k+3));
            
            Ind(1,4*(k-1)+5) = k+1; %derivative of kth sine oscillator
            Ind(2,4*(k-1)+5) = 1; %wrt beta(1), ie the time zero.
            dPhi(:,4*(k-1)+5) = dsexpdm(t,beta(k+2),beta(1),beta(2),beta(k+3));
            
            Ind(1,4*(k-1)+6) = k+1; %derivative of kth sine oscillator
            Ind(2,4*(k-1)+6) = 2; %wrt beta(2), ie the response fwhm
            dPhi(:,4*(k-1)+6) = dsexpdd(t,beta(k+2),beta(1),beta(2),beta(k+3));
            
            Ind(1,4*(k-1)+7) = k+1; %derivative of kth sine oscillator
            Ind(2,4*(k-1)+7) = k+2; %wrt beta(k+2), ie its rate
            dPhi(:,4*(k-1)+7) = dsexpdk(t,beta(k+2),beta(1),beta(2),beta(k+3));
            
            Ind(1,4*(k-1)+8) = k+1; %derivative of kth sine oscillator
            Ind(2,4*(k-1)+8) = k+3; %wrt beta(k+3), ie its frequency
            dPhi(:,4*(k-1)+8) = dsexpdf(t,beta(k+2),beta(1),beta(2),beta(k+3));
        end
        if NdPhi>8*Nosc;
            for k=1:Nr;               
                Ind(1,8*Nosc+3*(k-1)+1) = 2*Nosc+k; %derv of 2*Nosc + kth function (a pure exponential decay)
                Ind(2,8*Nosc+3*(k-1)+1) = 1; %wrt beta(1), ie its time zero
                dPhi(:,8*Nosc+3*(k-1)+1) = dpexpdm(t,beta(2*Nosc+k+2),beta(1),beta(2));
                
                Ind(1,8*Nosc+3*(k-1)+2) = 2*Nosc+k; %derv of 2*Nosc + kth function (a pure exponential decay)
                Ind(2,8*Nosc+3*(k-1)+2) = 2; %wrt beta(2), ie the response fwhm 
                dPhi(:,8*Nosc+3*(k-1)+2) = dpexpdd(t,beta(2*Nosc+k+2),beta(1),beta(2));
                
                Ind(1,8*Nosc+3*(k-1)+3) = 2*Nosc+k; %derv of 2*Nosc + kth function (a pure exponential decay)
                Ind(2,8*Nosc+3*(k-1)+3) = 2*Nosc+k+2; %wrt 2*Nosc + kth + 2 value of beta (its rate)
                dPhi(:,8*Nosc+3*(k-1)+3) = dpexpdk(t,beta(2*Nosc+k+2),beta(1),beta(2));
            end
        end
    end
end


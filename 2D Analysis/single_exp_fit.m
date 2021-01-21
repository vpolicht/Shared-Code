function [residual,fit_to_data, parameters] = single_exp_fit(T, x)
%CREATEFIT(T,Z6_15)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : T
%      Y Output: x
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 13-Jun-2013 17:16:41


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( T, x );

% Set up fittype and options.
ft = @(x,xdata) (x(1)*exp(-x(2)*xdata)+x(3));
x0(1) = yData(1);
x0(2) = 0;
x0(3) = yData(end);
% Fit model to data.
% [fitresult, ~] = fit( xData, yData, ft, opts );
options = optimset('TolFun',1E-17,'TolX',1E-12,'Display','off');
p = lsqcurvefit(ft,x0,xData,yData,[],[],options);


residual = yData-p(1)*exp(-p(2)*xData)+p(3);
fit_to_data = p(1)*exp(-p(2)*xData)+p(3);
parameters = p;

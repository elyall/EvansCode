function [fitfull, fullgof] = FitFull(fulldriven)
%CREATEFIT(FULLDRIVEN)
%  Create a fit.
%
%  Data for 'fullwhisker' fit:
%      Y Output: fulldriven
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 12-Oct-2012 13:44:51


%% Fit: 'fullwhisker'.
[xData, yData] = prepareCurveData( [], fulldriven );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
%opts.StartPoint = [19.403958660351 0.204124145231932 0.278969503376246 14.9891752158003 -0.204124145231932 0.510482646646947];
opts.Upper = [100 100 100];
opts.Normalize = 'off';

% Fit model to data.
[fitfull, fullgof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'fullwhisker' );
% plot( fitfull, xData, yData );
% % Label axes
% ylabel( 'fulldriven' );
% grid off



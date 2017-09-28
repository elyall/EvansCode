function p = exp_nonlin(b,x)
% b (1:end-1) = linear filter coefficients
% b (end) = threshold coeff
% x = design matrix

y = x*b(1:end-1);
y = y + b(end); % threshold shift
% p = exp(y);     % exponential
p = y.^2;       % quadratic
% y(y<0) = 0;     % rectified
% p = y;          % linear


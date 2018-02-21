function p = div_stim_norm(b,x)
% b (1:end-1) = linear filter coefficients
% b (end) = normalization scaling factor
% x = design matrix

y = x*b(1:end-1);

p = y./(1 + b(end) * sum(x(:,2:6),2).^2); 

% p = exp(y);     % exponential
% p = y.^2;       % quadratic
% y(y<0) = 0;     % rectified
% p = y;          % linear

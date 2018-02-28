function p = div_pop_norm(b,x)
% b (1:end-1) = linear filter coefficients
% b (end) = normalization scaling factor
% x = design matrix

y = x(:,1:end-1)*b(1:end-1);

p = y./(1 + b(end) * x(:,end).^2); 

% p = exp(y);     % exponential
% p = y.^2;       % quadratic
% y(y<0) = 0;     % rectified
% p = y;          % linear

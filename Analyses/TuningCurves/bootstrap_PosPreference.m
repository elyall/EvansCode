function [s_diff, p] = bootstrap_PosPreference(Pref, Sel, num_samples)

% Compute variance
Prb = bsxfun(@rdivide, Sel, sum(Sel));      % calculate weight
Exp = sum(Prb.*Pref);                       % compute expected value (weighted mean)
Var = sum(Prb.*(Pref.^2)) - Exp.^2;         % calculate variance

% Compute change in variance
s_diff = Var(1)-Var(2);                     

% Bootstrap
numROIs = size(Pref,1);
null_dist = zeros(num_samples, 1);
for k = 1:num_samples
    
    % Select data
    [~, idx] = datasample(Pref(:), 2*numROIs, 'Replace', true);
    pref = Pref(reshape(idx,numROIs,2));
    sel = Sel(reshape(idx,numROIs,2));
    
    % Compute variance
    Prb = bsxfun(@rdivide, sel, sum(sel)); % calculate probability
    Exp = sum(Prb.*pref); % compute
    Var = sum(Prb.*(pref.^2)) - Exp.^2;
    
    % Save change to distribution
    null_dist(k) = Var(1) - Var(2);
    
end

% Calculate p-value off of bootstrap
p = (sum(null_dist > abs(s_diff)) + sum(null_dist < -abs(s_diff)))/num_samples;


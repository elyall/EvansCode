function Sel = computeSelectivity(Curves)

% % flip inhibited curves
index = max(Curves,[],2)<max(abs(Curves),[],2);
Curves(index,:) = -Curves(index,:);

% % shift bottom of curves to 0
% Curves = bsxfun(@minus, Curves, min(Curves, [], 2));

% normalize curves
Curves = bsxfun(@rdivide, Curves, max(Curves,[],2));

% compute integral
Sel = trapz(Curves,2);
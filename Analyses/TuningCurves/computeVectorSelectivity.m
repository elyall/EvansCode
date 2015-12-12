function Sel = computeVectorSelectivity(Curves)

% shift curves to 0
Curves = bsxfun(@minus, Curves, min(Curves,[],2));

% compute vector selectivity
Sel = 1 - (sqrt(sum(abs(Curves).^2,2))./max(Curves,[],2)-1)/(sqrt(size(Curves,2))-1);
function Corr = computeCurveCorrelation(Curves)

% % shift bottom of curves to 0
% Curves = bsxfun(@minus, Curves, min(Curves, [], 2));

% normalize curves
Curves = bsxfun(@rdivide, Curves, max(Curves,[],2));

% subtract off mean
Curves = bsxfun(@minus, Curves, mean(Curves,1));

% compute dot product
numROIs = size(Curves, 1);
Corr = nan(numROIs);
for rindex = 1:numROIs
    Corr(rindex,rindex+1:end) = sum(bsxfun(@times, Curves(rindex,:), Curves(rindex+1:end,:)),2);
%     Corr(rindex,rindex+1:end) = sum(abs(bsxfun(@minus, Curves(rindex,:), Curves(rindex+1:end,:))),2);
end
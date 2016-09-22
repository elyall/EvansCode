function [BW,vertices,ellipse] = segmentImage(X)

invert = true;
iterations = 300;

verbose = true;
pauseDur = .2;


% Threshold image - global threshold
BW = imbinarize(X);
if verbose; figure; imagesc(BW); pause(pauseDur); end

% Invert mask
if invert
    BW = imcomplement(BW);
    if verbose; imagesc(BW); pause(pauseDur); end
end

% % Clear border
% BW = imclearborder(BW);

% Fill holes
BW = imfill(BW, 'holes');
if verbose; imagesc(BW); pause(pauseDur); end

% Expand regions if too little is labeled
while nnz(BW)<.1*numel(BW)
    radius = 5;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    BW = imdilate(BW, se);
    if verbose; imagesc(BW); pause(pauseDur); end
end

% Open mask with disk
radius = 50;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);
if verbose; imagesc(BW); pause(pauseDur); end

% Erode mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);
if verbose; imagesc(BW); pause(pauseDur); end

% Active contour using Chan-Vese
BW = activecontour(X, BW, iterations, 'Chan-Vese');
if verbose; imagesc(BW); pause(pauseDur); end

% Erode mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);
if verbose; imagesc(BW); pause(pauseDur); end

% Open mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);
if verbose; imagesc(BW); pause(pauseDur); end

% Remove all but the largest blob
labeledImage = bwlabel(BW);
blobMeasurements = regionprops(labeledImage, 'area', 'Centroid');
[~,ind] = max([blobMeasurements.Area]);
BW(labeledImage~=ind) = false;
if verbose; imagesc(BW); pause(pauseDur); end

% % Fill any holes that might exist
% BW = imfill(BW, 'holes');
% if verbose; imagesc(BW); pause(pauseDur); end

% Find vertices of blob
vertices = bwboundaries(BW);
vertices = vertices{1};
vertices = fliplr(vertices);
if verbose
    imagesc(X); hold on;
    plot(vertices([1:end,1],1),vertices([1:end,1],2),'g--');
end

% Fit ellipse
if ~verbose
    ellipse = fit_ellipse(vertices(:,1),vertices(:,2));
else
    ellipse = fit_ellipse(vertices(:,1),vertices(:,2),gca);
end


function [BW,vertices,ellipse] = segmentImage(X)

verbose = true;

% Threshold image - global threshold
BW = imbinarize(X);
if verbose; figure; imagesc(BW); drawnow; end

% Invert mask
BW = imcomplement(BW);
if verbose; imagesc(BW); drawnow; end

% Fill holes
BW = imfill(BW, 'holes');
if verbose; imagesc(BW); drawnow; end

% Open mask with disk
radius = 50;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);
if verbose; imagesc(BW); drawnow; end

% Erode mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);
if verbose; imagesc(BW); drawnow; end

% Active contour using Chan-Vese over 495 iterations
iterations = 300;
BW = activecontour(X, BW, iterations, 'Chan-Vese');
if verbose; imagesc(BW); drawnow; end

% Erode mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);
if verbose; imagesc(BW); drawnow; end

% Open mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);
if verbose; imagesc(BW); drawnow; end

% Remove all but the largest blob
labeledImage = bwlabel(BW);
blobMeasurements = regionprops(labeledImage, 'area', 'Centroid');
[~,ind] = max([blobMeasurements.Area]);
BW(labeledImage~=ind) = false;
if verbose; imagesc(BW); drawnow; end

% Fill holes
BW = imfill(BW, 'holes');
if verbose; imagesc(BW); drawnow; end

% Find vertices
vertices = bwboundaries(BW);
vertices = vertices{1};

if verbose
    
    % Display output
    imagesc(X); hold on;
    plot(vertices([1:end,1],2),vertices([1:end,1],1),'g--');
    
    % Fit ellipse
    ellipse = fit_ellipse(vertices(:,2),vertices(:,1),gca);
    
else
    
    % Fit ellipse
    ellipse = fit_ellipse(vertices(:,2),vertices(:,1));
    
end
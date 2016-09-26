function [BW,vertices,ellipse] = segmentImage(X,varargin)

invert = true;
iterations = 300;

verbose = false;
pauseDur = .2;
hA = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'invert'
                invert = ~invert;
                index = index + 1;
            case 'iterations'
                iterations = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = ~verbose;
                index = index + 1;
            case 'hA'
                hA = varargin{index+1};
                index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if verbose && isempty(hA)
    figure;
    hA = axes();
end


%% Segment image
% Threshold image - global threshold
BW = imbinarize(X);
plotAxis(BW);

% Invert mask
if invert
    BW = imcomplement(BW);
    plotAxis(BW);
end

% % Clear border
% BW = imclearborder(BW);

% Fill holes
BW = imfill(BW, 'holes');
plotAxis(BW);

% Expand regions if too little is labeled
while nnz(BW)<.1*numel(BW)
    radius = 5;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    BW = imdilate(BW, se);
    plotAxis(BW);
end

% Open mask with disk
radius = 50;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);
plotAxis(BW);

% Erode mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);
plotAxis(BW);

% Active contour using Chan-Vese
BW = activecontour(X, BW, iterations, 'Chan-Vese');
plotAxis(BW);

% Erode mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);
plotAxis(BW);

% Open mask with disk
radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);
plotAxis(BW);

% Fill holes
BW = imfill(BW, 'holes');
plotAxis(BW);

% Remove all but the largest blob
labeledImage = bwlabel(BW);
blobMeasurements = regionprops(labeledImage, 'area', 'Centroid');
[~,ind] = max([blobMeasurements.Area]);
BW(labeledImage~=ind) = false;
plotAxis(BW);

% % Fill any holes that might exist
% BW = imfill(BW, 'holes');
% plotAxis(BW);

% Find vertices & smooth vertices of blob
vertices = bwboundaries(BW);
vertices = vertices{1};
vertices = fliplr(vertices);
vertices(:,1) = smooth(vertices(:,1),20,'sgolay');
vertices(:,2) = smooth(vertices(:,2),20,'sgolay');
if verbose
    axes(hA);
    imagesc(X); hold on;
    plot(vertices([1:end,1],1),vertices([1:end,1],2),'g--');
end

% Fit ellipse
if ~verbose
    ellipse = fit_ellipse(vertices(:,1),vertices(:,2));
else
    ellipse = fit_ellipse(vertices(:,1),vertices(:,2),hA);
end


%% Plotting Subfunction
    function plotAxis(img)
        if verbose
            axes(hA);
            imagesc(img);
            pause(pauseDur);
        end
    end
end
function [ROIs,Centroids,Labels] = UIroi(N,Image,Labels)

directory = cd;

%% Parse input arguments
if ~exist('N','var') || isempty(N)
    N = 1;
end

if ~exist('Image','var') || isempty(Image)
    [Image,p] = uigetfile({'*.tif'},'Choose iOS image',directory);
    if isnumeric(Image)
        return
    end
    Image = fullfile(p,Image);
end

if ~exist('Labels','var') || isempty(Labels)
    Labels = num2cell(num2str((1:N)'));
end

%% Load image
if ischar(Image)
    Image = imread(Image);
end

%% Define ROIs
ROIs = cell(1,N);
Centroids = nan(N,2);

hF = figure('NumberTitle','off','CloseRequestFcn',@DeleteFcn);
imagesc(Image); colormap(gray);
fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
for rindex = 1:N
    go = false;
    set(hF,'Name',sprintf('Select ROI %s (exit figure when finished)',Labels{rindex}));
    h = impoly(gca,'Closed',1,'PositionConstraintFcn',fcn);
    while ~go
        pause(.2)
    end
end


    function DeleteFcn(hObject,eventdata)
        ROIs{rindex} = getPosition(h);
        temp = regionprops(createMask(h), 'Centroid');
        Centroids(rindex,:) = temp.Centroid;
        delete(h);
        go = true;
        if rindex==N
            delete(hF);
        end
    end
end
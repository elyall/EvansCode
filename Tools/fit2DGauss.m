function Fit = fit2DGauss(Image,ROI,init)

% Filter = false;
Filter = fspecial('gaussian',5,1);
verbose = false;

%% Parse input arguments
if ~exist('Image','var') || isempty(Image)
    [Image,p] = uigetfile({'*.tif'},'Choose iOS image');
    if isnumeric(Image)
        return
    end
    Image = fullfile(p,Image);
end
if ischar(Image)
    Image = imread(Image);
end
if ~isa(Image,'double')
    Image = double(Image);
end

if ~exist('ROI','var')
    ROI = [];
end

if ~exist('init','var')
    init = [];
end

%% UI

% Define ROI
if isequal(ROI,true)
    hF = figure('NumberTitle','off','Name','Select ROI');
    imagesc(Image); colormap(gray);
    fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
    h = impoly(gca,'Closed',1,'PositionConstraintFcn',fcn);
    ROI = getPosition(h);
    close(hF);
end

% Make guess at centroid
if isequal(init,true)
    hF = figure('NumberTitle','off','Name','Select guess of centroid');
    imagesc(Image); colormap(gray);
    init = round(flip(ginput(1)));
    close(hF);
end


%% Fit 2D Gaussian
temp = Image;
if ~isequal(Filter,false)
    temp = imfilter(temp,Filter);
end
if ~isempty(ROI)
    mask = poly2mask(ROI(:,1),ROI(:,2),size(Image,1),size(Image,2));
    temp(~mask) = mean(temp(:));
end
[xx,yy] = meshgrid(1:size(Image,1),1:size(Image,2));
[Fit, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(xx,yy,temp,init);


%% Display output
if verbose
    figure;
    imagesc(Image); hold on;
    plotGauss(Fit(5:6), Fit(2), Fit(3:4));
    if ~isempty(ROI)
        plot(ROI([1:end,1],1),ROI([1:end,1],2),'r--','LineWidth',2);
    end
end

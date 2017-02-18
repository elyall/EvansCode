function Fit = fit2DGauss(Image,ROI,init)


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

if ~exist('ROI','var')
    ROI = [];
end

if ~exist('init','var')
    init = [];
end

%% Define ROI
if isequal(ROI,true)
    hF = figure;
    imagesc(Image); colormap(gray);
    fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
    h = impoly(gca,'Closed',1,'PositionConstraintFcn',fcn);
    ROI = getPosition(h);
    close(hF);
end

if ~isempty(ROI)
    mask = poly2mask(ROI(:,1),ROI(:,2),size(Image,1),size(Image,2));
    Image(~mask) = intmax(class(Image));
end

%% Make guess
if isempty(init)
    
end

%% Fit 2D Gaussian
Image = double(intmax(class(Image)) - Image);
Fit = D2GaussFit(Image,guess);

%% Display output
figure;
imagesc(image); hold on;
% overlay gaussian
plot(ROI([1:end,1],1),ROI([1:end,1],2),'r-','LineWidth',2);

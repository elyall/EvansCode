function [Projection,theta,rsquared] = projectOntoCoMAxis(Centroids, Data, varargin)

insidePWC = [];

% Display and video
ImageSize = [512,796];
verbose = false;
saveVideo = false;
VideoFile = '';

% % Saving output
% saveOut = false;
% saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'insidePWC'
                insidePWC = varargin{index+1};
                index = index + 2;
            case 'ImageSize'
                ImageSize = varargin{index+1};
                index = index + 2;
            case 'saveVideo'
                saveVideo = true;
                index = index + 1;
            case 'VideoFile'
                VideoFile = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
%             case {'Save', 'save'}
%                 saveOut = true;
%                 index = index + 1;
%             case {'SaveFile', 'saveFile'}
%                 saveFile = varargin{index+1};
%                 index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Centroids','var') || isempty(Centroids)
    [Centroids, p] = uigetfile({'*.rois'},'Choose ROI file(s)', directory,'MultiSelect','on');
    if isnumeric(Centroids)
        return
    end
    Centroids = fullfile(p,Centroids);
end
if ischar(Centroids)
    Centroids = {Centroids};
end

if ~exist('Data','var')
    Data = [];
end

% if saveOut && isempty(saveFile)
%     saveOut = false;
% end
% if saveVideo && isempty(VideoFile)
%     saveVideo = false;
% end


%% Load in data

% Load in ROIdata
if iscellstr(Centroids)
    ROIFiles = Centroids;
    for findex = 1:numFiles
        temp = load(ROIFiles{findex}, 'ROIdata', '-mat');
        Centroids{findex} = temp.ROIdata;
        clear temp;
    end
end

if iscell(Centroids) || isstruct(Centroids)
    
    % Load in data
    if isempty(Data)
        Curves = gatherROIdata(Centroids, 'curve', ':', 'none');
        Data = computeCenterOfMass(Curves,2);
    end
    
    % Load in centroids
    Centroids = gatherROIdata(Centroids, 'centroid', ':', 'none');

end

if size(Centroids,1)~=numel(Data)
    error('Number of centroids input does not equal number of data points input!');
end

    
%% Whiten centroids
Mean = mean(Centroids);                             % calculate mean
ZeroedCentroids = bsxfun(@minus, Centroids, Mean);  % subtract off mean (center data)
A = ZeroedCentroids'*ZeroedCentroids;
[V,D,~] = svd(A);                                   % compute SVD
if any(diag(D)==0)                                  % ensure all points on diagonal have some weight
    ind = find(diag(D)==0);
    D(ind,ind) = 0.00001;
end
whMat = sqrt(size(ZeroedCentroids,1)-1)*V*sqrtm(inv(D))*V'; % determine whitening transform
WhiteCentroids = ZeroedCentroids*whMat;                     % whiten data
invMat = pinv(whMat);                                       % determine inverse of transform

% if verbose % Check whitening maintains underlying structure
%     figure;
%     subplot(1,2,1);
%     plot(WhiteCentroids(:,1),WhiteCentroids(:,2),'.');
%     subplot(1,2,2);
%     plot(Centroids(:,1),Centroids(:,2),'.');
% end


%% Compute correlation for all angles

if verbose || saveVideo
    figure('Position',[300,300,1000,400]);
    hA = subplot(2,3,[1,4]);
    if saveVideo
        hV = VideoWriter(VideoFile); open(hV);
    end
end

if verbose || saveVideo
    rsquared = nan(180,1);
    p = nan(180,1);
    theta = 1:180;
else
    rsquared = nan(1791,1);
    p = nan(1791,1);
    theta = 1:0.1:180;
end
for tindex = 1:numel(theta)
        
    % Project data onto angle
    dist = ProjectOntoAngle(WhiteCentroids,theta(tindex));
    
    % Compute fit
    [~, order] = sort(dist);
    coeff = [ones(numel(dist),1), dist(order)]\Data(order);
    Yfit = [ones(numel(dist),1), dist(order)]*coeff;
    rsquared(tindex) = 1 - sum((Data(order) - Yfit).^2)/sum((Data - mean(Data)).^2);
    
    % Compute correlation
    [~,p(tindex)] = corr(dist,Data);
    
    if (verbose) || saveVideo  % && tindex==numel(theta)
        
%         % Check gradient aligns with axis
%         subplot(2,3,[1,4]);
%         image(zeros(ImageSize));
%         [temp,~,N] = scaleContinuousData(dist);
%         CMap = parula(N);
%         overlayROIs(Centroids, 'axes', hA, 'Radius', 10, 'Color', CMap(temp,:));
%         [x,y] = pol2cart(theta(tindex)*2*pi/360,0.5);
%         vector = [x,y;-x,-y];
%         vector = bsxfun(@plus,vector*invMat,Mean);
%         hold on; plot(vector(:,1),vector(:,2),'LineWidth',2,'Color',[1,192/255,203/255]); hold off;
        
        % Plot fit
        subplot(2,3,[2,5]);
        plot(dist, Data, 'b.');
        hold on;
        plot(dist(order), Yfit, 'r--'); hold off;
        xlim([-3,3]); ylim([4,20]);
        xlabel('Distance along axis');
        ylabel('Preferred pole position');
        text(1,1, sprintf('r^2 = %.4f', rsquared(tindex)), 'Units', 'Normalized', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top');
        
        % Plot rsquared
        subplot(2,3,3);
        plot(theta(1:tindex),rsquared(1:tindex),'b-');
        ylim([0,.5]);
        ylabel('r^2');
        xlabel('Angle');
        xlim([1,180]);
        subplot(2,3,6);
        plot(theta(1:tindex),p(1:tindex),'r-');
        ylim([0,.055]);
        ylabel('p');
        xlabel('Angle');
        xlim([1,180]);
%         yyaxis left
%         plot(theta(1:tindex),rsquared(1:tindex),'b-'); 
%         ylabel('r^2');
%         ylim([0,.5]);
%         yyaxis right
%         plot(theta(1:tindex),p(1:tindex),'r-');
%         ylabel('p');
%         xlabel('Angle'); 
%         xlim([1,180]); 
        
        % Display images
        drawnow; 
                
        % Save to video or pause display
        if saveVideo
            img = getframe(gcf); 
            writeVideo(hV,img);
        else
%             pause(.1);
        end
        
    end %verbose
end %theta

if saveVideo
    close(hV); %close video
end


%% Project data onto best angle

% Determine best angle
[~,temp] = max(rsquared);
theta = theta(temp);

% Compute projection along best fit axis, centered over PWC
if ~isempty(insidePWC)
    PWCcenter = mean(WhiteCentroids(insidePWC,:));                          % find PWC center
else
    PWCcenter = mean(WhiteCentroids);
end
[dist,Y] = ProjectOntoAngle(bsxfun(@minus,WhiteCentroids,PWCcenter),theta); % project on axis centered on PWC
Projection = [dist,Y]*invMat;                                               % transform back to non-white axis

% Calculate vector along the best angle
[x,y] = pol2cart(theta*2*pi/360,0.5);
vector = [x,y;-x,-y]*invMat;
theta = cart2pol(vector(1,1),vector(1,2))*180/pi;
vector = bsxfun(@plus,vector,Mean);


% %% Checks
% if verbose
%     
%     % Check gradient falls along axis
%     figure('Position', [50, 50, 1400, 800]);
%     hA = axes();
%     spatialOverlay(Centroids, Maps, ROIindex, FileIndex, dist, [], [], [],...
%         'DataType', 'continuous', 'Image', zeros(512,796), 'axes', hA);
%     hold on; plot(vector(:,1),vector(:,2),'LineWidth',2,'Color',[1,192/255,203/255])
%     
%     % Visually compare CoM to axis
%     figure('Position', [50, 50, 1400, 800]);
%     hA = axes();
%     spatialOverlay(Centroids, Maps, ROIindex, FileIndex, Data, [], [], [],...
%         'DataType', 'continuous', 'Image', zeros(512,796), 'axes', hA,'showColorBar','colorbarLabel','CoM (mm)');
%     hold on; plot(vector(:,1),vector(:,2),'LineWidth',2,'Color',[1,192/255,203/255])
%     
%     % [temp,~,numSamps] = scaleContinuousData(Data);
%     % img = embedROIs(Vertices(Index,1),'Color',temp,'show','CMap',parula(numSamps));
%     
%     % Check inside PWC
%     figure('Position', [50, 50, 1400, 800]);
%     hA = axes();
%     spatialOverlay(ROIs, Maps, ROIindex, FileIndex, SubInsidePWC+1, [], [0,0,1;1,0,0], [],...
%     'DataType', 'continuous', 'Image', zeros(512,796), 'axes', hA);
%     hold on; plot(vector(:,1),vector(:,2),'LineWidth',2,'Color',[1,192/255,203/255])
%     
% end
% 
% 
% %% Save to file
% if saveOut && ~isempty(saveFile)
%     for findex = 1:numFiles
%         dealROIdata(Centroids{findex}, 'CoMFitDist', Projection, 'ROIindex', ROIindex, 'Save', 'SaveFile', saveFile{findex}); % send to first file
%     end
% end


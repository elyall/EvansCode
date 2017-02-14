function [Projection,bestAngle,rsquared,theta] = projectOntoCoMAxis(Centroids, Data, varargin)

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

    
%% Whiten data

% Whiten centroids
[WhiteCentroids, invMat, ~, Mean] = whiten(Centroids);


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
    theta = 0.1:0.1:180;
    rsquared = nan(numel(theta),1);
    p = nan(numel(theta),1);
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

% De-whiten angles
% figure;
for tindex = 1:numel(theta)
    [x,y] = pol2cart(theta(tindex)*2*pi/360,50);% convert to cartesian vector
    v = [x,y]*invMat;                           % apply inverse whitening transform
    theta(tindex) = cart2pol(v(1),v(2))*180/pi; % convert to angle
%     plot([x,-x],[y,-y],'b-'); hold on;
%     plot([v(1),-v(1)],[v(2),-v(2)],'r-'); xlim([-100,100]); ylim([-100,100]); drawnow; hold off;
end


%% Project data onto best angle

% Determine best angle
[~,bestAngle] = max(rsquared);
bestAngle = theta(bestAngle);

% Determine PWC center
if ~isempty(insidePWC)
    m1 = min(Centroids(insidePWC,:));
    m2 = max(Centroids(insidePWC,:));
else
    m1 = min(Centroids);
    m2 = max(Centroids);
end
PWCcenter = (m2-m1)/2+m1;

% Compute projection along best fit axis, centered over PWC
[dist,Y] = ProjectOntoAngle(Centroids,bestAngle,PWCcenter); % project on axis centered on PWC
Projection = [dist,Y];



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


function [insidePWC, ROIindex, FileIndex] = determinePWC(ROIs, ExperimentFiles, varargin)

ROIindex = [1 inf];
FileIndex = [1 inf];

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
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

if ~exist('ROIs','var') || isempty(ROIs)
    [ROIs, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory,'MultiSelect','on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p,ROIs);
end

if ~exist('ExperimentFiles','var') || isempty(ExperimentFiles)
    [ExperimentFiles, p] = uigetfile({'*.exp;*.mat'},'Choose Experiment file',directory,'MultiSelect','on');
    if isnumeric(ExperimentFiles)
        return
    end
    ExperimentFiles = fullfile(p,ExperimentFiles);
end
if ~iscell(ExperimentFiles)
    ExperimentFiles = {ExperimentFiles};
end
numFiles = numel(ExperimentFiles);


%% Load in data

% Load in centroids
[Centroids, FileIndex, ROIindex, ~, ~, ROIs] = gatherROIdata(ROIs, 'centroid', ':', 'none', ROIindex, FileIndex);

% Load in StimResponse
Images = cell(numFiles, 1);
Maps = imref2d;
for findex = 1:numFiles
    load(ExperimentFiles{findex}, 'StimResponse', '-mat');
%     Images{findex} = max(StimResponse.avg,[],4);
    Images{findex} = StimResponse.pvalue(:,:,1,7);
    if numFiles > 1
        load(ExperimentFiles{findex}, 'Map', '-mat');
        Maps(findex) = Map;
    end
end


%% Build up Map and translate ROIs
if numFiles > 1
    offsets = mapFoVs(Maps, 'Type', 'index');
    
    % Translate ROIs
    for findex = 1:numFiles
        Centroids(FileIndex==findex,:) = bsxfun(@plus, Centroids(FileIndex==findex,:), offsets(findex, 1:2));
    end
end


%% UI-select PWC
Image = createImage(Images, Maps);
hF = figure('NumberTitle', 'off', 'Name', 'Select Principle Whisker Column');
imagesc(Image);
hP = impoly(gca);
for index = 10:-1:1
    set(hF, 'Name', sprintf('Closing in %d seconds', index));
    pause(1);
end
PWBoundary = getPosition(hP);
PWBoundary = [PWBoundary; PWBoundary(1,:)];
close(gcf);


%% Determine which ROIs are in the PWC
insidePWC = logical(inpolygon(Centroids(:,1), Centroids(:,2), PWBoundary(:,1), PWBoundary(:,2)));

% Check: overlay ROIs
figure; imagesc(Image); hold on; plot(PWBoundary(:,1), PWBoundary(:,2), 'w', 'LineWidth', 2);
hold on;
plot(Centroids(insidePWC==1,1), Centroids(insidePWC==1,2), 'r*'); 
plot(Centroids(insidePWC==0,1), Centroids(insidePWC==0,2), 'b*'); 
plot(PWBoundary(:,1), PWBoundary(:,2), 'g', 'LineWidth', 2);




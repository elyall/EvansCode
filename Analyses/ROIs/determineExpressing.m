function [likelyExpressing,unknownPhenotype,meanFluorescence,ROIindex] = determineExpressing(Images, ROIdata, varargin)

ROIindex = [1 inf];         % vector of ROIs to analyze
MotionCorrect = false;      % false, filename to load MCdata from, or true to prompt for file selection
FrameIndex = 2:501;         % vector of frame indices
borderLims = [0,0,32,32];   % number of pixels to remove from edges when computing ROI means (top, bottom, left, right)
% Channel = 2;                % channel to extract data from

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'ROIs','ROI','ROIindex'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'MotionCorrect'
                MotionCorrect = varargin{index+1};
                index = index + 2;
            case {'Frames', 'frames', 'FrameIndex'}
                FrameIndex = varargin{index+1};
                index = index + 2;
            case {'borderLims','Border','border'}
                borderLims = varargin{index+1};
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

if ~exist('Images', 'var') 
    Images = [];
end

if ~exist('ROIdata', 'var')
    ROIdata = [];
end


%% Determine which ROIs are putatively expressing
[~, GData] = extractSignals(Images, ROIdata, ROIindex,...
    'Frames',           FrameIndex,...
    'MotionCorrect',    MotionCorrect,...
    'Channel',          1,...
    'Border',           borderLims);
[~, RData] = extractSignals(Images, ROIdata, ROIindex,...
    'Frames',           FrameIndex,...
    'MotionCorrect',    MotionCorrect,...
    'Channel',          2,...
    'Border',           borderLims);
% RData = TrueRed + k*GData (leak)
likelyExpressing = min(RData./GData,[],2)>1;        % the red can't be solely explained by leak from the green channel
meanFluorescence = mean(RData(:,FrameIndex),2);     % average over time
threshold = prctile(meanFluorescence,40);           % threshold should be slightly above noise floor from PMT
likelyExpressing = all([likelyExpressing,meanFluorescence>threshold],2); % red is likely not explained by noise either


%% Determine which ROIs have an unknown phenotype
unknownPhenotype = all([meanFluorescence>threshold,~likelyExpressing],2);

% gid = kmeans(meanFluorescence,2);
% temp = nan(1,2);
% for gindex = 1:2
%     temp(gindex) = mean(meanFluorescence(gid==gindex));
% end
% [~,id] = max(temp);
% unknowns = gid==id;
% unknowns(putativelyExpressing) = false;
% 
% figure; hold on;
% [~,edges] = histcounts(meanFluorescence);
% for gindex=1:2
%     histogram(meanFluorescence(gid==gindex),edges);
% end



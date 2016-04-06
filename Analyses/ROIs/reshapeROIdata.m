function Out = reshapeROIdata(ROIdata, varargin)

ROIindex = [1 inf];
TrialIndex = [1 inf];
FrameIndex = [];

type = 'MatMean';
dataType = 'dFoF';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'FrameIndex', 'Frames', 'frames'}
                FrameIndex = varargin{index+1};
                index = index + 2;
            case {'TrialIndex', 'Trials', 'trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'dataType'
                dataType = varargin{index+1};
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

if ~exist('ROIdata', 'var')
    [ROIdata, p] = uigetfile({'*.rois;*.mat'}, 'Select ROI file', directory);
    if ~ROIdata
        return
    else
        ROIdata = fullfile(p, ROIdata);
    end
end


%% Load file
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
end


%% Determine parameters
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end

if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:numel(ROIdata.DataInfo.StimID)];
end

if isempty(FrameIndex) % default is stim period
    FrameIndex = ROIdata.DataInfo.numFramesBefore+1:ROIdata.DataInfo.numFramesBefore+mode(ROIdata.DataInfo.numStimFrames(TrialIndex));
elseif FrameIndex(end) == inf
    FrameIndex = [FrameIndex(1:end-1):FrameIndex(end-1)+1:ROIdata.DataInfo.numFramesBefore+1+ROIdata.DataInfo.numFramesAfter];
end


%% Reshape ROIdata
[T,F] = size(ROIdata.rois(1).(dataType));
Out = [ROIdata.rois(:).(dataType)];             % pull out data
Out = reshape(Out,T,F,numel(ROIdata.rois));     % reshape it
Out = permute(Out,[1,3,2]);                     % permute it

Out = Out(TrialIndex,ROIindex,FrameIndex);      % keep only data requested

switch type
    case 'MatMean'
        Out = mean(Out,3);
    case 'MatMedian'
        Out = median(Out,3);
end


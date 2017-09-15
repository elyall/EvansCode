function [Data, StimID] = ROIs2AvgMat(ROIs, varargin)


FrameIndex = [];
TrialIndex = [];
ROIindex = [];
FileIndex = [];

saveOut = false;
saveFile = '';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'FrameIndex', 'Frames', 'frames'}
                FrameIndex = varargin{index+1};
                index = index + 2;
            case {'TrialIndex', 'Trials', 'trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
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

if ~exist('ROIs', 'var')
    [ROIs, p] = uigetfile({'*.rois;*.mat'},'Select ROI file',directory,'MultiSelect','on');
    if ~ROIs
        return
    end
    ROIs = fullfile(p, ROIs);
end


%% Grab trial data
Data = gatherROIdata(ROIs, 'dFoF', 'ROIindex', ROIindex, 'FileIndex', FileIndex); % pull out data


%% Perform average and keep only specified trials

% Determine parameters
if isempty(TrialIndex)
    TrialIndex = 1:numel(ROIs{1}.DataInfo.StimID);
elseif TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:numel(ROIs{1}.DataInfo.StimID)]; % set trials to save
end
if isempty(FrameIndex)
    FrameIndex = [ROIs{1}.DataInfo.numFramesBefore+1, ROIs{1}.DataInfo.numFramesBefore+mode(ROIs{1}.DataInfo.numStimFrames(TrialIndex))];
end

% Perform average across specified frames
Data = mean(Data(:,TrialIndex,FrameIndex(1):FrameIndex(2)), 3);
StimID = ROIs{1}.DataInfo.StimID(TrialIndex)'; % pull out StimID


%% Save output
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'Data', 'StimID', '-mat', '-v7.3');
    else
        save(saveFile, 'Data', 'StimID', '-mat', '-append');
    end
end


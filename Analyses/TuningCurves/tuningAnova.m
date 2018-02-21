function [PValues, ROIdata] = tuningAnova(Data, varargin)


StimIndex = [2 inf];
ROIindex = [1 inf];

StimShape = [];

saveOut = false;
saveFile = '';


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'StimShape'
                StimShape = varargin{index+1};
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

if ~exist('Data', 'var') || isempty(Data)
    [Data,p] = uigetfile({'*.rois;*.mat'}, 'Select ROI files:', directory);
    if isnumeric(Data)
        return
    end
    Data = {fullfile(p, Data)};
end


%% Load data
if ischar(Data)
    ROIFile = Data;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
    Data = ROIdata;
elseif isstruct(Data)
    ROIdata = Data;
    Data = gatherROIdata(ROIdata,'Raw');
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end


%% Determine data to analyze
[numROIs,numStims] = size(Data);

% Determine ROIs to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numROIs];
end
numROIs = numel(ROIindex);
Data = Data(ROIindex,:);

% Determine stimuli to analyze
if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1), StimIndex(end-1)+1:numStims];
end
numStims = numel(StimIndex);
Data = Data(:,StimIndex);


%% Compute anova

N = cellfun(@numel,Data); % determine # of trials for each stimulus

% Compute significance
PValues = nan(numROIs, 1+numel(StimShape));
% fid = parfor_progress(numROIs);
parfor rindex = 1:numROIs
    
    try
        current = cat(1,Data{rindex,:});          % gather data
    catch
        current = cat(2,Data{rindex,:})';         % gather data
    end
    group = repelem(1:numStims,N(rindex,:))'; % define grouping (each stimulus is unique)
    if ~isempty(StimShape) % create other groupings
        group = repmat(group,1,3);
        [X,Y] = meshgrid(1:StimShape(2),1:StimShape(1));
%         X = X+numStims;  % offset group ID
%         Y = Y+max(X(:)); % offest group ID
        group(:,2) = repelem(X(:),N(rindex,:)); % each column is unique
        group(:,3) = repelem(Y(:),N(rindex,:)); % each row is unique
    end
%     PValues(rindex,:) = anovan(current, group, 'model', 'linear', 'display', 'off'); % compute N-way ANOVA (DOESN'T WORK AS ADVERTISED)
    p = nan(1,1+numel(StimShape));
    for gindex = 1:size(group,2)
        p(gindex) = anova1(current, group(:,gindex), 'off');
    end
    PValues(rindex,:) = p;
    
%     parfor_progress(fid);
end
% parfor_progress(fid,0);


%% Distribute to struct
if nargin > 1 || saveOut
    for rindex = 1:numROIs
        ROIdata.rois(ROIindex(rindex)).TunedPValue = PValues(rindex);
    end
end


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end


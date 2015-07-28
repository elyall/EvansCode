function [AnalysisInfo, frames] = behaviorPostProcess2(info, DataIn, data, Config, varargin)

numFramesBefore = 0; % numeral
saveInputData = false;
saveOut = false;
SaveFile = '';
samplingFrequency = 30000;

% Default
existInputData = false;
seriesType = {};
seriesNames = {};

%% Parse input arguments
if ~exist('info', 'var') || isempty(info)
    directory = cd;
    [info,p] = uigetfile({'*.mat'}, 'Choose scanbox info file', directory);
    if isnumeric(info)
        return
    end
    info = fullfile(p,info);
end

if ~exist('DataIn', 'var') || isempty(DataIn)
    directory = cd;
    [DataIn,p] = uigetfile({'*.bin'}, 'Choose NI-DAQ DataIn File', directory);
    if isnumeric(DataIn)
        return
    end
    DataIn = fullfile(p,DataIn);
end

if ~exist('data', 'var') || isempty(data)
    directory = cd;
    [data,p] = uigetfile({'*.txt'}, 'Choose stimulation text file', directory);
    if isnumeric(data)
        return
    end
    data = fullfile(p,data);
end

if ~exist('Config', 'var') || isempty(Config)
    directory = cd;
    [Config,p] = uigetfile({'*.sbx;*.tiff;*.imgs'}, 'Choose images file', directory, 'MultiSelect', 'on');
    if isnumeric(Config)
        return
    elseif iscell(Config)
        for findex = 1:numel(Config)
            Config{findex} = fullfile(p, Config);
        end
    elseif ischar(Config)
        Config = {fullfile(p, Config)};
    end
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numFramesBefore'
                numFramesBefore = varargin{index+1};
                index = index + 2;
            case 'Save'
                saveOut = true;
                index = index + 1;
            case 'SaveFile'
                SaveFile = varargin{index+1};
                index = index + 2;
            case 'seriesType'
                seriesType = varargin{index+1};
                index = index + 2;
            case 'seriesNames'
                seriesNames = varargin{index+1};
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

if saveOut && isempty(SaveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end

%% Load in data

% Image info
if ischar(info)
    InfoFile = info;
    load(InfoFile, 'info');
end
if ischar(Config)
    Config = load2PConfig(Config);
end
numFrames = sum(Config(:).Frames);
FrameRate = Config(1).FrameRate;

% NI-DAQ data
if ischar(DataIn)
    existInputData = true;
    FID = fopen(DataIn, 'r');
    DataIn = fread(FID, 'uint16');
    DataIn = reshape(DataIn, 9, numel(DataIn)/9)';
end
[numScans, numSeriesVariables] = size(DataIn);
if isempty(seriesNames)
    seriesNames = cellstr(strcat('var',num2str([1:numSeriesVariables]')));
end
if isempty(seriesType)
    seriesType = repmat({'digital'}, numSeriesVariables, 1);
end

% Stim info
if ischar(data)
    TextFile = data;
    data = csvread(TextFile);
end
% numScans = size(data, 1);


%% Determine number of trials
stimulusIndex = DataIn(:,strcmp(seriesNames, 'TriggerPeriod'));
stimulusIndex = stimulusIndex - [0;stimulusIndex(1:end-1)];
stimulusIndex = find(stimulusIndex);
% nTrials = numel(find(stimulusIndex>0));

[TrialIndex, temp] = unique(data(:,6));
TrialIndex(1) = []; temp(1) = []; % remove first trial
nTrials = numel(TrialIndex);

nImagingTrials = numel(info.frame);


%% Initialize output
AnalysisInfo = table(nan(nTrials,1),nan(nTrials,1),nan(nTrials,1),cell(nTrials,1),nan(nTrials,1),nan(nTrials,2),nan(nTrials,2),nan(nTrials,2),nan(nTrials,2),nan(nTrials,2),...
    'VariableNames', {'StimID', 'TrialIndex', 'ImgIndex', 'ImgFilename', 'nFrames', 'ExpStimFrames', 'ExpFrames', 'StimFrameLines', 'EventIDs', 'TrialStimFrames'});
frames = struct('Stimulus', nan(numFrames,1), 'Trial', nan(numFrames,1));

if existInputData
    AnalysisInfo = [AnalysisInfo, table(nan(nTrials,1), nan(nTrials,2), nan(nTrials,2), nan(nTrials,2),...
        'VariableNames',{'nScans', 'ExpStimScans', 'ExpScans', 'TrialStimScans'})];
    if saveInputData
        AnalysisInfo = [AnalysisInfo, table(cell(nTrials,1),'VariableNames',{'DataIn'})];
    end
    for sindex = 1:numSeriesVariables
        if any(strcmp(seriesType{sindex}, {'analog', 'digital', 'counter', 'logical'}))
            frames.(seriesNames{sindex}) = nan(numFrames, 1);
        end
    end
end


%% Process and parse experiment

% AnalysisInfo.ImgFilename = ImageFile;
AnalysisInfo.TrialIndex = TrialIndex;
AnalysisInfo.StimID = data(temp, 11);

for tindex = 1:nTrials
    
    % Imaging data
    if tindex <= nImagingTrials %assumes first trial is recorded but possibly not last trial(s)
        AnalysisInfo.ImgIndex(tindex) = tindex;
        
        AnalysisInfo.ExpStimFrames(tindex, :) = info.frame(2*(tindex-1)+1:2*tindex); % two TTLs: one at beginning and one at end
        AnalysisInfo.StimFrameLines(tindex, :) = info.line(2*(tindex-1)+1:2*tindex);
        AnalysisInfo.EventIDs(tindex, :) = info.event_id(2*(tindex-1)+1:2*tindex);
%         AnalysisInfo.ExpStimFrames(tindex, :) = [info.frame(tindex), info.frame(tindex)+round(1*FrameRate)]; % one TTL: only at beginning
%         AnalysisInfo.StimFrameLines(tindex, 1) = info.line(tindex);
%         AnalysisInfo.EventIDs(tindex, 1) = info.event_id(tindex);

        frames.Stimulus(AnalysisInfo.ExpStimFrames(tindex,1):AnalysisInfo.ExpStimFrames(tindex,2)) = AnalysisInfo.StimID(tindex);
        AnalysisInfo.TrialStimFrames(tindex, :) = [numFramesBefore + 1, numFramesBefore + diff(AnalysisInfo.ExpStimFrames(tindex, :)) + 1];
        AnalysisInfo.ExpFrames(tindex, 1) = AnalysisInfo.ExpStimFrames(tindex, 1) - numFramesBefore;
        if tindex > 1
            AnalysisInfo.ExpFrames(tindex-1, 2) = AnalysisInfo.ExpFrames(tindex, 1) - 1;
            AnalysisInfo.nFrames(tindex-1) = diff(AnalysisInfo.ExpFrames(tindex-1, :))+1;
            frames.Trial(AnalysisInfo.ExpFrames(tindex-1,1):AnalysisInfo.ExpFrames(tindex-1,2)) = tindex-1;
        end
    end
    if tindex == nImagingTrials
        AnalysisInfo.ExpFrames(tindex, 2) = numFrames;
        AnalysisInfo.nFrames(tindex) = diff(AnalysisInfo.ExpFrames(tindex, :))+1;
        frames.Trial(AnalysisInfo.ExpFrames(tindex,1):AnalysisInfo.ExpFrames(tindex,2)) = tindex;
    end
    
    % Behavioral Data
    if existInputData
        AnalysisInfo.ExpStimScans(tindex, :) = stimulusIndex(2*(tindex-1)+1:2*tindex) + [0;-1];
        AnalysisInfo.TrialStimScans(tindex, :) = [1, diff(AnalysisInfo.ExpStimScans(tindex, :))+1];
        if tindex > 1
            AnalysisInfo.ExpScans(tindex-1, :) = [AnalysisInfo.ExpStimScans(tindex-1, 1), AnalysisInfo.ExpStimScans(tindex, 1)-1];
            AnalysisInfo.nScans(tindex-1) = diff(AnalysisInfo.ExpScans(tindex-1,:)) + 1;
        end
        if tindex == nTrials
            AnalysisInfo.ExpScans(tindex, :) = [AnalysisInfo.ExpStimScans(tindex, 1), numScans];
            AnalysisInfo.nScans(tindex) = diff(AnalysisInfo.ExpScans(tindex,:)) + 1;
        end
        if saveInputData
            AnalysisInfo.saveInputData{tindex} = DataIn(AnalysisInfo.ExpScans(tindex, 1):AnalysisInfo.ExpScans(tindex, 2),:);
        end
    end
    
end

% Parse continuous data to be frame-wise
if existInputData
    theoreticalScansPerFrame = samplingFrequency/FrameRate;
    for tindex = 1:nImagingTrials
        %     hF = figure; hold on
        if tindex ~= nImagingTrials
            scansPerFrame = AnalysisInfo.nScans(tindex)/AnalysisInfo.nFrames(tindex); % miniscule difference between "samplingFrequency/FrameRate"
        else % allows for imaging to cut off before NI-DAQ
            scansPerFrame = theoreticalScansPerFrame;
        end
        if abs(scansPerFrame - theoreticalScansPerFrame)/theoreticalScansPerFrame > 0.05*theoreticalScansPerFrame % if number of scans differs by 5% assume imaging stopped in the middle of a trial
            scansPerFrame = theoreticalScansPerFrame;
        end
        for findex = 1:AnalysisInfo.nFrames(tindex)
            scans = AnalysisInfo.ExpScans(tindex,1)+round(scansPerFrame*(findex-1)):min(AnalysisInfo.ExpScans(tindex,1)+round(scansPerFrame*findex), numScans);
            for sindex = 1:numSeriesVariables
                switch seriesType{sindex}
                    case 'analog'
                        frames.(seriesNames{sindex})(AnalysisInfo.ExpFrames(tindex,1)+findex-1) = mean(DataIn(scans, sindex));
                    case 'digital'
                        frames.(seriesNames{sindex})(AnalysisInfo.ExpFrames(tindex,1)+findex-1) = sum(DataIn(scans, sindex));
                    case 'logical'
                        frames.(seriesNames{sindex})(AnalysisInfo.ExpFrames(tindex,1)+findex-1) = any(DataIn(scans, sindex));
                    case 'counter'
                        frames.(seriesNames{sindex})(AnalysisInfo.ExpFrames(tindex,1)+findex-1) = DataIn(scans(end), sindex);
                end
            end %series
        end %frames
        %     fprintf('\nEnded on scan %d rather than scan %d: difference of %.2f seconds', scans(end), AnalysisInfo.ExpScans(tindex,2),(AnalysisInfo.ExpScans(tindex,2)-scans(end))/samplingFrequency);
        %     for sindex = 1:numSeriesVariables
        %         plot(frames.(seriesNames{findex})(AnalysisInfo.ExpFrames(tindex,1):AnalysisInfo.ExpFrames(tindex,2)));
        %     end
    end %trials
end

if saveOut
    save(SaveFile, 'AnalysisInfo', 'frames', '-append');
end

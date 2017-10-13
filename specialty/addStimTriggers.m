function info = addStimTriggers(info, ExperimentFile, FirstID)

frameTriggerScanOffset = false;
saveOut = true;

%% Grab experiments and look at other first frame IDs
% [~,list]=system('find /media/elyall/Data/ -type f -name "*.exp"');
% c=strsplit(list,'\n');
% c(find(cellfun(@isempty,strfind(c,'.exp')))) = [];
% c(find(cellfun(@isempty,strfind(c,'/17')))) = [];
% c(find(~cellfun(@isempty,strfind(c,'5701')))) = [];
% 
% c = cellfun(@(x) [x(1:end-3),'mat'],c,'UniformOutput',false);
% firstFrame = nan(numel(c),1);
% for findex = 1:numel(c)
%     try
%         load(c{findex},'info','-mat')
%         firstFrame(findex) = info.frame(1);
%     end
% end


%% Load DataIn
DataInFile = [ExperimentFile(1:end-3),'bin'];
if ~exist(DataInFile,'file')
    error('Cannot locate DataInFile: %s',DataInFile);
end
load(ExperimentFile,'DAQChannels','Experiment','TrialInfo','numDelayScans','-mat');
if ~exist('numDelayScans','var')
    numDelayScans = 0;
end
InputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
DataInFID = fopen(DataInFile, 'r');
DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
fclose(DataInFID);
nInputChannels = numel(InputNames);
DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';    % reshape to be nScans x nInputChannels
DataIn = DataIn(:,strcmp(InputNames,'I_FrameCounter'));                     % keep only frame trigger
DataIn = DataIn-[0;DataIn(1:end-1)]>0;                                      % locate onset of trigger
if frameTriggerScanOffset<0                                                 % offset frame trigger to coincide with start of frame
    DataIn = cat(1,DataIn(abs(frameTriggerScanOffset)+1:end),false(abs(frameTriggerScanOffset),1));
elseif frameTriggerScanOffset>0
    DataIn = cat(1,false(frameTriggerScanOffset,1),DataIn(1:end-frameTriggerScanOffset));
end

%% Determine scan ID of each trigger output
numTrials = numel(TrialInfo.StimID);
OutputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));

% Create output triggers
Trigger = false(numDelayScans,1);
Trial = logical(Experiment.Triggers(:,strcmp(OutputNames,'O_EventTrigger'),1));
for tindex = 1:numTrials
    Trigger = cat(1,Trigger,Trial);
    if isfield(TrialInfo,'numRandomScansPost')
        Trigger = cat(1,Trigger,false(TrialInfo.numRandomScansPost(tindex),1));
    end
end

% Determine scan ID of each trigger
TriggerScans = find(Trigger); 
% clear Trigger;


%% Determine which frame each trigger occurred on
if ischar(info)
    load(info,'info','-mat')
end

% Determine # of triggers
numTrigsPerTrial = nnz((Trial-[0;Trial(1:end-1)])>0);
numTriggers = numTrials*numTrigsPerTrial;

% Initialize output
info.frame = nan(numTriggers,1);
info.line = nan(numTriggers,1);
info.event_id = ones(numTriggers,1);

% Determine frame ID of each trigger
if numDelayScans==0
    info.frame(1) = FirstID; % set first frame of first stimulus from user input
else % assumes holdStart is active -> all frames are accounted for
    info.frame(1) = sum(DataIn(1:TriggerScans(1)));
end
for tindex = 2:numTriggers
    info.frame(tindex) = sum(DataIn(TriggerScans(tindex-1)+1:TriggerScans(tindex)))+info.frame(tindex-1); % add # of frames initiated since last trigger
end
info.frame = info.frame - 1; % change 1 indexing to 0 indexing


%% Save to file
if saveOut
    save([ExperimentFile(1:end-3),'mat'],'info','-mat');
    fprintf('Saved ''info'' to: %s\n',[ExperimentFile(1:end-3),'mat']);
end



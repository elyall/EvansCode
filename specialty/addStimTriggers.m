function info = addStimTriggers(info, ExperimentFile, FirstID)

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
load(ExperimentFile,'DAQChannels','Experiment','TrialInfo','-mat');
InputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
DataInFID = fopen(DataInFile, 'r');
DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
fclose(DataInFID);
nInputChannels = numel(InputNames);
DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';    % reshape to be nScans x nInputChannels
DataIn = DataIn(:,strcmp(InputNames,'I_FrameCounter'));                     % keep only frame trigger
DataIn = DataIn-[0;DataIn(1:end-1)]>0;                                      % locate onset of trigger


%% Determine scan ID of each trigger output
numTrials = numel(TrialInfo.StimID);
OutputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));

% Create output triggers
Trigger = logical([]);
Trial = logical(Experiment.Triggers(:,strcmp(OutputNames,'O_2PTrigger'),1));
for tindex = 1:numTrials
    Trigger = cat(1,Trigger,Trial);
    if isfield(TrialInfo,'numRandomScansPost')
        Trigger = cat(1,Trigger,false(TrialInfo.numRandomScansPost(tindex),1));
    end
end

% Determine scan ID of each trigger
ScanID = find(Trigger); clear Trigger;


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
info.frame(1) = FirstID; % set first frame from user input
for tindex = 2:numTriggers
    info.frame(tindex) = sum(DataIn(ScanID(tindex-1):ScanID(tindex)))+info.frame(tindex-1);
end

%% Save to file
if saveOut
    save([ExperimentFile(1:end-3),'mat'],'info','-mat');
    fprintf('Saved ''info'' to: %s\n',[ExperimentFile(1:end-3),'mat']);
end



function [TrialInfo, SaveFile] = stimNubulate(SaveFile, varargin)


%% Configure Settings
gd.Internal.ImagingType = 'sbx';                % 'sbx' or 'scim'
gd.Internal.ImagingComp.ip = '128.32.173.30';   % SCANBOX ONLY: for UDP
gd.Internal.ImagingComp.port = 7000;            % SCANBOX ONLY: for UDP
% gd.Internal.LinearStage.port = 'com3';

gd.Internal.numStimuliRepetitions = 10;          % default number of times to repeat each stimulus (can be changed via GUI)


Display.units = 'pixels';
Display.position = [200, 200, 1000, 600];

%% Initialize Experiment

% Saving info
gd.Internal.save.path = 'C:\Users\Resonant-2\OneDrive\StimData';
gd.Internal.save.base = '0000';
gd.Internal.save.depth = '000';
gd.Internal.save.index = '000';

% Experiment info
gd.Experiment.saving.save = false;
gd.Experiment.saving.SaveFile = fullfile(gd.Internal.save.path, gd.Internal.save.base);
gd.Experiment.saving.DataFile = '';
gd.Experiment.saving.dataPrecision = 'uint16';

gd.Experiment.params.samplingFrequency = 30000;
gd.Experiment.params.numTrials = 180;
gd.Experiment.params.frameRateWT = 125;

gd.Experiment.timing.stimDuration = .75; % in seconds
gd.Experiment.timing.ITI = .25; % in seconds

gd.Experiment.stim.setup = table(...
    {'C1';'C2';'B1';'D1';'beta';'gamma';'';''},...
    {'port0/line8';'port0/line9';'port0/line10';'port0/line11';'port0/line12';'port0/line13';'port0/line14';'port0/line15'},...
    [true;true;true;true;true;false;false;false],...
    'VariableNames',{'Name','Port','Active'});
gd.Experiment.stim.pistonCombinations = {};
gd.Experiment.stim.control = false; % true or false (give control stimulus)
gd.Experiment.stim.blockShuffle = false; % true or false (shuffle trials)
gd.Experiment.stim.repeatBadTrials = true; % true or false (repeat non-running trials)
gd.Experiment.stim.speedThreshold = 100; % speed threshold for good running trials (deg/s)

% Properties for display or processing input data
gd.Internal.buffer.numTrials = 5; %4*gd.Experiment.params.samplingFrequency * (gd.Experiment.timing.stimDuration+gd.Experiment.timing.ITI);
gd.Internal.buffer.downSample = 20;

% Place holder for offline control of stimulus
gd.Internal.daq = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'stimDuration'
                gd.Experiment.timing.stimDuration = varargin{index+1};
                index = index + 2;
            case 'ITI'
                gd.Experiment.timing.ITI = varargin{index+1};
                index = index + 2;
            case 'control'
                gd.Experiment.stim.control = varargin{index+1};
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

if exist('SaveFile', 'var')
    gd.Experiment.saving.SaveFile = SaveFile;
    gd.Internal.save = true;
end


%% Create & Populate Figure

% Create figure
gd.fig = figure(...
    'NumberTitle',          'off',...
    'Name',                 'Stimulus: Nubulate',...
    'ToolBar',              'none',...
    'Units',                Display.units,...
    'Position',             Display.position);

% SAVING DATA
% panel
gd.Saving.panel = uipanel(...
    'Title',                'Save File',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [0, .8, 1, .2]);
% save button
gd.Saving.save = uicontrol(...
    'Style',                'radiobutton',...
    'String',               'Save?',...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,.1,1],...
    'BackgroundColor',      [1,0,0],...
    'Value',                gd.Experiment.saving.save,...
    'Callback',             @(hObject,eventdata)toggleSave(hObject, eventdata, guidata(hObject)));
% directory selection
gd.Saving.dir = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Dir',...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.1,.4,.1,.6],...
    'Callback',             @(hObject,eventdata)ChooseDir(hObject, eventdata, guidata(hObject)));
gd.Saving.FullFilename = uicontrol(...
    'Style',                'text',...
    'String',               '',...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.2,.05,.8,.2],...
    'Callback',             @(hObject,eventdata)CreateFilename(hObject, eventdata, guidata(hObject)));
% image acq type
gd.Saving.imagingType = uicontrol(...
    'Style',                'popupmenu',...
    'String',               {'scim';'sbx'},...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.1,0,.1,.4]);
switch gd.Internal.ImagingType
    case 'scim'
        gd.Saving.imagingType.Value = 1;
    case 'sbx'
        gd.Saving.imagingType.Value = 2;
end
% basename input
gd.Saving.base = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Internal.save.base,...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.3,.3,.2,.5],...
    'Callback',             @(hObject,eventdata)SetFilename(hObject, eventdata, guidata(hObject)));
gd.Saving.baseText = uicontrol(...
    'Style',                'text',...
    'String',               'Basename',...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.3,.8,.2,.2]);
% depth input
gd.Saving.depth = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Internal.save.depth,...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.3,.2,.5],...
    'Callback',             @(hObject,eventdata)SetDepth(hObject, eventdata, guidata(hObject)));
gd.Saving.depthText = uicontrol(...
    'Style',                'text',...
    'String',               'Depth',...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.8,.2,.2]);
% file index
gd.Saving.index = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Internal.save.index,...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.7,.3,.2,.5],...
    'Callback',             @(hObject,eventdata)SetFileIndex(hObject, eventdata, guidata(hObject)));
gd.Saving.indexText = uicontrol(...
    'Style',                'text',...
    'String',               'File Index',...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.7,.8,.2,.2]);

% Stimuli
% panel
gd.Stimuli.panel = uipanel(...
    'Title',                'Stimuli',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [0, 0, .5, .8]);
% ports list
numPorts = size(gd.Experiment.stim.setup,1);
gd.Stimuli.ports = uitable(...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'Normalized',...
    'Position',             [0,.2,.6,.8],...
    'ColumnName',           {'Name','Port','Active','Toggle'},...
    'ColumnFormat',         {'char','char','logical','logical'},...
    'ColumnEditable',       [true,true,true,true],...
    'Data',                 [table2cell(gd.Experiment.stim.setup), mat2cell(false(numPorts, 1), ones(numPorts,1), 1)],...
    'CellEditCallback',     @(hObject,eventdata)EditPorts(hObject, eventdata, guidata(hObject)));
% basis for combination creation
gd.Stimuli.editCombinations = uicontrol(...
    'Style',                'edit',...
    'String',               num2str(1:nnz(gd.Experiment.stim.setup.Active)),...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [0,.1,.4,.1],...
    'Callback',             @(hObject,eventdata)EditCombinations(hObject, eventdata, guidata(hObject)));
% create port combinations
gd.Stimuli.generateCombinations = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Generate combinations',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,.4,.1],...
    'Callback',             @(hObject,eventdata)GenerateCombinations(hObject, eventdata, guidata(hObject)));
% load combinations
gd.Stimuli.load = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Load',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.4,.1,.2,.1],...
    'Callback',             @(hObject,eventdata)LoadStimuli(hObject, eventdata, guidata(hObject)));
% save combinations
gd.Stimuli.save = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Save',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.4,0,.2,.1],...
    'Callback',             @(hObject,eventdata)SaveStimuli(hObject, eventdata, guidata(hObject)));
% stimuli list
gd.Stimuli.list = uitable(...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'Normalized',...
    'Position',             [.6,0,.4,1],...
    'ColumnName',           {'Combination','Delete'},...
    'ColumnFormat',         {'char','logical'},...
    'ColumnEditable',       [true,true],...
    'CellEditCallback',     @(hObject,eventdata)EditStimuli(hObject, eventdata, guidata(hObject)));

% EXPERIMENT
% panel
gd.Run.panel = uipanel(...
    'Title',                'Run Experiment',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [.5, 0, .5, .8]);
% length of inter-trial interval
gd.Run.ITI = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.timing.ITI,...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.25,.9,.25,.1]);
gd.Run.ITIText = uicontrol(...
    'Style',                'text',...
    'String',               'ITI duration',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlignment',  'right',...
    'Units',                'normalized',...
    'Position',             [0,.925,.25,.05]);
% length of stimulus
gd.Run.stimDur = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.timing.stimDuration,...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.25,.8,.25,.1]);
gd.Run.stimDurText = uicontrol(...
    'Style',                'text',...
    'String',               'stim duration',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlignment',  'right',...
    'Units',                'normalized',...
    'Position',             [0,.825,.25,.05]);
% number of trials
gd.Run.numTrials = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.numTrials,...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.25,.7,.25,.1]);
gd.Run.numTrialsText = uicontrol(...
    'Style',                'text',...
    'String',               '# Trials',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlignment',  'right',...
    'Units',                'normalized',...
    'Position',             [0,.725,.25,.05]);
% control toggle
gd.Run.control = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Control Trial?',...
    'Parent',               gd.Run.panel,...
    'Value',                gd.Experiment.stim.control,...
    'Units',                'normalized',...
    'Position',             [.65,.9,.3,.1]);
% block shuffle toggle
gd.Run.shuffle = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Block Shuffle?',...
    'Parent',               gd.Run.panel,...
    'Value',                gd.Experiment.stim.blockShuffle,...
    'Units',                'normalized',...
    'Position',             [.65,.8,.3,.1]);
% repeat bad trials toggle
gd.Run.repeatBadTrials = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Repeat bad trials?',...
    'Parent',               gd.Run.panel,...
    'Value',                gd.Experiment.stim.repeatBadTrials,...
    'Units',                'normalized',...
    'Position',             [.65,.7,.3,.1]);
% run button
gd.Run.run = uicontrol(...
    'Style',                'togglebutton',...
    'String',               'Run?',...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [0,.55,1,.15],...
    'Callback',             @(hObject,eventdata)RunExperiment(hObject, eventdata, guidata(hObject)));
% running speed axes
gd.Run.runSpeedAxes = axes(...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,1,.55]);

guidata(gd.fig, gd); % save guidata
CreateFilename(gd.Saving.FullFilename, [], gd);
end


%% SAVING CALLBACKS

function ChooseDir(hObject, eventdata, gd)
temp = uigetdir(gd.Internal.save.path, 'Choose directory to save to');
if ischar(temp)
    gd.Internal.save.path = temp;
    guidata(hObject, gd);
end
CreateFilename(gd.Saving.FullFilename, [], gd);
end

function SetFilename(hObject, eventdata, gd)
gd.Internal.save.base = hObject.String;
CreateFilename(gd.Saving.FullFilename, [], gd);
end

function SetDepth(hObject, eventdata, gd)
if numel(hObject.String)>3
    hObject.String = hObject.String(1:3);
end
CreateFilename(gd.Saving.FullFilename, [], gd);
end

function SetFileIndex(hObject, eventdata, gd)
if numel(hObject.String)>3
    hObject.String = hObject.String(end-2:end);
end
CreateFilename(gd.Saving.FullFilename, [], gd);
end

function toggleSave(hObject, eventdata, gd)
if hObject.Value
    hObject.BackgroundColor = [0,1,0];
    hObject.String = 'Saving';
else
    hObject.BackgroundColor = [1,0,0];
    hObject.String = 'Save?';
end
end

function CreateFilename(hObject, eventdata, gd)
gd.Experiment.saving.SaveFile = fullfile(gd.Internal.save.path, strcat(gd.Saving.base.String, '_', gd.Saving.depth.String, '_', gd.Saving.index.String, '.exp'));
hObject.String = gd.Experiment.saving.SaveFile;
guidata(hObject, gd);
if exist(gd.Experiment.saving.SaveFile, 'file')
    hObject.BackgroundColor = [1,0,0];
else
    hObject.BackgroundColor = [.94,.94,.94];
end
end


%% STIMULI CALLBACKS
% function linearMotorPos(hObject, eventdata, gd)
% if str2num(hObject.String) < 0
%     hObject.String = '0';
% elseif str2num(hObject.String) > 25
%     hObject.String = '25';
% end
% end
%
% function linearMotorMove(hObject, eventdata, gd)
% H_LinearStage = serial(gd.Internal.LinearStage.port, 'BaudRate', 9600);
% fopen(H_LinearStage);
% moveLinearMotor(str2num(gd.Controls.linearMotorPos.String), H_LinearStage); %move linear motor
% fclose(H_LinearStage);
% fprintf('Motor: Linear motor moved to position %.1f mm down track\n', str2num(gd.Controls.linearMotorPos.String));
% clear H_LinearStage;
% end

function EditPorts(hObject, eventdata, gd)
if eventdata.Indices(2)==4 % Trigger port
    if isempty(gd.Internal.daq)
        gd.Internal.daq = daq.createSession('ni');
        gd.Internal.daq.addDigitalChannel('Dev1', gd.Experiment.stim.setup.Port, 'OutputOnly');
    end
    gd.Internal.daq.outputSingleScan([hObject.Data{:,4}]);
else % Save user changes
    gd.Experiment.stim.setup = cell2table(hObject.Data(:,1:3),'VariableNames',{'Name','Port','Active'}); % save edits
    EditCombinations(gd.Stimuli.editCombinations, 'edit', gd);
end

guidata(hObject,gd); % update guidata
end

function EditStimuli(hObject, eventdata, gd)

% Delete stimulus
if eventdata.Indices(2)==2
    hObject.Data(eventdata.Indices(1),:) = [];
end

% Update registry
gd.Experiment.stim.pistonCombinations = cellfun(@str2num, hObject.Data(:,1), 'UniformOutput',false);
guidata(hObject, gd);

end

function EditCombinations(hObject, eventdata, gd)

if ischar(eventdata) % Update to all possible combinations
    hObject.String = num2str(1:nnz(gd.Experiment.stim.setup.Active));
else % Remove not possible combinations
    combinations = str2num(hObject.String);
    combinations(combinations > nnz(gd.Experiment.stim.setup.Active)) = [];
    combinations(combinations < 1) = [];
    hObject.String = num2str(combinations);
end

end

function GenerateCombinations(hObject, eventdata, gd)

% Generate combinations
activePorts = find(gd.Experiment.stim.setup.Active);
stimuli = [];
for cindex = str2num(gd.Stimuli.editCombinations.String)
    current = nchoosek(activePorts,cindex);
    stimuli = cat(1,stimuli,mat2cell(current, ones(size(current,1),1), size(current,2)));
end

% Display combinations
gd.Stimuli.list.Data = cellfun(@num2str, stimuli, 'UniformOutput',false);

% Update registry
gd.Experiment.stim.pistonCombinations = stimuli;

% Update number of trials
gd.Run.numTrials.String = num2str(size(stimuli,1)*gd.Internal.numStimuliRepetitions);

guidata(hObject, gd);

end

function LoadStimuli(hObject, eventdata, gd)

% Select and load file
[f,p] = uigetfile({'*.stim';'*.mat'},'Select stim file to load',cd);
if isnumeric(f)
    return
end
load(fullfile(p,f), 'stimuli', '-mat');
fprintf('Loaded stimuli from: %s\n', fullfile(p,f));

% Display combinations
gd.Stimuli.list.Data = cellfun(@num2str, stimuli, 'UniformOutput',false);

% Save loaded stimuli
gd.Experiment.stim.pistonCombinations = stimuli;

% Update number of trials
gd.Run.numTrials.String = num2str(size(stimuli,1)*10);

guidata(hObject, gd);
end

function SaveStimuli(hObject, eventdata, gd)
% Determine file to save to
[f,p] = uiputfile({'*.stim';'*.mat'},'Save stimuli as?',cd);
if isnumeric(f)
    return
end

% Save stimuli
stimuli = gd.Experiment.stim.pistonCombinations;
save(fullfile(p,f), 'stimuli', '-mat', '-v7.3');
fprintf('Stimuli saved to: %s\n', fullfile(p,f));
end


%% RUN EXPERIMENT
function RunExperiment(hObject, eventdata, gd)

if hObject.Value
    if isempty(gd.Experiment.stim.pistonCombinations)
        hObject.Value = 0;
        error('Load some stimulus combinations first');
    else
        try
            
            %% Record date & time information
            gd.Experiment.timing.init = datestr(now);
            
            %% Determine filenames to save to
            if gd.Saving.save.Value
                gd.Experiment.saving.save = true;
                % mat file
                if exist(gd.Experiment.saving.SaveFile, 'file')
                    answer = questdlg(sprintf('File already exists! Continue?\n%s', gd.Experiment.saving.SaveFile), 'Overwrite file?', 'Yes', 'No', 'No');
                    if strcmp(answer, 'No')
                        hObject.Value = false;
                        return
                    end
                end
                SaveFile = gd.Experiment.saving.SaveFile;
                % bin file
                gd.Experiment.saving.DataFile = strcat(gd.Experiment.saving.SaveFile(1:end-4), '.bin');
            else
                gd.Experiment.saving.save = false;
            end
            
            %% Initialize button
            hObject.BackgroundColor = [0,0,0];
            hObject.ForegroundColor = [1,1,1];
            hObject.String = 'Stop';
            
            %% Initialize NI-DAQ session
            gd.Internal.daq = [];
            
            DAQ = daq.createSession('ni'); % initialize session
            DAQ.IsContinuous = true; % set session to be continuous (call's 'DataRequired' listener)
            DAQ.Rate = gd.Experiment.params.samplingFrequency; % set sampling frequency
            gd.Experiment.params.samplingFrequency = DAQ.Rate; % the actual sampling frequency is rarely perfect from what is input
            
            % Add ports
            % Pistons
            activePorts = find(gd.Experiment.stim.setup.Active);
            for index = activePorts'
                [~,id] = DAQ.addDigitalChannel('Dev1',gd.Experiment.stim.setup.Port(index),'OutputOnly');
                DAQ.Channels(id).Name = strcat('O_Piston',index);
            end
            % Imaging Computer Trigger (for timing)
            [~,id] = DAQ.addDigitalChannel('Dev1','port0/line0','OutputOnly');
            DAQ.Channels(id).Name = 'O_2PTrigger';
            % Running Wheel
            [~,id] = DAQ.addDigitalChannel('Dev1','port0/line5:7','InputOnly');
            DAQ.Channels(id(1)).Name = 'I_RunWheelA';
            DAQ.Channels(id(2)).Name = 'I_RunWheelB';
            DAQ.Channels(id(3)).Name = 'I_RunWheelIndex';
%             % Whisker tracking
%             [~,id] = DAQ.addDigitalChannel('Dev1','port0/line1:2','OutputOnly');
%             DAQ.Channels(id(1)).Name = 'O_WhiskerTracker';
%             DAQ.Channels(id(2)).Name = 'O_WhiskerIllumination';
            [~,id] = DAQ.addDigitalChannel('Dev1','port0/line3','InputOnly');
            DAQ.Channels(id).Name = 'I_WhiskerTracker';
            % Cleanup
            DAQChannels = {DAQ.Channels(:).Name};
            OutChannels = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));
            InChannels = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
            
            % Add clock
            daqClock = daq.createSession('ni');
            daqClock.addCounterOutputChannel('Dev1',0,'PulseGeneration');
            clkTerminal = daqClock.Channels(1).Terminal;
            daqClock.Channels(1).Frequency = DAQ.Rate;
            daqClock.IsContinuous = true;
            daqClock.startBackground;
            DAQ.addClockConnection('External',['Dev1/' clkTerminal],'ScanClock');
            
            % Add QueueData callback
            DAQ.addlistener('DataRequired', @QueueData); % create listener for queueing trials
            DAQ.NotifyWhenScansQueuedBelow = DAQ.Rate-1; % queue more data when less than a second of data left
            % Add DataIn callback
            DAQ.addlistener('DataAvailable', @SaveDataIn);
            % DAQ.NotifyWhenDataAvailableExceeds = DAQ.Rate/100;
            
            %% Determine stimuli
            
            % Determine stimulus IDs
            gd.Experiment.StimID = 1:numel(gd.Experiment.stim.pistonCombinations);
            
            % Determine if presenting control stimulus
            gd.Experiment.stim.control = gd.Run.control.Value;
            if gd.Experiment.stim.control
                gd.Experiment.StimID = [0, gd.Experiment.StimID];
            end
            
            % Grab UI variables
            gd.Experiment.stim.blockShuffle = gd.Run.shuffle.Value;
            gd.Experiment.stim.repeatBadTrials = gd.Run.repeatBadTrials.Value;
            
            %% Create triggers
            
            % Compute timing of each trial
            gd.Experiment.timing.ITI = str2num(gd.Run.ITI.String);
            gd.Experiment.timing.stimDuration = str2num(gd.Run.stimDur.String);
            gd.Experiment.timing.trialDuration = gd.Experiment.timing.stimDuration + gd.Experiment.timing.ITI;
            gd.Experiment.timing.numScansPerTrial = ceil(gd.Experiment.params.samplingFrequency * gd.Experiment.timing.trialDuration);
            if gd.Experiment.timing.numScansPerTrial < DAQ.NotifyWhenScansQueuedBelow
                DAQ.NotifyWhenScansQueuedBelow = gd.Experiment.timing.numScansPerTrial - 1;
            end
            
            % Initialize blank triggers
            gd.Experiment.blankTriggers = zeros(gd.Experiment.timing.numScansPerTrial, numel(OutChannels));
            
            % Trigger pistons
            startTrig = max(floor(gd.Experiment.params.samplingFrequency * gd.Experiment.timing.ITI),1);    % start after ITI
            endTrig = gd.Experiment.timing.numScansPerTrial-1;                                              % end on last trigger of trial
            gd.Experiment.PistonTrigger = zeros(gd.Experiment.timing.numScansPerTrial, 1);
            gd.Experiment.PistonTrigger(startTrig:endTrig) = 1;
            
            % Trigger imaging computer on every single trial
            gd.Experiment.blankTriggers([startTrig, endTrig], strcmp(OutChannels, 'O_2PTrigger')) = 1; % trigger at beginning and end of stimulus
            
            % Trigger whisker tracking camera on every single trial
            if gd.Experiment.timing.ITI >= 0.01;
                gd.Experiment.blankTriggers(startTrig-ceil(DAQ.Rate/100):endTrig, strcmp(OutChannels, 'O_WhiskerIllumination')) = 1; % start LED a little before the start of imaging
                gd.Experiment.blankTriggers(startTrig:ceil(DAQ.Rate/gd.Experiment.params.frameRateWT):endTrig, strcmp(OutChannels, 'O_WhiskerTracker')) = 1; % image during stimulus period
            else
                gd.Experiment.blankTriggers(:, strcmp(OutChannels, 'O_WhiskerIllumination')) = 1; % image during entire time
                gd.Experiment.blankTriggers(1:ceil(DAQ.Rate/gd.Experiment.params.frameRateWT):endTrig, strcmp(OutChannels, 'O_WhiskerTracker')) = 1; % image during entire time
            end
            % gd.Experiment.blankTriggers(1, strcmp(OutChannels, 'O_WhiskerTracker')) = 1;        % mode 15 limited to 255 frames
            % gd.Experiment.blankTriggers(stopMove1, strcmp(OutChannels, 'O_WhiskerTracker')) = 1;
            
            % Build up vector to display when stimulus is present
            gd.Experiment.Stimulus = zeros(size(gd.Experiment.blankTriggers,1), 1);
            gd.Experiment.Stimulus(startTrig:endTrig) = 1;
            
            
            %% Initialize imaging session (scanbox only)
            ImagingType = gd.Saving.imagingType.String{gd.Saving.imagingType.Value};
            if strcmp(ImagingType, 'sbx')
                H_Scanbox = udp(gd.Internal.ImagingComp.ip, 'RemotePort', gd.Internal.ImagingComp.port); % create udp port handle
                fopen(H_Scanbox);
                fprintf(H_Scanbox,sprintf('A%s',gd.Saving.base.String));
                fprintf(H_Scanbox,sprintf('U%s',gd.Saving.depth.String));
                fprintf(H_Scanbox,sprintf('E%s',gd.Saving.index.String));
            end
            
            
            %% Initialize saving
            Experiment = gd.Experiment;
            if gd.Experiment.saving.save
                save(SaveFile, 'DAQChannels', 'Experiment', '-mat', '-v7.3');
                H_DataFile = fopen(gd.Experiment.saving.DataFile, 'w');
            end
            
            
            %% Initialize shared variables (only share what's necessary)
            
            % Necessary variables
            numTrialsObj = gd.Run.numTrials;
            numStimuli = numel(Experiment.StimID);
            numStimuliCurrentBlock = numStimuli;
            BaseTriggers = Experiment.blankTriggers;
            PistonTrigger = Experiment.PistonTrigger;
            ControlTrial = Experiment.stim.control;
            BlockShuffle = Experiment.stim.blockShuffle;
            currentBlockOrder = Experiment.StimID;
            currentTrial = 0;
            PistonCombinations = Experiment.stim.pistonCombinations;
            TrialInfo = struct('StimID', [], 'Running', [], 'RunSpeed', []);
            saveOut = Experiment.saving.save;
            Stimulus = Experiment.Stimulus;
            ExperimentReachedEnd = false; % boolean to see if max trials has been reached
            numScansPerTrial = Experiment.timing.numScansPerTrial;
            
            % Variables if saving input data
            if saveOut
                Precision = Experiment.saving.dataPrecision;
            end
            
            % Variables for calculating and displaying running speed
            RunChannelIndices = [find(strcmp(InChannels, 'I_RunWheelB')),find(strcmp(InChannels,'I_RunWheelA'))];
            numBufferScans = gd.Internal.buffer.numTrials*numScansPerTrial;
            DataInBuffer = zeros(numBufferScans, 2);
            dsamp = gd.Internal.buffer.downSample;
            dsamp_Fs = Experiment.params.samplingFrequency / dsamp;
            smooth_win = gausswin(dsamp_Fs, 23.5/2);
            smooth_win = smooth_win/sum(smooth_win);
            sw_len = length(smooth_win);
            d_smooth_win = [0;diff(smooth_win)]/(1/dsamp_Fs);
            hAxes = gd.Run.runSpeedAxes;
            
            % Variables for displaying stim info
            numScansReturned = DAQ.NotifyWhenDataAvailableExceeds;
            BufferStim = zeros(numBufferScans, 1);
                        
            % Variables for determing if mouse was running
            RepeatBadTrials = Experiment.stim.repeatBadTrials;
            StimuliToRepeat = [];
            SpeedThreshold = Experiment.stim.speedThreshold;
            RunIndex = 1;
            
            %% Start Experiment
            % Start imaging
            if strcmp(ImagingType, 'sbx')
                fprintf(H_Scanbox,'G'); %go
                pause(5);
            end
            
            % Start experiment
            Experiment.timing.start = datestr(now);
            QueueData();
            DAQ.startBackground;
            
            
            %% During Experiment
            while DAQ.IsRunning
                pause(0.1);
            end
            Experiment.timing.finish = datestr(now);
            
            %% End Experiment
            if strcmp(ImagingType, 'sbx')
                fprintf(H_Scanbox,'S'); %stop
                fclose(H_Scanbox);
            end
            if saveOut
                save(SaveFile, 'Experiment', '-append'); % update with "Experiment.timing.finish" info
                fclose(H_DataFile);                      % close binary file
                gd.Saving.index.String = sprintf('%03d',str2double(gd.Saving.index.String) + 1); % update file index for next experiment
                CreateFilename(gd.Saving.FullFilename, [], gd); % update filename for next experiment
            end
            
            % Reset button properties
            hObject.Value = false;
            hObject.BackgroundColor = [.94,.94,.94];
            hObject.ForegroundColor = [0,0,0];
            hObject.String = 'Run';
            
        catch ME
            warning('Running experiment failed');
            
            % Reset button properties
            hObject.Value = false;
            hObject.BackgroundColor = [.94,.94,.94];
            hObject.ForegroundColor = [0,0,0];
            hObject.String = 'Run';
            
            % Close any open connections
            try
                if strcmp(ImagingType, 'sbx')
                    fclose(H_Scanbox);
                end
            end
            clear DAQ H_Scanbox
            
            % Rethrow error
            rethrow(ME);
        end
    end  
    
else % user quit experiment (hObject.Value = false)
    
    % Change button properties to reflect change in state
    hObject.BackgroundColor = [1,1,1];
    hObject.ForegroundColor = [0,0,0];
    hObject.String = 'Stopping...';
    
end

%% Callback: DataIn
    function SaveDataIn(src,eventdata)
        
        % Save input data
        if saveOut
            fwrite(H_DataFile, eventdata.Data', Precision);
        end
        
        % Refresh buffer
        DataInBuffer = cat(1, DataInBuffer(numScansReturned+1:end,:), eventdata.Data(:,RunChannelIndices));   % concatenate new data and remove old data
        BufferStim = BufferStim(numScansReturned+1:end);
        
        % Convert entire buffer of pulses to run speed
        Data = [0;diff(DataInBuffer(:,1))>0];       % gather pulses' front edges
        Data(all([Data,DataInBuffer(:,2)],2)) = -1; % set backwards steps to be backwards
        Data = cumsum(Data);                        % convert pulses to counter data
        x_t = downsample(Data, dsamp);              % downsample data to speed up computation
        x_t = padarray(x_t, sw_len, 'replicate');   % pad for convolution
        dx_dt = conv(x_t, d_smooth_win, 'same');    % perform convolution
        dx_dt([1:sw_len,end-sw_len+1:end]) = [];    % remove values produced by padding the data
        % dx_dt = dx_dt * 360/360;                  % convert to degrees (360 pulses per 360 degrees)
        
        % Record average running speed during stimulus period
        currentStim = downsample(BufferStim(1:numBufferScans), dsamp);
        if any(diff(BufferStim(numBufferScans-numScansReturned:numBufferScans)) == -RunIndex) % stimulus ended during current DataIn call
            % max(BufferStim(numBufferScans-numScansReturned+1:numBufferScans)) > RunIndex || diff(BufferStim([numBufferScans-numScansReturned,numBufferScans])) == -RunIndex     % next stimulus is already running
            TrialInfo.RunSpeed(RunIndex) = mean(dx_dt(currentStim==RunIndex));  % calculate average running speed during stimulus
            fprintf('\t\t\tT%d S%d RunSpeed= %.2f', RunIndex, TrialInfo.StimID(RunIndex), TrialInfo.RunSpeed(RunIndex));
            
            % Save to file
            if saveOut
                save(SaveFile, 'TrialInfo', '-append');
            end
            
            % Determine if mouse was running during previous trial
            if TrialInfo.RunSpeed(RunIndex) < SpeedThreshold
                TrialInfo.Running(RunIndex) = false;                                        % record that trial was bad
                if RepeatBadTrials                                                          % repeate trial
                    fprintf(' (trial to be repeated)');
                    StimuliToRepeat = [StimuliToRepeat, TrialInfo.StimID(RunIndex)];        % add trial to repeat queue
                end
            else % mouse was running
                TrialInfo.Running(RunIndex) = true;                                         % record that trial was good
            end
            RunIndex = RunIndex+1; % increment index
            
        end
        
        % Display new data and stimulus info
        plot(hAxes, numBufferScans:-dsamp:1, dx_dt, 'b-', numBufferScans:-dsamp:1, 400*(currentStim>0), 'r-');
        ylim([0,500]);
        
    end %SaveDateIn

%% Callback: QueueOutputData
    function QueueData(src,eventdata)
        
        % Queue trial
        if currentTrial < str2double(numTrialsObj.String) && hObject.Value
            
            % Update index
            currentTrial = currentTrial + 1;
            ExperimentReachedEnd = false;
        
            blockIndex = rem(currentTrial-1, numStimuliCurrentBlock)+1;
            % If starting new block, shuffle the stimuli order
            if BlockShuffle && blockIndex == 1
                numStimuliCurrentBlock = numStimuli;                            % reset size of block (changes if repeating bad trials)
                currentBlockOrder = currentBlockOrder(randperm(numStimuli));    % shuffle block
            end
            
            % Queue triggers
            TrialInfo.StimID(currentTrial) = currentBlockOrder(blockIndex);     % determine StimID for current trial
            if ControlTrial && TrialInfo.StimID(currentTrial) == 0              % current trial is control trial
                DAQ.queueOutputData(BaseTriggers);
            else                                                                % current trial is not control trial
                CurrentTriggers = BaseTriggers;
                CurrentTriggers(:,PistonCombinations{TrialInfo.StimID(currentTrial)+ControlTrial}) = repmat(PistonTrigger, 1, numel(PistonCombinations{TrialInfo.StimID(currentTrial)}+ControlTrial));
                DAQ.queueOutputData(CurrentTriggers);
            end
            BufferStim = cat(1, BufferStim, Stimulus*currentTrial);
            
            % Update information
            fprintf('\nQueued trial %d: stimulus %d', currentTrial, TrialInfo.StimID(currentTrial));
            if saveOut
                save(SaveFile, 'TrialInfo', '-append');
            end
            
        elseif ~isempty(StimuliToRepeat) % Repeat bad trials
            
            % Update index
            currentTrial = currentTrial + 1;
            ExperimentReachedEnd = false;
            
            % Queue triggers
            TrialInfo.StimID(currentTrial) = StimuliToRepeat(1);                % determine StimID for current trial
            StimuliToRepeat(1) = [];
            if ControlTrial && TrialInfo.StimID(currentTrial) == 0              % current trial is control trial
                DAQ.queueOutputData(BaseTriggers);
            else                                                                % current trial is not control trial
                CurrentTriggers = BaseTriggers;
                CurrentTriggers(:,PistonCombinations{TrialInfo.StimID(currentTrial)}) = repmat(PistonTrigger, 1, numel(PistonCombinations{TrialInfo.StimID(currentTrial)}));
                DAQ.queueOutputData(CurrentTriggers);
            end
            BufferStim = cat(1, BufferStim, Stimulus*currentTrial);
            
            % Update information
            fprintf('\nQueued trial %d: stimulus %d', currentTrial, TrialInfo.StimID(currentTrial));
            if saveOut
                save(SaveFile, 'TrialInfo', '-append');
            end
            
        elseif currentTrial >= str2double(numTrialsObj.String) && ~ExperimentReachedEnd
            % Queue blank trial to ensure last trial does not have to be
            % repeated and to ensure no imaging frames get clipped
            DAQ.queueOutputData(zeros(numScansPerTrial, numel(OutChannels)));
            BufferStim = cat(1, BufferStim, zeros(numScansPerTrial, 1));
            ExperimentReachedEnd = true;
                
        elseif ~hObject.Value
                fprintf('\nComplete: finished %d trial(s) (user quit)\n', currentTrial);
            
        elseif currentTrial >= str2double(numTrialsObj.String)
            fprintf('\nComplete: finished %d trials(s) (max trials reached)\n', currentTrial);
            
            
        end
    end

end
function [TrialInfo, SaveFile] = stimNubulate(SaveFile, varargin)


%% Configure Settings
gd.Internal.ImagingType = 'sbx';                % 'sbx' or 'scim'
gd.Internal.ImagingComp.ip = '128.32.173.30';   % SCANBOX ONLY: for UDP
gd.Internal.ImagingComp.port = 7000;            % SCANBOX ONLY: for UDP
% gd.Internal.LinearStage.port = 'com3';

Display.units = 'pixels';
Display.position = [100, 100, 800, 400];

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

gd.Experiment.timing.trialDuration = 6; % in seconds
gd.Experiment.timing.preDuration = 1.5; % in seconds
gd.Experiment.timing.stimDuration = 1.5; % in seconds

gd.Experiment.stim.setup = table(...
    {'C1';'C2';'B1';'D1';'beta';'gamma';'';''},...
    {'port0/line8';'port0/line9';'port0/line10';'port0/line11';'port0/line12';'port0/line13';'port0/line14';'port0/line15'},...
    [true;true;true;true;true;true;false;false],...
    'VariableNames',{'Name','Port','Active'});
gd.Experiment.stim.pistonCombinations = {};
gd.Experiment.stim.control = true; % true or false (give control stimulus)

% Properties for display or processing input data
gd.Internal.buffer.numscans = gd.Experiment.params.samplingFrequency * 6;
gd.Internal.buffer.downSample = 20;

% Place holder for offline control of stimulus
gd.Internal.daq = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'trialDuration'
                gd.Experiment.timing.trialDuration = varargin{index+1};
                index = index + 2;
            case 'preDuration'
                gd.Experiment.timing.preDuration = varargin{index+1};
                index = index + 2;
            case 'stimDuration'
                gd.Experiment.timing.stimDuration = varargin{index+1};
                index = index + 2;
            case 'numPositions'
                gd.Experiment.stim.numPositions = varargin{index+1};
                index = index + 2;
            case 'distBetween'
                gd.Experiment.stim.distBetweenPositions = varargin{index+1};
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
    'Name',                 'Stimulus: Vertical Pole along Azimuth',...
    'ToolBar',              'none',...
    'Units',                Display.units,...
    'Position',             Display.position,...
    'KeyPressFcn',          @(hObject,eventdata)KeyPressCallback(hObject, eventdata, guidata(hObject)));

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
% number of trials
gd.Run.numTrials = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.numTrials,...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.9,.5,.1]);
gd.Run.numTrialsText = uicontrol(...
    'Style',                'text',...
    'String',               '# Trials',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlignment',  'right',...
    'Units',                'normalized',...
    'Position',             [0,.925,.5,.05]);
% control toggle
gd.Run.control = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Control Trial?',...
    'Parent',               gd.Run.panel,...
    'Value',                gd.Experiment.stim.control,...
    'Units',                'normalized',...
    'Position',             [.5,.8,.5,.1]);
% run button
gd.Run.run = uicontrol(...
    'Style',                'togglebutton',...
    'String',               'Run?',...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [0,.6,1,.2],...
    'Callback',             @(hObject,eventdata)RunExperiment(hObject, eventdata, guidata(hObject)));
% running speed axes
gd.Run.runSpeedAxes = axes(...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,1,.6]);

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
gd.Experiment.stim.pistonCombinations = cellfun(@str2num, hObject.Data(:,1:3), 'UniformOutput',false);

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
gd.Run.numTrials.String = num2str(size(stimuli,1)*10);

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
        % Whisker tracking
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line1:2','OutputOnly');
        DAQ.Channels(id(1)).Name = 'O_WhiskerTracker';
        DAQ.Channels(id(2)).Name = 'O_WhiskerIllumination';
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line3','InputOnly');
        DAQ.Channels(id).Name = 'I_WhiskerTracker';
        % Cleanup
        DAQChannels = {DAQ.Channels(:).Name};
        OutChannels = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));
        InChannels = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
        
        % Add clock
        daqClock = daq.createSession('ni');
        daqClock.addCounterOutputChannel('Dev1',0,'PulseGeneration')
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
        
        
        %% Create triggers
        
        % Determine stimulus IDs
        gd.Experiment.StimID = 1:numel(gd.Experiment.stim.pistonCombinations);

        % Determine if presenting control stimulus
        gd.Experiment.stim.control = gd.Run.control.Value;
        if gd.Experiment.stim.control
            gd.Experiment.StimID = [0, gd.Experiment.StimID];
        end
        
        % Initial
        gd.Experiment.Triggers = zeros(round(gd.Experiment.params.samplingFrequency*gd.Experiment.timing.trialDuration), numel(OutChannels));
        
        % Trigger pistons
        startTrig = floor(gd.Experiment.params.samplingFrequency*gd.Experiment.timing.preDuration);
        endTrig = ceil(gd.Experiment.params.samplingFrequency*(gd.Experiment.timing.preDuration+gd.Experiment.timing.stimDuration));
        gd.Experiment.PistonTrigger = zeros(size(gd.Experiment.Triggers,1), 1);
        gd.Experiment.PistonTrigger(startTrig:endTrig) = 1;
        
        % Trigger imaging computer
        gd.Experiment.Triggers([startTrig, endTrig], strcmp(OutChannels, 'O_2PTrigger')) = 1;
        
        % Trigger whisker tracking camera
        gd.Experiment.Triggers(startTrig-ceil(DAQ.Rate/100):endTrig, strcmp(OutChannels, 'O_WhiskerIllumination')) = 1;
        gd.Experiment.Triggers(startTrig:ceil(DAQ.Rate/gd.Experiment.params.frameRateWT):endTrig, strcmp(OutChannels, 'O_WhiskerTracker')) = 1;
        % gd.Experiment.Triggers(1, strcmp(OutChannels, 'O_WhiskerTracker')) = 1;        % mode 15 limited to 255 frames
        % gd.Experiment.Triggers(stopMove1, strcmp(OutChannels, 'O_WhiskerTracker')) = 1;
        
        % Stimulus info
        gd.Experiment.Stimulus = zeros(size(gd.Experiment.Triggers,1), 1);
        gd.Experiment.Stimulus(startTrig:endTrig) = 1;
        
        
%         %% Initialize linear motor
%         H_LinearStage = serial(gd.Internal.LinearStage.port, 'BaudRate', 9600);
%         fopen(H_LinearStage);
        

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
        numTrials = str2num(gd.Run.numTrials.String);
        numStimuli = numel(gd.Experiment.StimID);
        BaseTriggers = gd.Experiment.Triggers;
        PistonTrigger = gd.Experiment.PistonTrigger;
        ControlTrial = gd.Experiment.stim.control;
        currentBlockOrder = gd.Experiment.StimID;
        currentTrial = 0;
        PistonCombinations = gd.Experiment.stim.pistonCombinations;
        TrialInfo = struct('StimID', []);
        saveOut = gd.Experiment.saving.save;

        % Variables if saving input data
        if saveOut
            Precision = gd.Experiment.saving.dataPrecision;
        end
        
        % Variables for calculating and displaying running speed
        DataInBuffer = zeros(gd.Internal.buffer.numscans, 1);
        numBufferScans = gd.Internal.buffer.numscans;
        dsamp_Fs = gd.Experiment.params.samplingFrequency / gd.Internal.buffer.downSample;
        smooth_win = gausswin(dsamp_Fs, 23.5/2);
        smooth_win = smooth_win/sum(smooth_win);
        sw_len = length(smooth_win);
        d_smooth_win = [0;diff(smooth_win)]/(1/dsamp_Fs);
        hAxes = gd.Run.runSpeedAxes;
        
        % Variables for displaying stim info
        numScansPerTrial = size(Experiment.Triggers, 1);
        scanCount = 0;
        numScansReturned = DAQ.NotifyWhenDataAvailableExceeds;
        Stimulus = 400*repmat(Experiment.Stimulus, 2, 1);
        
        
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
            numTrials = str2num(gd.Run.numTrials.String); % update in case request changes
        end
        Experiment.timing.finish = datestr(now);
        
        %% End Experiment
%         fclose(H_LinearStage);
        if strcmp(ImagingType, 'sbx')
            fprintf(H_Scanbox,'S'); %stop
            fclose(H_Scanbox);
        end
        if saveOut
            save(SaveFile, 'Experiment', '-append');
            fclose(H_DataFile);
            gd.Saving.index.String = sprintf('%03d',str2num(gd.Saving.index.String) + 1);
            CreateFilename(gd.Saving.FullFilename, [], gd);
        end

        hObject.Value = false;
        hObject.BackgroundColor = [.94,.94,.94];
        hObject.ForegroundColor = [0,0,0];
        hObject.String = 'Run';      
        
    catch ME
        warning('Running experiment failed');
        hObject.Value = false;
        hObject.BackgroundColor = [.94,.94,.94];
        hObject.ForegroundColor = [0,0,0];
        hObject.String = 'Run';
        try
            fclose(H_LinearStage);
        end
        try
            if strcmp(ImagingType, 'sbx')
                fclose(H_Scanbox);
            end
        end
        clear DAQ H_LinearStage H_Scanbox
        rethrow(ME);
    end
    
    
else % hObject.Value == false -> user quit experiment
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
        
        % Update buffer
%         temp = eventdata.Data;
        DataInBuffer = cat(1, DataInBuffer, eventdata.Data(:,strcmp(InChannels, 'I_RunWheelA')));
        DataInBuffer = DataInBuffer(end-numBufferScans+1:end,:);
        
        % Update stim location
        scanCount = scanCount + numScansReturned;
        currentScan = rem(scanCount, numScansPerTrial) + numScansPerTrial;
        
        % Determine running speed
        Data = cumsum(diff(DataInBuffer)>0);
        x_t = downsample(Data, 20); %downsample data to speed up computation
        x_t = padarray(x_t, sw_len, 'replicate');
        dx_dt = conv(x_t, d_smooth_win, 'same');
        dx_dt([1:sw_len,end-sw_len+1:end]) = []; %remove values produced by convolving kernel with padded values
        % RunningSpeed = dx_dt*360/360; % convert to degrees (360 pulses per 360 degrees)
       
        % Display input data
        plot(hAxes, numBufferScans:-20:1, dx_dt, 'b-', numBufferScans:-1:1, Stimulus(currentScan-numBufferScans+1:currentScan), 'r-');
        ylim([0,500]);
        
%         % Display stimulus data
%         numCurrentScans = size(eventdata.Data, 1);
%         numScans = numScans + numCurrentScans;
%         scanIndex = rem(numScans-1, numScansPerTrial)+1;
%         hold on;
%         plot(hAxes, 250*Stimulus(scanIndex - numBufferScans + 1 : scanIndex), 'r-o');
%         hold off;
    end

%% Callback: QueueOutputData
    function QueueData(src,eventdata)
        
        % Update indices
        currentTrial = currentTrial + 1;
        blockIndex = rem(currentTrial-1, numStimuli)+1;
        
        % Queue trial
        if currentTrial <= numTrials && hObject.Value
            
            % If starting new block, shuffle the stimuli order
            if blockIndex == 1
                currentBlockOrder = currentBlockOrder(randperm(numStimuli));                    % shuffle block
            end
            
            % Determine stimulus for current trial
            TrialInfo.StimID(currentTrial) = currentBlockOrder(blockIndex);                     % determine StimID for current trial
            
            % Queue output data
            if ControlTrial && TrialInfo.StimID(currentTrial) == 0                              % current trial is control trial
                DAQ.queueOutputData(BaseTriggers);
            else                                                                                % current trial is not control trial
                CurrentTriggers = BaseTriggers;
                CurrentTriggers(:,PistonCombinations{TrialInfo.StimID(currentTrial)}) = repmat(PistonTrigger, 1, numel(PistonCombinations{TrialInfo.StimID(currentTrial)}));
                DAQ.queueOutputData(CurrentTriggers);
            end
            
            % Update information
            fprintf('\nQueued trial %d: stimulus %d', currentTrial, TrialInfo.StimID(currentTrial));
            if saveOut
                save(SaveFile, 'TrialInfo', '-append');
            end
            
        elseif ~hObject.Value
            fprintf('\nComplete: finished %d trial(s) (user quit)\n', currentTrial-1);
            
        elseif currentTrial > numTrials
            fprintf('\nComplete: finished %d trials(s) (max trials reached)\n', currentTrial-1);
        end
    end

end
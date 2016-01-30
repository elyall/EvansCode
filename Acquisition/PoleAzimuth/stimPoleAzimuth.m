function [TrialInfo, SaveFile] = stimPoleAzimuth(SaveFile, varargin)


%% Configure Settings
gd.Internal.ImagingType = 'sbx'; % 'sbx' or 'scim'
gd.Internal.ImagingComp.ip = '128.32.173.30'; % SCANBOX ONLY: for UDP
gd.Internal.ImagingComp.port = 7000; % SCANBOX ONLY: for UDP
gd.Internal.LinearStage.APport = 'com3';
gd.Internal.LinearStage.MLport = 'com4';

Display.units = 'pixels';
Display.position = [100, 100, 800, 400];

%% Initialize Experiment

gd.Internal.save.path = 'C:\Users\Resonant-2\OneDrive\StimData';
gd.Internal.save.base = '0000';
gd.Internal.save.depth = '000';
gd.Internal.save.index = '000';

gd.Experiment.saving.save = false;
gd.Experiment.saving.SaveFile = fullfile(gd.Internal.save.path, gd.Internal.save.base);
gd.Experiment.saving.DataFile = '';
gd.Experiment.saving.dataPrecision = 'uint16';

gd.Experiment.params.samplingFrequency = 30000;
gd.Experiment.params.numTrials = 180;
gd.Experiment.params.angleMove = 30;
gd.Experiment.params.frameRateWT = 125;

gd.Experiment.timing.trialDuration = 6; % in seconds
gd.Experiment.timing.preDuration = 1.5; % in seconds
gd.Experiment.timing.stimDuration = 1.5; % in seconds

gd.Experiment.stim.numPositions = 8; % number of stimulus positions
gd.Experiment.stim.distBetweenPositions = 2.5; % in mm
gd.Experiment.stim.control = true; % true or false (give control stimulus)

gd.Internal.buffer.numscans = gd.Experiment.params.samplingFrequency * 6;
gd.Internal.buffer.downSample = 20;

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

% CONTROLS
% panel
gd.Controls.panel = uipanel(...
    'Title',                'Controls',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [0, 0, .5, .8]);
% stepper motor CCW
gd.Controls.stepCCW = uicontrol(...
    'Style',                'pushbutton',...
    'String',               '<',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [0,.8,.35,.1],...
    'Callback',             @(hObject,eventdata)stepCCW(hObject, eventdata, guidata(hObject)));
% stepper motor degrees
gd.Controls.stepAngle = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.angleMove,...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.35,.8,.3,.1]);
% stepper motor CW
gd.Controls.stepCW = uicontrol(...
    'Style',                'pushbutton',...
    'String',               '>',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.65,.8,.35,.1],...
    'Callback',             @(hObject,eventdata)stepCW(hObject, eventdata, guidata(hObject)));
% stepper motor text
gd.Controls.stepperMotorText = uicontrol(...
    'Style',                'text',...
    'String',               'Move stepper motor (relative)',...
    'Parent',               gd.Controls.panel,...
    'HorizontalAlignment',  'center',...
    'Units',                'normalized',...
    'Position',             [0,.925,1,.05]);
% anterior-posterior linear motor position
gd.Controls.apLinearMotorPos = uicontrol(...
    'Style',                'edit',...
    'String',               '0',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.35,.6,.3,.1],...
    'Callback',             @(hObject,eventdata)linearMotorPos(hObject, eventdata, guidata(hObject)));
% ap linear motor move
gd.Controls.apMotorMove = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'move',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'UserData',             'ap',...
    'Position',             [.65,.6,.35,.1],...
    'Callback',             @(hObject,eventdata)linearMotorMove(hObject, eventdata, guidata(hObject)));
% ap linear motor text
gd.Controls.mlMotorText = uicontrol(...
    'Style',                'text',...
    'String',               'Move anterior-posterior motor (absolute)',...
    'Parent',               gd.Controls.panel,...
    'HorizontalAlignment',  'center',...
    'Units',                'normalized',...
    'Position',             [0,.725,1,.05]);
% medial-lateral linear motor position
gd.Controls.mlLinearMotorPos = uicontrol(...
    'Style',                'edit',...
    'String',               '0',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'UserData',             'ml',...
    'Position',             [.35,.4,.3,.1],...
    'Callback',             @(hObject,eventdata)linearMotorPos(hObject, eventdata, guidata(hObject)));
% ml linear motor move
gd.Controls.mlMotorMove = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'move',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.65,.4,.35,.1],...
    'Callback',             @(hObject,eventdata)linearMotorMove(hObject, eventdata, guidata(hObject)));
% ml linear motor toggle
gd.Controls.mlMotorToggle = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Activate?',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.01,.4,.34,.1],...
    'BackgroundColor',      [1,0,0],...
    'UserData',             {[1,0,0;0,1,0],'Activate?','Active'},...
    'Callback',             @(hObject,eventdata)set(hObject,'BackgroundColor',hObject.UserData{1}(hObject.Value+1,:),'String',hObject.UserData{hObject.Value+2}));
% ml linear motor text
gd.Controls.linearMotorText = uicontrol(...
    'Style',                'text',...
    'String',               'Move medial-lateral motor (absolute)',...
    'Parent',               gd.Controls.panel,...
    'HorizontalAlignment',  'center',...
    'Units',                'normalized',...
    'Position',             [0,.525,1,.05]);
% position table
gd.Controls.positionTable = uitable(...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.01,.1,1,.25],...
    'RowName',              {'AP','ML'},...
    'ColumnName',           {'First Pos (mm)','Last Pos (mm)','# Steps'},...
    'ColumnEditable',       [true, true, true],...
    'ColumnFormat',         {'numeric','numeric','numeric'},...
    'ColumnWidth',          {100,100,100},...
    'Data',                 [zeros(2,1),[20;5].*ones(2,1),nan(2,1)]);

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
% angle
gd.Run.angle = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.angleMove,...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.8,.5,.1]);
gd.Parameters.numTrialsText = uicontrol(...
    'Style',                'text',...
    'String',               'Angle to Move',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlignment',  'right',...
    'Units',                'normalized',...
    'Position',             [0,.825,.5,.05]);
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

%% Callbacks
function ImagingSelection(hObject, eventdata, gd)
if isequal(gd.Saving.tabgroup.SelectedTab, gd.Saving.scim.tab)
    gd.Internal.ImagingType = 'scim';
elseif isequal(gd.Saving.tabgroup.SelectedTab, gd.Saving.sbx.tab)
    gd.Internal.ImagingType = 'sbx';
end
guidata(hObject, gd);
end

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

%% CONTORLS CALLBACKS
function linearMotorPos(hObject, eventdata, gd)
if str2num(hObject.String) < 0
    hObject.String = '0';
elseif str2num(hObject.String) > 25
    hObject.String = '25';
end
end

function linearMotorMove(hObject, eventdata, gd)
switch hObject.UserData
    case 'ap'
        H_LinearStage = serial(gd.Internal.LinearStage.APport, 'BaudRate', 9600);
        fopen(H_LinearStage);
        moveLinearMotor(str2num(gd.Controls.apMotorPos.String), H_LinearStage); %move linear motor
        fclose(H_LinearStage);
        fprintf('Anterior-Posterior motor: moved to position %.1f mm down track\n', str2num(gd.Controls.linearMotorPos.String));
        clear H_LinearStage;
    case 'ml'
        H_LinearStage = serial(gd.Internal.LinearStage.MLport, 'BaudRate', 9600);
        fopen(H_LinearStage);
        moveLinearMotor(str2num(gd.Controls.mlMotorPos.String), H_LinearStage); %move linear motor
        fclose(H_LinearStage);
        fprintf('Medial-Lateral motor: Linear motor moved to position %.1f mm down track\n', str2num(gd.Controls.linearMotorPos.String));
        clear H_LinearStage;
end
end

function stepCCW(hObject, eventdata, gd)
moveStepperMotor(str2num(gd.Controls.stepAngle.String), gd.Experiment.params.samplingFrequency, 1, 2, true, 2);
fprintf('Motor: Stepper motor moved %d degrees CCW\n', str2num(gd.Controls.stepAngle.String));
end

function stepCW(hObject, eventdata, gd)
moveStepperMotor(-str2num(gd.Controls.stepAngle.String), gd.Experiment.params.samplingFrequency, 1, 2, true, 2);
fprintf('Motor: Stepper motor moved %d degrees CW\n', str2num(gd.Controls.stepAngle.String));
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
        DAQ = daq.createSession('ni'); % initialize session
        DAQ.IsContinuous = true; % set session to be continuous (call's 'DataRequired' listener)
        DAQ.Rate = gd.Experiment.params.samplingFrequency; % set sampling frequency
        gd.Experiment.params.samplingFrequency = DAQ.Rate; % the actual sampling frequency is rarely perfect from what is input
        
        % Add clock
        % daqClock = daq.createSession('ni');
        % daqClock.addCounterOutputChannel('Dev1',0,'PulseGeneration')
        % clkTerminal = daqClock.Channels(1).Terminal;
        % daqClock.Channels(1).Frequency = DAQ.Rate;
        % daqClock.IsContinuous = true;
        % daqClock.startBackground;
        % DAQ.addClockConnection('External',['Dev1/' clkTerminal],'ScanClock');
        
        % Add ports
        % Imaging Computer Trigger (for timing)
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line0','OutputOnly');
        DAQ.Channels(id).Name = 'O_2PTrigger';
        % Motor Rotation
        [~,id] = DAQ.addAnalogOutputChannel('Dev1',0:1,'Voltage');
        % [~,id] = DAQ.addDigitalChannel('Dev1','port0/line4:5','OutputOnly');
        DAQ.Channels(id(1)).Name = 'O_MotorStep';
        DAQ.Channels(id(2)).Name = 'O_MotorDir';
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
        
        % Add QueueData callback
        DAQ.addlistener('DataRequired', @QueueData); % create listener for queueing trials
        DAQ.NotifyWhenScansQueuedBelow = DAQ.Rate-1; % queue more data when less than a second of data left
        % Add DataIn callback
        DAQ.addlistener('DataAvailable', @SaveDataIn);
        % DAQ.NotifyWhenDataAvailableExceeds = DAQ.Rate/100;
        
        
        %% Create triggers

        
        % Gather GUI data
        gd.Experiment.params.angleMove = str2num(gd.Run.angle.String);
        
        % Initial
        gd.Experiment.Triggers = zeros(round(gd.Experiment.params.samplingFrequency*gd.Experiment.timing.trialDuration), numel(OutChannels), 1+gd.Experiment.stim.control);
        % gd.Experiment.params.frameRateWT
        % gd.Experiment.stim.numPositions
        % gd.Experiment.stim.distBetweenPositions
        
        % Stepper Motor triggers
        stepTriggers = moveStepperMotor(gd.Experiment.params.angleMove, gd.Experiment.params.samplingFrequency, 1, 2, false, 2);
        numStepTriggers = size(stepTriggers, 1);
        
        % Move In Stepper Motor
        stopMove1 = round(gd.Experiment.timing.preDuration * gd.Experiment.params.samplingFrequency);
        beginMove1 = stopMove1 - numStepTriggers + 1; 
        gd.Experiment.Triggers(beginMove1:stopMove1, strcmp(OutChannels, 'O_MotorStep'), 1) = stepTriggers(:,1);
        if gd.Experiment.stim.control
            gd.Experiment.Triggers(beginMove1:stopMove1, strcmp(OutChannels, 'O_MotorStep'), 2) = stepTriggers(:,1);
            gd.Experiment.Triggers(beginMove1:stopMove1-1, strcmp(OutChannels, 'O_MotorDir'), 2) = 5;
        end
        
        % Move Out Stepper Motor
        beginMove2 = round((gd.Experiment.timing.preDuration + gd.Experiment.timing.stimDuration) * gd.Experiment.params.samplingFrequency) + 1;
        stopMove2 = beginMove2 + numStepTriggers - 1;
        gd.Experiment.Triggers(beginMove2:stopMove2, strcmp(OutChannels, 'O_MotorStep'), 1) = stepTriggers(:,1);
        gd.Experiment.Triggers(beginMove2:stopMove2-1, strcmp(OutChannels, 'O_MotorDir'), 1) = 5;
        if gd.Experiment.stim.control
            gd.Experiment.Triggers(beginMove2:stopMove2, strcmp(OutChannels, 'O_MotorStep'), 2) = stepTriggers(:,1);
        end
        
        % Trigger imaging computer
        gd.Experiment.Triggers([beginMove1, stopMove1, beginMove2, stopMove2], strcmp(OutChannels, 'O_2PTrigger'), :) = 1;
        
        % Trigger whisker tracking camera
        gd.Experiment.Triggers(beginMove1-ceil(DAQ.Rate/100):beginMove2, strcmp(OutChannels, 'O_WhiskerIllumination'), :) = 1;
        gd.Experiment.Triggers(beginMove1:ceil(DAQ.Rate/gd.Experiment.params.frameRateWT):beginMove2, strcmp(OutChannels, 'O_WhiskerTracker'), :) = 1;
        % gd.Experiment.Triggers(1, strcmp(OutChannels, 'O_WhiskerTracker'), :) = 1;        % mode 15 limited to 255 frames
        % gd.Experiment.Triggers(stopMove1, strcmp(OutChannels, 'O_WhiskerTracker'), :) = 1;
        
        % Stimulus info
        gd.Experiment.Stimulus = zeros(size(gd.Experiment.Triggers,1), 1);
        gd.Experiment.Stimulus(stopMove1:beginMove2-1) = 1;
        
        
        %% Initialize stimuli and linear motors
        temp = get(gd.Controls.positionTable, 'Data');
        ap = temp(:,1); ap(isnan(ap)) = [];
        if ~gd.Controls.mlMotorToggle.Value
            out = genDistribution(0,ap);
        else
            ml = temp(:,2); ml(isnan(ml)) = [];
            out = genDistribution(0,ap,ml);
        end
        % gd.Experiment.StimID = 1:gd.Experiment.stim.numPositions;
        % gd.Experiment.Position = 0:gd.Experiment.stim.distBetweenPositions:gd.Experiment.stim.numPositions*gd.Experiment.stim.distBetweenPositions;
        if gd.Experiment.stim.control
            gd.Experiment.StimID = [0, gd.Experiment.StimID];
            gd.Experiment.Position = [nan, gd.Experiment.Position];
        end
        gd.Experiment.Position
        H_LinearStage = serial(gd.Internal.LinearStage.APport, 'BaudRate', 9600);
        fopen(H_LinearStage);
        
        
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
        Triggers = gd.Experiment.Triggers;
        Positions = gd.Experiment.Position;
        ControlTrial = gd.Experiment.stim.control;
        currentBlockOrder = gd.Experiment.StimID;
        currentTrial = 0;
        TrialInfo = struct('StimID', [], 'Position', []);
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
        fclose(H_LinearStage);
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
                TrialInfo.Position(currentTrial) = Positions(randi(numStimuli-1))+1;            % thus choose random position from stimuli positions
                DAQ.queueOutputData(Triggers(:,:,2));
            else                                                                                % current trial is not control trial
                TrialInfo.Position(currentTrial) = Positions(TrialInfo.StimID(currentTrial));   % thus choose stimulus's position
                DAQ.queueOutputData(Triggers(:,:,1));
            end
            
            % Move linear motor
            moveLinearMotor(TrialInfo.Position(currentTrial), H_LinearStage); %move linear motor
            
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
function [TrialInfo, SaveFile] = stimPoleAzimuth(SaveFile, varargin)


%% Configure Settings
gd.Internal.ImagingType = 'scim'; % 'sbx' or 'scim'
gd.Internal.ImagingComp.ip = '128.32.173.30'; % SCANBOX ONLY: for UDP
gd.Internal.ImagingComp.port = 7000; % SCANBOX ONLY: for UDP
gd.Internal.LinearStage.APport = 'com3';
gd.Internal.LinearStage.MLport = 'com4';

Display.units = 'pixels';
Display.position = [100, 100, 1200, 500];

%% Initialize Experiment

gd.Internal.save.path = 'C:\Users\Resonant-2\OneDrive\StimData';
gd.Internal.save.base = '0000';
gd.Internal.save.depth = '000';
gd.Internal.save.index = '000';

gd.Experiment.saving.save = true;
gd.Experiment.saving.SaveFile = fullfile(gd.Internal.save.path, gd.Internal.save.base);
gd.Experiment.saving.DataFile = '';
gd.Experiment.saving.dataPrecision = 'uint16';

gd.Experiment.params.samplingFrequency = 30000;
gd.Experiment.params.numTrials = 5;
gd.Experiment.params.angleMove = 30;
gd.Experiment.params.frameRateWT = 200;

gd.Experiment.timing.stimDuration = 1.5; % in seconds
gd.Experiment.timing.ITI = 4.5; % in seconds
gd.Experiment.timing.randomITImax = 2; % in seconds

gd.Experiment.stim.control = true; % true or false (give control stimulus)
gd.Experiment.stim.blockShuffle = true; % true or false (shuffle trials)
gd.Experiment.stim.repeatBadTrials = false; % true or false (repeat non-running trials)
gd.Experiment.stim.speedThreshold = 100; % speed threshold for good running trials (deg/s)

% Properties for display or processing input data
gd.Internal.buffer.numTrials = 2; %4*gd.Experiment.params.samplingFrequency * (gd.Experiment.timing.stimDuration+gd.Experiment.timing.ITI);
gd.Internal.buffer.downSample = 30;


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
    'Name',                 'Stimulus: Vertical Pole',...
    'ToolBar',              'none',...
    'Units',                Display.units,...
    'Position',             Display.position);%,...
% 'KeyPressFcn',          @(hObject,eventdata)KeyPressCallback(hObject, eventdata, guidata(hObject)));

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
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,.1,1],...
    'Value',                gd.Experiment.saving.save,...
    'Callback',             @(hObject,eventdata)toggleSave(hObject, eventdata, guidata(hObject)));
if gd.Experiment.saving.save
    set(gd.Saving.save,'String','Saving','BackgroundColor',[0,1,0]);
else
    set(gd.Saving.save,'String','Save?','BackgroundColor',[1,0,0]);
end
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
    'Position',             [0, 0, .35, .8]);
% stepper motor CCW
gd.Controls.stepCCW = uicontrol(...
    'Style',                'pushbutton',...
    'String',               '<',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [0,.85,.35,.1],...
    'Callback',             @(hObject,eventdata)stepCCW(hObject, eventdata, guidata(hObject)));
% stepper motor degrees
gd.Controls.stepAngle = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.angleMove,...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.35,.85,.3,.1]);
% stepper motor CW
gd.Controls.stepCW = uicontrol(...
    'Style',                'pushbutton',...
    'String',               '>',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.65,.85,.35,.1],...
    'Callback',             @(hObject,eventdata)stepCW(hObject, eventdata, guidata(hObject)));
% stepper motor text
gd.Controls.stepperMotorText = uicontrol(...
    'Style',                'text',...
    'String',               'Move stepper motor (relative)',...
    'Parent',               gd.Controls.panel,...
    'HorizontalAlignment',  'center',...
    'Units',                'normalized',...
    'Position',             [0,.95,1,.05]);
% anterior-posterior linear motor position
gd.Controls.apMotorPos = uicontrol(...
    'Style',                'edit',...
    'String',               '0',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.35,.7,.3,.1],...
    'Callback',             @(hObject,eventdata)linearMotorPos(hObject, eventdata, guidata(hObject)));
% ap linear motor move
gd.Controls.apMotorMove = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'move',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'UserData',             'ap',...
    'Position',             [.65,.7,.35,.1],...
    'Callback',             @(hObject,eventdata)linearMotorMove(hObject, eventdata, guidata(hObject)));
% ap linear motor text
gd.Controls.mlMotorText = uicontrol(...
    'Style',                'text',...
    'String',               'Move anterior-posterior motor (absolute)',...
    'Parent',               gd.Controls.panel,...
    'HorizontalAlignment',  'center',...
    'Units',                'normalized',...
    'Position',             [0,.8,1,.05]);
% medial-lateral linear motor position
gd.Controls.mlMotorPos = uicontrol(...
    'Style',                'edit',...
    'String',               '0',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.35,.55,.3,.1],...
    'Callback',             @(hObject,eventdata)linearMotorPos(hObject, eventdata, guidata(hObject)));
% ml linear motor move
gd.Controls.mlMotorMove = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'move',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'UserData',             'ml',...
    'Position',             [.65,.55,.35,.1],...
    'Callback',             @(hObject,eventdata)linearMotorMove(hObject, eventdata, guidata(hObject)));
% ml linear motor text
gd.Controls.linearMotorText = uicontrol(...
    'Style',                'text',...
    'String',               'Move medial-lateral motor (absolute)',...
    'Parent',               gd.Controls.panel,...
    'HorizontalAlignment',  'center',...
    'Units',                'normalized',...
    'Position',             [0,.65,1,.05]);
% whisker imaging toggle
gd.Controls.wtToggle = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Imaging Whiskers',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.01,.4,.34,.1],...
    'UserData',             {[.94,.94,.94;0,1,0],'Image Whiskers?','Imaging Whiskers'},...
    'Callback',             @(hObject,eventdata)set(hObject,'BackgroundColor',hObject.UserData{1}(hObject.Value+1,:),'String',hObject.UserData{hObject.Value+2}));
% whisker imaging frame rate
gd.Controls.wtFrameRate = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.frameRateWT,...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.65,.4,.3,.1],...
    'Callback',             @(hObject,eventdata)ChangeWTFrameRate(hObject,eventdata,guidata(hObject)));
gd.Controls.wtFrameRateText = uicontrol(...
    'Style',                'text',...
    'String',               'Frame Rate:',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'HorizontalAlignment',  'right',...
    'Position',             [.35,.375,.29,.1]);
% random interval toggle
gd.Controls.randomITIToggle = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Random ITI?',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.01,.25,.34,.1],...
    'UserData',             {[.94,.94,.94;0,1,0],'Random ITI?','Randomizing ITI'},...
    'Callback',             @(hObject,eventdata)set(hObject,'BackgroundColor',hObject.UserData{1}(hObject.Value+1,:),'String',hObject.UserData{hObject.Value+2}));
% random interval max
gd.Controls.randomITImax = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.timing.randomITImax,...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'Position',             [.65,.25,.3,.1]);
gd.Controls.randomITIText = uicontrol(...
    'Style',                'text',...
    'String',               'Max Duration:',...
    'Parent',               gd.Controls.panel,...
    'Units',                'normalized',...
    'HorizontalAlignment',  'right',...
    'Position',             [.35,.225,.29,.1]);

% STIMULI
% panel
gd.Stimuli.panel = uipanel(...
    'Title',                'Run Experiment',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [.35, 0, .3, .8]);
% load stimuli
gd.Stimuli.load = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Load',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [0,.9,.5,.1],...
    'Callback',             @(hObject,eventdata)LoadStimuli(hObject, eventdata, guidata(hObject)));
% save stimuli
gd.Stimuli.save = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Save',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.9,.5,.1],...
    'Enable',               'off',...
    'Callback',             @(hObject,eventdata)SaveStimuli(hObject, eventdata, guidata(hObject)));
% position table
gd.Stimuli.positionTable = uitable(...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [0,.7,1,.2],...
    'RowName',              {'AP','ML'},...
    'ColumnName',           {'Active?','First (mm)','Last (mm)','# Steps'},...
    'ColumnEditable',       [true,true,true,true],...
    'ColumnFormat',         {'logical','numeric','numeric','numeric'},...
    'ColumnWidth',          {50,80,80,80},...
    'Data',                 [{true;false},num2cell([zeros(2,1),[17.5;7.5],[8;4]])],...
    'CellEditCallback',     @(hObject,eventdata)EditStimInputs(hObject, eventdata, guidata(hObject)));
% stimuli list
gd.Stimuli.list = uitable(...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'Normalized',...
    'Position',             [0,0,.6,.7],...
    'ColumnName',           {'AP Pos','ML Pos','Delete'},...
    'ColumnFormat',         {'numeric','numeric','logical'},...
    'ColumnEditable',       [true,true,true],...
    'ColumnWidth',          {75,75,50},...
    'CellEditCallback',     @(hObject,eventdata)EditStimuli(hObject, eventdata, guidata(hObject)));
% choose type
gd.Stimuli.type = uicontrol(...
    'Style',                'popupmenu',...
    'String',               {'Grid';'Random'},...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.6,.65,.4,.05]);
% toggle weighting
gd.Stimuli.weightToggle = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Weight sampling?',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.6,.6,.4,.05],...
    'UserData',             {'Weight sampling?','Weighting samples'},...
    'Callback',             @(hObject,eventdata)ChangeStimuliSampling(hObject,eventdata,guidata(hObject)));
% number of samples
gd.Stimuli.numSamples = uicontrol(...
    'Style',                'edit',...
    'String',               '# Samples',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.6,.5,.4,.1]);
% set AP sampling weights
gd.Stimuli.apDistWeights = uicontrol(...
    'Style',                'edit',...
    'String',               'AP Weights',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Enable',               'off',...
    'Position',             [.6,.4,.4,.1]);
% set AP sampling weights
gd.Stimuli.mlDistWeights = uicontrol(...
    'Style',                'edit',...
    'String',               'ML Weights',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Enable',               'off',...
    'Position',             [.6,.3,.4,.1]);
% create stimuli
gd.Stimuli.generateStimuli = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Generate Stimuli',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.6,.1,.4,.2],...
    'Callback',             @(hObject,eventdata)GenerateStimuli(hObject, eventdata, guidata(hObject)));
% view stimuli
gd.Stimuli.view = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'View',...
    'Parent',               gd.Stimuli.panel,...
    'Units',                'normalized',...
    'Position',             [.6,0,.4,.1],...
    'Enable',               'off',...
    'Callback',             @(hObject,eventdata)ViewStimuli(hObject, eventdata, guidata(hObject)));

% EXPERIMENT
% panel
gd.Run.panel = uipanel(...
    'Title',                'Run Experiment',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [.65, 0, .35, .8]);
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
    'Position',             [.075,.075,.9,.425]);

guidata(gd.fig, gd); % save guidata
CreateFilename(gd.Saving.FullFilename, [], gd);

% For timing
% gd=guidata(gd.fig);
% gd.Run.run.Value = true;
% RunExperiment(gd.Run.run,[],gd); % timing debugging
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

%% CONTROLS CALLBACKS
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
        fprintf('Anterior-Posterior motor: moved to position %.1f mm down track\n', str2num(gd.Controls.apMotorPos.String));
        clear H_LinearStage;
    case 'ml'
        H_LinearStage = serial(gd.Internal.LinearStage.MLport, 'BaudRate', 9600);
        fopen(H_LinearStage);
        moveLinearMotor(str2num(gd.Controls.mlMotorPos.String), H_LinearStage); %move linear motor
        fclose(H_LinearStage);
        fprintf('Medial-Lateral motor: Linear motor moved to position %.1f mm down track\n', str2num(gd.Controls.mlMotorPos.String));
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

function ChangeWTFrameRate(hObject,eventdata,gd)
newValue = str2num(hObject.String);
if newValue <= 0
    newValue = .0000001;
    hObject.String = num2str(newValue);
end
gd.Experiment.params.frameRateWT = newValue;
guidata(hObject, gd);
end


%% STIMULI CALLBACKS
function LoadStimuli(hObject, eventdata, gd)
% Select and load file
[f,p] = uigetfile({'*.stim';'*.mat'},'Select stim file to load',cd);
if isnumeric(f)
    return
end
load(fullfile(p,f), 'stimuli', '-mat');
fprintf('Loaded stimuli from: %s\n', fullfile(p,f));

% Load stimuli
gd.Stimuli.list.Data = num2cell(stimuli); % display stimuli
gd.Experiment.Position = stimuli; % Save loaded stimuli
guidata(hObject, gd);
end

function SaveStimuli(hObject, eventdata, gd)
% Determine file to save to
[f,p] = uiputfile({'*.stim';'*.mat'},'Save stimuli as?',cd);
if isnumeric(f)
    return
end
saveFile = fullfile(p,f);

% Save stimuli
stimuli = gd.Experiment.Position;
if ~exist(saveFile,'file')
    save(fullfile(p,f), 'stimuli', '-mat', '-v7.3');
else
    save(fullfile(p,f), 'stimuli', '-mat', '-append');
end
fprintf('Stimuli saved to: %s\n', fullfile(p,f));
end

function EditStimInputs(hObject, eventdata, gd)
val = hObject.Data{eventdata.Indices(1),eventdata.Indices(2)};
if eventdata.Indices(2)==1
    temp = {gd.Stimuli.apDistWeights, gd.Stimuli.mlDistWeights};
    if val && gd.Stimuli.weightToggle.Value
        set(temp{eventdata.Indices(1)}, 'Enable', 'on');
    else
        set(temp{eventdata.Indices(1)}, 'Enable', 'off');
    end
elseif eventdata.Indices(2)==2 && val < 0
    hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = 0;
elseif eventdata.Indices(2)==3 && val > 25
    hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = 25;
end
end

function ChangeStimuliSampling(hObject,eventdata,gd)
objs = [gd.Stimuli.apDistWeights, gd.Stimuli.mlDistWeights];
if hObject.Value
    toggle = [gd.Stimuli.positionTable.Data{:,1}];
    set(objs(toggle),'Enable','on');
else
    set(objs,'Enable','off');
end
end

function GenerateStimuli(hObject, eventdata, gd)
% Gather inputs
toggle = [gd.Stimuli.positionTable.Data{:,1}];
if ~any(toggle)
    error('Need at least one active axis!');
end
stimParams = cell2mat(gd.Stimuli.positionTable.Data(:,2:end));

% Determine # of samples
numSamples = str2num(gd.Stimuli.numSamples.String);
if isempty(numSamples)
    numSamples = 0;
end

% Determine weights
if gd.Stimuli.weightToggle.Value
    weights = {str2num(gd.Stimuli.apDistWeights.String), str2num(gd.Stimuli.mlDistWeights.String)};
    weights = weights(toggle);
else
    weights = {};
end

% Generate stimuli
if gd.Stimuli.type.Value == 1
    stimuli = genDistribution(numSamples,stimParams(toggle',:),'Distribution','grid','Weights',weights);
elseif gd.Stimuli.type.Value == 2
    stimuli = genDistribution(numSamples,stimParams(toggle',:),'Distribution','random','Weights',weights);
end
out = nan(size(stimuli,1),numel(toggle));
out(:,toggle) = stimuli;

gd.Experiment.Position = out; % Update registry
guidata(hObject, gd);
gd.Stimuli.list.Data = out; % display stimuli
gd.Stimuli.view.Enable = 'on';
end

function ViewStimuli(hObject, eventdata, gd)
stimuli = gd.Experiment.Position;
if any(isnan(stimuli(:)))
    stimuli(:,any(isnan(stimuli),1))=[];
    numDims = 1;
else
    numDims=2;
end
nbins = 20;
figure;
if numDims == 1
    subplot(1,2,1); plot(zeros(size(stimuli,1),1),stimuli,'k.'); ylabel('Value');
    subplot(1,2,2); histogram(stimuli,nbins); xlabel('Position'); ylabel('count');
elseif numDims == 2
    subplot(1,3,1); plot(stimuli(:,2),stimuli(:,1),'k.'); xlabel('ML'); ylabel('AP');
    subplot(1,3,2); histogram(stimuli(:,2),nbins); xlabel('ML Position'); ylabel('count'); 
    subplot(1,3,3); histogram(stimuli(:,1),nbins); xlabel('AP Position'); ylabel('count'); 
end
end

function EditStimuli(hObject, eventdata, gd)
if eventdata.Indices(2)==3
    gd.Experiment.Position(eventdata.Indices(1),:) = []; % remove stimulus
    guidata(hObject, gd);
    hObject.Data(eventdata.Indices(1),:) = []; % update display
end
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
        % Imaging Computer Trigger (for timing)
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line0','OutputOnly');
        DAQ.Channels(id).Name = 'O_2PTrigger';
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line1','InputOnly');
        DAQ.Channels(id).Name = 'I_FrameCounter';
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
        if gd.Controls.wtToggle.Value
            [~,id] = DAQ.addDigitalChannel('Dev1','port0/line17','OutputOnly');
            DAQ.Channels(id).Name = 'O_WhiskerTracker';
            [~,id] = DAQ.addDigitalChannel('Dev1','port0/line2','OutputOnly');
            DAQ.Channels(id).Name = 'O_WhiskerIllumination';
            [~,id] = DAQ.addDigitalChannel('Dev1','port0/line18','InputOnly');
            DAQ.Channels(id).Name = 'I_WhiskerTracker';
        end
        % Cleanup
        DAQChannels = {DAQ.Channels(:).Name};
        OutChannels = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));
        InChannels = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
        
%         % Add clock
%         daqClock = daq.createSession('ni');
%         daqClock.addCounterOutputChannel('Dev1',0,'PulseGeneration');
%         clkTerminal = daqClock.Channels(1).Terminal;
%         daqClock.Channels(1).Frequency = DAQ.Rate;
%         daqClock.IsContinuous = true;
%         daqClock.startBackground;
%         DAQ.addClockConnection('External',['Dev1/' clkTerminal],'ScanClock');
        
        % Add Callbacks
        DAQ.addlistener('DataRequired', @QueueData);    % create listener for queueing trials
        DAQ.NotifyWhenScansQueuedBelow = DAQ.Rate-1; % queue more data when less than a second of data left
        DAQ.addlistener('DataAvailable', @SaveDataIn);  % create listener for when data is returned
        % DAQ.NotifyWhenDataAvailableExceeds = round(DAQ.Rate/100);

        % Compute timing of each trial
        gd.Experiment.timing.ITI = str2num(gd.Run.ITI.String);
        gd.Experiment.timing.stimDuration = str2num(gd.Run.stimDur.String);
        if gd.Controls.randomITIToggle.Value
            gd.Experiment.timing.randomITImax = str2num(gd.Controls.randomITImax.String);
        else
            gd.Experiment.timing.randomITImax = false;
        end
        gd.Experiment.timing.trialDuration = gd.Experiment.timing.stimDuration + gd.Experiment.timing.ITI;
        gd.Experiment.timing.numScansPerTrial = ceil(gd.Experiment.params.samplingFrequency * gd.Experiment.timing.trialDuration);
        
        
        %% Determine stimuli
        
        % Determine axes
        gd.Experiment.stim.activeAxes = false(1,2);
        for aindex = 1:2
            if ~any(isnan(gd.Experiment.Position(:,aindex)))
                gd.Experiment.stim.activeAxes(aindex) = true;
            end
        end
        
        % Determine stimulus IDs
        gd.Experiment.StimID = 1:size(gd.Experiment.Position,1);
        
        % Determine if presenting control stimulus
        gd.Experiment.stim.control = gd.Run.control.Value;
        if gd.Experiment.stim.control
            gd.Experiment.StimID = [0, gd.Experiment.StimID];
            gd.Experiment.Position = [nan(1,2); gd.Experiment.Position];
        end
        
        % Grab UI variables
        gd.Experiment.stim.blockShuffle = gd.Run.shuffle.Value;
        gd.Experiment.stim.repeatBadTrials = gd.Run.repeatBadTrials.Value;
        
        
        %% Create triggers
        
        % Initialize triggers
        gd.Experiment.Triggers = zeros(gd.Experiment.timing.numScansPerTrial, numel(OutChannels), gd.Experiment.stim.control+1);
        stepTriggers = moveStepperMotor(gd.Experiment.params.angleMove, gd.Experiment.params.samplingFrequency, 1, 2, false, 2);
        numStepTriggers = size(stepTriggers, 1);
        
        % Determine timing
        numITITriggers = floor(gd.Experiment.params.samplingFrequency * gd.Experiment.timing.ITI);    
        if numITITriggers < 2*numStepTriggers
            error('Not enough time for bar to move in and out');
        end
        startTrig = numITITriggers-numStepTriggers+1;                       % move during end of ITI but leave room for bar moving out
        endTrig = gd.Experiment.timing.numScansPerTrial-numStepTriggers+1;  % bar moving our occurs during ITI but is queued with previous trial
        
        % Adjust Callback timing so next trial is queued right after previous trial starts
        DAQ.NotifyWhenScansQueuedBelow = gd.Experiment.timing.numScansPerTrial - startTrig;
        
        % Move In Stepper Motor
        stopMove1 = startTrig-1;
        beginMove1 = stopMove1-numStepTriggers+1;
        gd.Experiment.Triggers(beginMove1:stopMove1, strcmp(OutChannels,'O_MotorStep'), 1) = stepTriggers(:,1);
        if gd.Experiment.stim.control
            gd.Experiment.Triggers(beginMove1:stopMove1, strcmp(OutChannels,'O_MotorStep'), 2) = stepTriggers(:,1);
            gd.Experiment.Triggers(beginMove1:stopMove1-1, strcmp(OutChannels,'O_MotorDir'), 2) = 5;
        end
        
        % Move Out Stepper Motor
        beginMove2 = gd.Experiment.timing.numScansPerTrial-numStepTriggers+1;
        stopMove2 = gd.Experiment.timing.numScansPerTrial;
        gd.Experiment.Triggers(beginMove2:stopMove2, strcmp(OutChannels,'O_MotorStep'), 1) = stepTriggers(:,1);
        gd.Experiment.Triggers(beginMove2:stopMove2-1, strcmp(OutChannels,'O_MotorDir'), 1) = 5;
        if gd.Experiment.stim.control
            gd.Experiment.Triggers(beginMove2:stopMove2, strcmp(OutChannels,'O_MotorStep'), 2) = stepTriggers(:,1);
        end
        
        % Trigger imaging computer on every single trial
        gd.Experiment.Triggers([startTrig, endTrig], strcmp(OutChannels,'O_2PTrigger')) = 1; % trigger at beginning and end of stimulus
        % gd.Experiment.Triggers(startTrig:endTrig-1, strcmp(OutChannels,'O_2PTrigger')) = 1; % high during whole stimulus
        
        % Trigger whisker tracking camera on every single trial
        if gd.Controls.wtToggle.Value
            if gd.Experiment.timing.ITI >= 0.01;
                gd.Experiment.Triggers(startTrig-ceil(DAQ.Rate/100):endTrig, strcmp(OutChannels,'O_WhiskerIllumination')) = 1; % start LED a little before the start of imaging
                gd.Experiment.Triggers(startTrig:ceil(DAQ.Rate/gd.Experiment.params.frameRateWT):endTrig, strcmp(OutChannels,'O_WhiskerTracker')) = 1; % image during stimulus period
            else
                gd.Experiment.Triggers(:, strcmp(OutChannels,'O_WhiskerIllumination')) = 1; % image during entire time
                gd.Experiment.Triggers(1:ceil(DAQ.Rate/gd.Experiment.params.frameRateWT):endTrig, strcmp(OutChannels,'O_WhiskerTracker')) = 1; % image during entire time
            end
            % gd.Experiment.Triggers(1, strcmp(OutChannels, 'O_WhiskerTracker')) = 1;        % mode 15 limited to 255 frames
            % gd.Experiment.Triggers(stopMove1, strcmp(OutChannels, 'O_WhiskerTracker')) = 1;
        end
        
        % Build up vector to display when stimulus is present
        gd.Experiment.Stimulus = zeros(gd.Experiment.timing.numScansPerTrial,1);
        gd.Experiment.Stimulus(startTrig:endTrig-1) = 1;
        
        
        %% Initialize stimuli and linear motors
        H_LinearStage = [];
        if gd.Experiment.stim.activeAxes(1)
            H_LinearStage(1) = serial(gd.Internal.LinearStage.APport, 'BaudRate', 9600);
            fopen(H_LinearStage(1));
        end
        if gd.Experiment.stim.activeAxes(2)
            H_LinearStage(2) = serial(gd.Internal.LinearStage.MLport, 'BaudRate', 9600);
            fopen(H_LinearStage(2));
        end
        
        
        %% Initialize imaging session (scanbox only)
        gd.Experiment.ImagingType = gd.Saving.imagingType.String{gd.Saving.imagingType.Value};
        if strcmp(gd.Experiment.ImagingType, 'sbx')
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
        ActiveAxes = find(Experiment.stim.activeAxes);
        numStimuli = numel(Experiment.StimID);
        numStimuliCurrentBlock = numStimuli;
        Triggers = Experiment.Triggers;
        Positions = Experiment.Position;
        ControlTrial = Experiment.stim.control;
        BlockShuffle = Experiment.stim.blockShuffle;
        currentBlockOrder = Experiment.StimID;
        currentTrial = 0;
        TrialInfo = struct('StimID', [], 'Running', [], 'RunSpeed', []);
        saveOut = Experiment.saving.save;
        Stimulus = Experiment.Stimulus;
        ExperimentReachedEnd = false; % boolean to see if max trials has been reached
        numScansPerTrial = Experiment.timing.numScansPerTrial;
        MaxRandomScans = floor(gd.Experiment.timing.randomITImax*Experiment.params.samplingFrequency)+1;
        
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
        if strcmp(gd.Experiment.ImagingType, 'sbx')
            fprintf(H_Scanbox,'G'); %go
            pause(5);
        end
        
        % Start experiment
        Experiment.timing.start = datestr(now);
        QueueData();
        DAQ.startBackground;
        
        
        %% During Experiment
        while DAQ.IsRunning
            pause(1);
        end
        Experiment.timing.finish = datestr(now);
        
        %% End Experiment
        if strcmp(gd.Experiment.ImagingType, 'sbx')
            fprintf(H_Scanbox,'S'); %stop
            fclose(H_Scanbox);
        end
        if saveOut
            save(SaveFile, 'Experiment', '-append'); % update with "Experiment.timing.finish" info
            fclose(H_DataFile);                      % close binary file
            gd.Saving.index.String = sprintf('%03d',str2double(gd.Saving.index.String) + 1); % update file index for next experiment
            CreateFilename(gd.Saving.FullFilename, [], gd); % update filename for next experiment
        end
        
        % Close connections
        for findex = 1:size(H_LinearStage,2)
            fclose(H_LinearStage(findex));
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
            if strcmp(gd.Experiment.ImagingType, 'sbx')
                fprintf(H_Scanbox,'S'); %stop
                fclose(H_Scanbox);
            end
            for findex = 1:size(H_LinearStage,2)
                fclose(H_LinearStage(findex));
            end
        end
        clear DAQ H_Scanbox H_LinearStage
        
        % Rethrow error
        rethrow(ME);
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
        BufferStim = BufferStim(numScansReturned+1:end); % remove corresponding data from stimulus buffer
        
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
            
            % Move motor(s) into position for next trial
            if numel(TrialInfo.StimID)>=RunIndex        % next trial has already been queued (necessary to move motors)
                if TrialInfo.StimID(RunIndex)==0        % next trial is a control trial
                    temp = randi([ControlTrial+1,numStimuli]);
                    for index=ActiveAxes
                        moveLinearMotor(Positions(temp,index), H_LinearStage(index)); %move motor
                    end
                else                                    % next trial is not a control trial
                    for index=ActiveAxes
                        moveLinearMotor(Positions(TrialInfo.StimID(RunIndex)+ControlTrial,index), H_LinearStage(index)); %move motor
                    end
                end
            end
            
        end
        
        % Display new data and stimulus info
        plot(hAxes, numBufferScans:-dsamp:1, dx_dt, 'b-', numBufferScans:-dsamp:1, 500*(currentStim>0), 'r-');
        ylim([-100,600]);
        
    end %SaveDateIn

%% Callback: QueueOutputData
    function QueueData(src,eventdata)
        
        % Queue next trial
        if  hObject.Value && (currentTrial < str2double(numTrialsObj.String) || ~isempty(StimuliToRepeat)) % user hasn't quit, and requested # of trials hasn't been reached or no trials need to be repeated
            
            % Update index
            currentTrial = currentTrial + 1;
            ExperimentReachedEnd = false;
            
            if currentTrial <= str2double(numTrialsObj.String) %haven't reached end of experiment
                blockIndex = rem(currentTrial-1, numStimuliCurrentBlock)+1;
                if BlockShuffle && blockIndex == 1                                  % if starting new block, shuffle the stimuli order
                    % numStimuliCurrentBlock = numStimuli;                            % reset size of block (changes if repeating bad trials)
                    currentBlockOrder = currentBlockOrder(randperm(numStimuli));    % shuffle block
                end
                TrialInfo.StimID(currentTrial) = currentBlockOrder(blockIndex);     % determine StimID for current trial
            else %trials need to be repeated
                TrialInfo.StimID(currentTrial) = StimuliToRepeat(1);                % present first trial in repeat queue
                StimuliToRepeat(1) = [];                                            % remove trial from repeat queue
            end
            
            % Move motors into position for first trial
            if currentTrial == 1
                if TrialInfo.StimID(currentTrial)==0        % first trial is a control trial
                    temp = randi([ControlTrial+1,numStimuli]);
                    for index=ActiveAxes
                        moveLinearMotor(Positions(temp,index), H_LinearStage(index)); %move motor
                    end
                else                                        % first trial is not a control trial
                    for index=ActiveAxes
                        moveLinearMotor(Positions(TrialInfo.StimID(1)+ControlTrial,index), H_LinearStage(index)); %move motor
                    end
                end
            end
            
            % Queue triggers
            if TrialInfo.StimID(currentTrial) ~= 0              % current trial is not control trial
                if MaxRandomScans == 1                          % do not add random ITI
                    DAQ.queueOutputData(Triggers(:,:,1));       % queue normal stim triggers
                else
                    TrialInfo.numRandomScansPost(currentTrial) = randi(MaxRandomScans)-1; % determine amount of extra scans to add
                    DAQ.queueOutputData(cat(1,Triggers(:,:,1),zeros(TrialInfo.numRandomScansPost(currentTrial),numel(OutChannels)))); % queue normal stim triggers with extra scans
                end
            else                                                % current trial is control trial
                if MaxRandomScans == 1                          % do not add random ITI
                    DAQ.queueOutputData(Triggers(:,:,2));       % queue normal control triggers
                else
                    TrialInfo.numRandomScansPost(currentTrial) = randi(MaxRandomScans)-1; % determine amount of extra scans to add
                    DAQ.queueOutputData(cat(1,Triggers(:,:,2),zeros(TrialInfo.numRandomScansPost(currentTrial),numel(OutChannels)))); % queue normal control triggers with extra scans
                end
            end
            
            % Update buffer
            if MaxRandomScans == 1
                BufferStim = cat(1, BufferStim, Stimulus*currentTrial);
            else
                BufferStim = cat(1, BufferStim, Stimulus*currentTrial, zeros(TrialInfo.numRandomScansPost(currentTrial),1));
            end
            
            % Update information
            fprintf('\nQueued trial %d: stimulus %d', currentTrial, TrialInfo.StimID(currentTrial));
            if saveOut
                save(SaveFile, 'TrialInfo', '-append');
            end
            
        elseif ~ExperimentReachedEnd
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
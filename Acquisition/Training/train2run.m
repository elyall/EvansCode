function [TrialInfo, SaveFile] = train2run(SaveFile, varargin)


%% Configure Settings
Display.units = 'pixels';
Display.position = [100, 100, 800, 400];

%% Initialize Experiment

gd.Internal.save.path = fullfile(cd,datestr(now,'yymmdd'));
gd.Internal.save.filename = '0000';

gd.Experiment.saving.save = false;
gd.Experiment.saving.SaveFile = fullfile(gd.Internal.save.path, gd.Internal.save.filename);
gd.Experiment.saving.DataFile = '';
gd.Experiment.saving.dataPrecision = 'uint16';

gd.Experiment.params.samplingFrequency = 5000;
gd.Experiment.params.threshold = 100; %deg/sec
gd.Experiment.params.duration = 60; %minutes

gd.Experiment.punish.delay = 10; %sec
gd.Experiment.punish.volume = 1;
gd.Experiment.punish.duration = 1; %sec

gd.Internal.buffer.seconds = 6;
gd.Internal.buffer.downSample = 20;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'duration','Duration'}
                gd.Experiment.params.duration = varargin{index+1};
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
    'Name',                 'Train 2 Run',...
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
    'Position',             [.1,0,.1,1],...
    'Callback',             @(hObject,eventdata)ChooseDir(hObject, eventdata, guidata(hObject)));
gd.Saving.FullFilename = uicontrol(...
    'Style',                'text',...
    'String',               '',...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.2,.05,.8,.2],...
    'Callback',             @(hObject,eventdata)CreateFilename(hObject, eventdata, guidata(hObject)));
% filename input
gd.Saving.filename = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Internal.save.filename,...
    'Parent',               gd.Saving.panel,...
    'Units',                'normalized',...
    'Position',             [.2,.3,.8,.7],...
    'Callback',             @(hObject,eventdata)SetFilename(hObject, eventdata, guidata(hObject)));
% % TRAINING
% panel
gd.Run.panel = uipanel(...
    'Title',                'Run Session',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [0, 0, 1, .8]);
% duration input
gd.Run.dur = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.duration,...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [0,.8,.2,.2]);
% run button
gd.Run.run = uicontrol(...
    'Style',                'togglebutton',...
    'String',               'Run?',...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.2,.8,.6,.2],...
    'Callback',             @(hObject,eventdata)RunExperiment(hObject, eventdata, guidata(hObject)));
% threshold input
gd.Run.thresh = uicontrol(...
    'Style',                'edit',...
    'String',               gd.Experiment.params.threshold,...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.8,.8,.2,.2]);
% running speed axes
gd.Run.runSpeedAxes = axes(...
    'Parent',               gd.Run.panel,...
    'Units',                'normalized',...
    'Position',             [.1,0,.6,.8]);
% average speed
gd.Run.meanText = uicontrol(...
    'Style',                'text',...
    'String',               'Mean run speed:',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlign',      'right',...
    'Units',                'normalized',...
    'Position',             [.7,.6,.2,.2]);
gd.Run.mean = uicontrol(...
    'Style',                'text',...
    'String',               '',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlign',      'left',...
    'Units',                'normalized',...
    'Position',             [.9,.6,.1,.2]);
% percent above threshold
gd.Run.prcAboveText = uicontrol(...
    'Style',                'text',...
    'String',               '% above threshold:',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlign',      'right',...
    'Units',                'normalized',...
    'Position',             [.7,.4,.2,.2]);
gd.Run.prcAbove = uicontrol(...
    'Style',                'text',...
    'String',               '',...
    'Parent',               gd.Run.panel,...
    'HorizontalAlign',      'left',...
    'Units',                'normalized',...
    'Position',             [.9,.4,.2,.2]);
% % average duration above threshold
% gd.Run.avgDurAboveText = uicontrol(...
%     'Style',                'text',...
%     'String',               'Avg Burst Duration:',...
%     'Parent',               gd.Run.panel,...
%     'HorizontalAlign',      'right',...
%     'Units',                'normalized',...
%     'Position',             [.7,.2,.2,.2]);
% gd.Run.avgDurAbove = uicontrol(...
%     'Style',                'text',...
%     'String',               '',...
%     'Parent',               gd.Run.panel,...
%     'HorizontalAlign',      'left',...
%     'Units',                'normalized',...
%     'Position',             [.9,.2,.2,.2]);
% % variance
% gd.Run.varText = uicontrol(...
%     'Style',                'text',...
%     'String',               'Variance:',...
%     'Parent',               gd.Run.panel,...
%     'HorizontalAlign',      'right',...
%     'Units',                'normalized',...
%     'Position',             [.7,0,.2,.2]);
% gd.Run.var = uicontrol(...
%     'Style',                'text',...
%     'String',               '',...
%     'Parent',               gd.Run.panel,...
%     'HorizontalAlign',      'left',...
%     'Units',                'normalized',...
%     'Position',             [.9,0,.2,.2]);

guidata(gd.fig, gd); % save guidata
CreateFilename(gd.Saving.FullFilename, [], gd);
end

%% Callbacks
function ChooseDir(hObject, eventdata, gd)
temp = uigetdir(gd.Internal.save.path, 'Choose directory to save to');
if ischar(temp)
    gd.Internal.save.path = temp;
    guidata(hObject, gd);
end
CreateFilename(gd.Saving.FullFilename, [], gd);
end

function SetFilename(hObject, eventdata, gd)
gd.Internal.save.filename = hObject.String;
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
gd.Experiment.saving.SaveFile = fullfile(gd.Internal.save.path, strcat(gd.Saving.filename.String, '.training'));
hObject.String = gd.Experiment.saving.SaveFile;
guidata(hObject, gd);
if exist(gd.Experiment.saving.SaveFile, 'file')
    hObject.BackgroundColor = [1,0,0];
else
    hObject.BackgroundColor = [.94,.94,.94];
end
end

%% CONTORLS CALLBACKS


%% RUN EXPERIMENT
function RunExperiment(hObject, eventdata, gd)

if hObject.Value
    try
        hObject.BackgroundColor = [0,0,0];
        hObject.ForegroundColor = [1,1,1];
        hObject.String = 'Stop';
        
        
        %% Determine filenames to save to
        if gd.Saving.save.Value
            gd.Experiment.saving.save = true;
            % mat file
            if exist(gd.Experiment.saving.SaveFile, 'file')
                answer = questdlg(sprintf('File already exists! Continue?\n%s', gd.Experiment.saving.SaveFile), 'Overwrite file?', 'Yes', 'No', 'No');
                if strcmp(answer, 'No')
                    return
                end
            end
            SaveFile = gd.Experiment.saving.SaveFile;
            % bin file
            gd.Experiment.DataFile = strcat(gd.Experiment.saving.SaveFile(1:end-4), '.bin');
            DataFile = gd.Experiment.DataFile;
        else
            gd.Experiment.saving.save = false;
        end
        
        
        %% Grab UI data
        gd.Experiment.params.duration = str2num(gd.Run.dur.String);
        gd.Experiment.params.threshold = str2num(gd.Run.thresh.String);
        
        %% Initialize NI-DAQ session
        gd.Internal.daq = [];
        
        DAQ = daq.createSession('ni'); % initialize session
        DAQ.IsContinuous = true; % set session to be continuous (call's 'DataRequired' listener)
        DAQ.Rate = gd.Experiment.params.samplingFrequency; % set sampling frequency
        gd.Experiment.params.samplingFrequency = DAQ.Rate; % the actual sampling frequency is rarely perfect from what is input
        
        % Add ports
        % Output port
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line16','OutputOnly');
        DAQ.Channels(id).Name = 'O_Speaker';
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line17','OutputOnly');
        DAQ.Channels(id).Name = 'O_Air';
        % Running Wheel
        [~,id] = DAQ.addDigitalChannel('Dev1','port0/line5:7','InputOnly');
        DAQ.Channels(id(1)).Name = 'I_RunWheelA';
        DAQ.Channels(id(2)).Name = 'I_RunWheelB';
        DAQ.Channels(id(3)).Name = 'I_RunWheelIndex';
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
        DAQ.NotifyWhenScansQueuedBelow = DAQ.Rate*gd.Experiment.punish.duration;
        % Add DataIn callback
        DAQ.addlistener('DataAvailable', @SaveDataIn);
        DAQ.NotifyWhenDataAvailableExceeds = DAQ.Rate/100;
        
        
        %% Create triggers
        gd.Experiment.Triggers = zeros(round(gd.Experiment.params.samplingFrequency), numel(OutChannels));
        
        gd.Experiment.PunishTriggers = 
        %% Initialize saving
        Experiment = gd.Experiment;
        if gd.Experiment.saving.save
            save(SaveFile, 'DAQChannels', 'Experiment', '-mat', '-v7.3');
            H_DataFile = fopen(gd.Experiment.saving.DataFile, 'w');
        end
        
        
        %% Initialize shared variables (only share what's necessary)
        
        % Necessary variables
        numScansQueued = 0;
        numScansReturned = 0;
        samplingFrequency = double(Experiment.params.samplingFrequency);
        numTotalScans = round(samplingFrequency*Experiment.params.duration*60);
        
        Triggers = Experiment.Triggers;
        saveOut = Experiment.saving.save;

        % Variables if saving input data
        if saveOut
            Precision = Experiment.saving.dataPrecision;
        end
        
        % Variables for calculating and displaying running speed
        RunChannelIndices = [find(strcmp(InChannels, 'I_RunWheelB')),find(strcmp(InChannels,'I_RunWheelA'))];
        numBufferScans = round(gd.Internal.buffer.seconds * samplingFrequency);
        DataInBuffer = zeros(numBufferScans, 2);
        dsamp = gd.Internal.buffer.downSample;
        dsamp_Fs = samplingFrequency / dsamp;
        smooth_win = gausswin(dsamp_Fs, 23.5/2);
        smooth_win = smooth_win/sum(smooth_win);
        sw_len = length(smooth_win);
        d_smooth_win = [0;diff(smooth_win)]/(1/dsamp_Fs);
        hAxes = gd.Run.runSpeedAxes;
        numScansPerReturn = double(DAQ.NotifyWhenDataAvailableExceeds);
        
        % Metrics
        threshold = Experiment.params.threshold;
        h.mean = gd.Run.mean;
        h.prcAbove = gd.Run.prcAbove;
%         h.avgDurAbove = gd.Run.avgDurAbove;
%         h.var = gd.Run.var;
        metrics.dist = 0;
        metrics.numScansAbove = 0;
        metrics.avgDurAbove = 0;
        metrics.var = 0;
        
        %% Start Experiment
        fprintf('\n000 / 000 minutes completed, 000 minutes remain');
        QueueData();
        DAQ.startBackground;
        
        
        %% During Experiment
        while DAQ.IsRunning
            pause(0.1);
            Experiment.params.duration = str2num(gd.Run.dur.String);
            numTotalScans = round(Experiment.params.samplingFrequency*Experiment.params.duration*60);
        end
        
        
        %% End Experiment
        if saveOut
            save(SaveFile, 'Experiment', '-append'); % update with "Experiment.timing.finish" info
            fclose(H_DataFile);                      % close binary file
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
        
        % Update state
        numScansReturned = numScansReturned + numScansPerReturn;
        fprintf(repmat('\b',1,47));
        fprintf('%3.0f / %3.0f minutes completed, %3.0f minutes remain', numScansReturned/samplingFrequency/60, numTotalScans/samplingFrequency/60, (numTotalScans-numScansReturned)/samplingFrequency/60);
                
        % Refresh buffer
        DataInBuffer = cat(1, DataInBuffer(numScansPerReturn+1:end,:), eventdata.Data(:,RunChannelIndices));   % concatenate new data and remove old data
        
        % Convert entire buffer of pulses to run speed
        Data = [0;diff(DataInBuffer(:,1))>0];       % gather pulses' front edges
        Data(all([Data,DataInBuffer(:,2)],2)) = -1; % set backwards steps to be backwards
        x_t = downsample(cumsum(Data), dsamp);      % convert pulses to counter data & downsample data to speed up computation
        x_t = padarray(x_t, sw_len, 'replicate');   % pad for convolution
        dx_dt = conv(x_t, d_smooth_win, 'same');    % perform convolution
        dx_dt([1:sw_len,end-sw_len+1:end]) = [];    % remove values produced by padding the data
        % dx_dt = dx_dt * 360/360;                  % convert to degrees (360 pulses per 360 degrees)
        
        % Display new data and stimulus info
        plot(hAxes, numBufferScans:-dsamp:1, dx_dt, 'b-');
        ylim([-100,800]);
        
        % Compute mean
        metrics.dist = metrics.dist + sum(Data(end-numScansPerReturn+1:end));
        h.mean.String = sprintf('%.0f deg/s', metrics.dist/(numScansReturned/samplingFrequency));
        % Compute percent above threshold
        metrics.numScansAbove = metrics.numScansAbove + sum(dx_dt(end-ceil(numScansPerReturn/dsamp)+1:end)>threshold);
        h.prcAbove.String = sprintf('%.1f%%',metrics.numScansAbove/ceil(numScansReturned/dsamp)*100);
        
    end %SaveDateIn

%% Callback: QueueOutputData
    function QueueData(src,eventdata)
        
        if numScansQueued < numTotalScans && hObject.Value
            
                DAQ.queueOutputData(Triggers);

                numScansQueued = numScansQueued + size(Triggers,1);
                
        elseif ~hObject.Value
%             fprintf('\nComplete: user quit\n');
            
        elseif numScansQueued >= numTotalScans
%             fprintf('\nComplete: max time reached\n');
        end
    end

end
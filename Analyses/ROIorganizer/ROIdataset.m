classdef ROIdataset < handle
    
    properties (Hidden = true)
        DataSetID                   % define dataset's unique ID
    end
    
    properties
        ImagingFile                 % filename of imaging data
        ExperimentFile              % filename of experiment data
        frameSize                   % [H, W]
        date                        % date
        depth                       % depth of imaging
        fileLocation = nan(1,3)     % [frame #, z-level, channel]
        Vertices                    % (Nvertices x 2) [X, Y]
        Mask                        % Pixels to average over
        NeuropilMask                % Pixels to average over for Neuropil
        RawData                     % Extracted flourescence
        RawNeuropil                 % Extracted neuropil flourescence
        StimID                      % (Ntrials x 1)
        Labels = {};                % (cellstr) User specified labels
%         ParsedStimIndices           % (Ntrials x 2) [before, after]
        ParsedFrameIndices          % (Ntrials x 2) [first frame, last frame]
        dFoF                        % (Ntrials x Nframes)
%         ParsedTrialID               % (Ntrials x 1)
        UserData                    % user specified information
    end %main properties
        
    properties (Dependent = true)
        Pixels
        Centroid
        Data            % (Ntrials x Nframes)
        Neuropil        % (Ntrials x Nframes)
%         Stimulus        % (Ntrials x Nframes) logical 1 when stimulus present, 0 when stimulus not present
    end %dependent properties
    
    methods
        
        % Constructor
        function obj = ROIdataset(ImagingFile, varargin)
            obj.DataSetID = now;                % define dataset's unique ID
            
            % Define imaging file dataset is associated with
            if nargin > 0
                obj.ImagingFile = ImagingFile;
            end
            
            % Define associated information
            for index = 1:2:numel(varargin)
                try
                    obj.(varargin{index}) = varargin{index+1};
                end
            end
            
        end %constructor
        
        function set.ImagingFile(obj,val)
            if ~exist(val, 'file')
                error('File does not exist: %s', val);
            end
            obj.ImagingFile = val;
            
            load([val(1:end-3),'mat'], 'info')
            
            Config = load2PConfig(val);
            obj.FrameSize = [Config.Height, Config.Width];
            obj.FileLocation(Config.size(3:end)==1) = 1;
        end
        
        function set.ExperimentFile(obj,val)
            if ~exist(val, 'file')
                error('File does not exist: %s', val);
            end
            obj.ExperimentFile = val;
            
            load(val, 'Experiment', 'TrialInfo', '-mat');
            obj.Date = Experiment.timing.init;
            obj.StimID = TrialInfo.StimID;
        end
                
        function set.Vertices(obj,val)
            obj.Vertices = val;
        end
        
        function val = get.Pixels(obj)
            val = poly2mask(obj.Vertices(:,1), obj.Vertices(:,2), obj.FrameSize(1), obj.FrameSize(2));
        end
        
        function val = get.Centroid(obj)
            temp = regionprops(obj.Pixels, 'centroid');
            val = temp.Centroid;
        end
        
        function obj = parse(obj, varargin)
            
            % Set Defaults
            NumTriggersPerTrial = 4;
            TrialID = [1 inf];
            NumFramesBefore = nan;
            NumFramesAfter = nan;
            
            % Parse inputs
            for index = 1:2:numel(varargin)
                switch varargin{index}
                    case 'Trials'
                        TrialID = varargin{index+1};
                    case 'NumFramesBefore'
                        NumFramesBefore = varargin{index+1};
                    case 'NumFramesAfter'
                        NumFramesAfter = varargin{index+1};
                    case 'NumTriggersPerTrial'
                        NumTriggersPerTrial = varargin{index+1};
                end
            end
            numTrials = numel(obj.StimID);
            
            % Load trigger information
            load([obj.ImagingFile(1:end-3),'mat'], 'info');
            FrameTriggers = reshape(info.frame+1, NumTriggersPerTrial, numTrials);
            
            % Determine trials to parse
            if TrialID(end) == inf
                TrialID = [TrialID(1:end-1), TrialID(end-1)+1:numTrials];
            end
            
            if numel(NumFramesBefore)==1
                NumFramesBefore = repmat(NumFramesBefore, numTrials, 1);
            end
            if numel(NumFramesAfter)==1
                NumFramesAfter = repmat(NumFramesAfter, numTrials, 1);
            end
            
            % Determine frame indices for each trial
            obj.ParsedFrameIndices = nan(numTrials, 2);
            for tindex = TrialID
                
                % Start of trial
                if isnan(NumFramesBefore(tindex))
                    if tindex == 1
                        firstframe = 1;
                    else
                        firstframe = FrameTriggers(4,tindex-1);
                    end
                else
                    firstframe = FrameTriggers(2,tindex) - NumFramesBefore(tindex);
                    if firstframe < 1
                        firstframe = 1;
                    end
                end
                obj.ParsedFrameIndices(tindex, 1) = firstframe;
                                
                % End of trial
                if isnan(NumFramesAfter(tindex))
                    obj.ParsedFrameIndices(tindex, 2) = FrameTriggers(3,tindex);
                else
                    obj.ParsedFrameIndices(tindex, 2) = obj.ParsedFrameIndices(tindex, 1) + NumFramesAfter(tindex);
                end

            end
        end %parse
        
        function Data = get.Data(obj)
            numFramesPerTrial = diff(obj.ParsedFrameIndices, 2);
            numTrials = size(obj.ParsedFrameIndices, 1);
            if range(numFramesPerTrial) == 0
                Data = zeros(numTrials, numFramesPerTrial(1));
                for tindex = 1:numTrials
                    Data(tindex, :) = obj.RawData(obj.ParsedFrameIndices(tindex,1):obj.ParsedFrameIndices(tindex,2));
                end
            else
                Data = cell(numTrials, 1);
                for tindex = 1:numTrials
                    Data{tindex} = obj.RawData(obj.ParsedFrameIndices(tindex,1):obj.ParsedFrameIndices(tindex,2));
                end
            end
        end
        
        function Neuropil = get.Neuropil(obj)
            numFramesPerTrial = diff(obj.ParsedFrameIndices, 2);
            numTrials = size(obj.ParsedFrameIndices, 1);
            if range(numFramesPerTrial) == 0
                Neuropil = zeros(numTrials, numFramesPerTrial(1));
                for tindex = 1:numTrials
                    Neuropil(tindex, :) = obj.RawNeuropil(obj.ParsedFrameIndices(tindex,1):obj.ParsedFrameIndices(tindex,2));
                end
            else
                Neuropil = cell(numTrials, 1);
                for tindex = 1:numTrials
                    Neuropil{tindex} = obj.RawNeuropil(obj.ParsedFrameIndices(tindex,1):obj.ParsedFrameIndices(tindex,2));
                end
            end
        end
        
%         function Stimulus = get.Stimulus(obj)
%             numFramesPerTrial = diff(obj.ParsedFrameIndices, 2);
%             numTrials = size(obj.ParsedFrameIndices, 1);
%             if range(numFramesPerTrial) == 0
%                 Stimulus = zeros(numTrials, numFramesPerTrial(1));
%                 for tindex = 1:numTrials
%                     Stimulus(tindex, obj.ParsedStimIndices(1):obj.ParsedStimIndices(2)) = 1;
%                 end
%             else
%                 Stimulus = cell(numTrials, 1);
%                 for tindex = 1:numTrials
%                     Stimulus{tindex} = zeros(1,diff(obj.ParsedFrameIndices(tindex,1))+1);
%                     Stimulus{tindex}(obj.ParsedStimIndices(1):obj.ParsedStimIndices(2)) = 1;
%                 end
%             end
%         end
                
    end %methods
end %classdef
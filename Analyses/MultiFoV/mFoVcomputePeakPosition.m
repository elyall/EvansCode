function [StimEvokedMax, MaxDFoF, Data, Map] = mFoVcomputePeakPosition(Data, Map, varargin)

blur = false;
filt = fspecial('gaussian', 5, 1);
mergetype = 'pretty'; % 'quick' or 'pretty'
% Crop = false;
Crop = repmat([32.51, 0, 729.98, 512], numel(Data), 1);

% Default settings
ExperimentFiles = [];

%% Parse input arguments
if ~exist('Data','var') || isempty(Data)
    directory = CanalSettings('ExperimentDirectory');
    [Data, p] = uigetfile({'*.mat'},'Choose Experiment file',directory,'MultiSelect','on');
    if isnumeric(Data)
        return
    end
    if iscellstr(Data)
        for findex = 1:numel(Data)
            Data{findex} = fullfile(p,Data{findex});
        end
    else
        error('Use ''vidAverageStim'' instead for a single FoV');
    end
end

if ~exist('Map', 'var')
    Map = [];
end

index = 1;
while index<=length(varargin)
    switch varargin{index}
        case 'blur'
            blur = true;
            index = index + 1;
        case 'filter'
            filter = varargin{index+1};
            index = index + 2;
        case 'Crop'
            Crop = varargin{index + 1};
            index = index + 2;
        case 'files'
            ExperimentFiles = varargin{index+1};
            index = index + 2;
        otherwise
            warning('Argument ''%s'' not recognized',varargin{index});
            index = index + 1;
    end
end

%% Load data
numFiles = numel(Data);
fprintf('mFoV Computing Peak Position\n');
if iscellstr(Data)
    fprintf(' \tloading %d files...', numFiles);
    ExperimentFiles = Data;
    clear Data;
    for findex = 1:numFiles
        Data(findex) = load(ExperimentFiles{findex}, 'AnalysisInfo', 'AvgTrialdFoF', 'Map', '-mat');
        fprintf('\n\tloaded %d: %s', findex, ExperimentFiles{findex});
    end
    for findex = 1:numFiles
        Data(findex).filename = ExperimentFiles{findex};
        Data(findex).origMap = Data(findex).Map;
        Data(findex).cropped = false;
    end
elseif isstruct(Data)
    if ~isfield(Data, 'AvgTrialdFoF') || ~isfield(Data, 'AnalysisInfo')
        if isempty(ExperimentFiles)
            error('First argument is a ''Data'' struct, but does not contain necessary fields. Please pass in the ''ExperimentFiles'' as an argument');
        else
            for findex = 1:numFiles
                load(ExperimentFiles{findex}, 'AnalysisInfo', 'AvgTrialdFoF', '-mat');
                Data(findex).AnalysisInfo = AnalysisInfo;
                Data(findex).AvgTrialdFoF = AvgTrialdFoF;
                fprintf('\n\tloaded %d: %s', findex, ExperimentFiles{findex});
            end
        end
    end
end

%% Determine stimulus periods
numStims = numel(Data(1).AvgTrialdFoF);
stimIndex = zeros(numStims, 2);
for sindex = 1:numStims
    stimIndex(sindex, :) = mode(Data(1).AnalysisInfo.TrialStimFrames(Data(1).AnalysisInfo.StimID==sindex-1,:), 1);
end

%% Initialize Map for Generating Images
temp = cell(numFiles,1);
if isempty(Map) && numFiles > 1
    for findex = 1:numFiles
        temp{findex} = Data(findex).AvgTrialdFoF{4}(:,:,:,30);
    end
    [MaxDFoF,~,Map,Data] = createMultiFoVImage(Data, temp, mergetype, Crop);
    fprintf('Completed map\n');
%     assignin('base', 'Map', Map);
%     assignin('base', 'Data', Data);
elseif ~isempty(Map)
    sz = size(Map);
    MaxDFoF = nan(sz(1), sz(2));
else
    sz = size(Data(1).AvgTrialdFoF{1});
    MaxDFoF = nan(sz(1), sz(2));
end

%% Determine peak response
StimEvokedMax = MaxDFoF;
MaxDFoF(:) = -inf;
for sindex = 1:numStims
    
    for findex = stimIndex(sindex, 1):stimIndex(sindex, 2)
        
        % Build whole frame
        for fileindex = 1:numFiles
            temp{fileindex} = Data(fileindex).AvgTrialdFoF{sindex}(:,:,:,findex);
            if blur
                temp{fileindex} = imfilter(temp{fileindex}, filt);
            end
        end
        if numFiles > 1
            Image = createMultiFoVImage(Data, temp, 'pretty', Crop, Map);
        else
            Image = temp{1};
        end
        
        % Determine if value in frame is larger than value stored
        [MaxDFoF, pos] = max(cat(3,MaxDFoF,Image),[],3);
        StimEvokedMax(pos==2) = sindex;
        
    end

end

function Files = splitIntoTrials(File,AnalysisInfo,varargin)
% Splits a single sbx file into one file for each trial. NUMFRAMESBEFORE
% sets how many frames is saved before each stimulus, while the number of
% frames saved after each stimulus extends until the NUMFRAMESBEFORE cut
% off for the next stimulus.

numFramesBefore = 16;
SaveBase = '';          % beginning of filename to save all files as
Depth = [1,inf];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numFramesBefore'
                numFramesBefore = varargin{index+1};
                index = index + 2;
            case 'SaveBase'
                SaveBase = varargin{index+1};
                index = index + 2;
            case 'Depth'
                Depth = varargin{index+1};
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

if isempty(SaveBase)
    SaveBase = File(1:end-4);
end
if ischar(AnalysisInfo)
    load(AnalysisInfo,'AnalysisInfo','-mat');
end


%% Determine frame indices to save for each trial
FrameIndex = AnalysisInfo.ExpStimFrames(:,1)-numFramesBefore;                         % set beginning of each trial based on NUMFRAMESBEFORE
FrameIndex = [FrameIndex,[FrameIndex(2:end)-1;0]];                                    % extend each trial up to first frame in next trial
Config = load2PConfig(File); % load config information
FrameIndex(end,2) = min(FrameIndex(end,1)+mode(diff(FrameIndex,[],2)),Config.Frames); % set last trial to have as many frames as the other trials (or extend to end of the file if that's shorter)
FrameIndex = cat(1,[2,FrameIndex(1,1)-1],FrameIndex);                                 % add init frames before 1st trial while ignoring 1st frame as it's incomplete
N = size(FrameIndex,1); % # of trials


%% Determine metadata for each file

% Load header
info = Config.header{1};
frame = info.frame;
line = info.line;
event_id = info.event_id;

% Determine which trigger occurs in which file
inds = cat(3,bsxfun(@ge,frame',FrameIndex(:,1)),bsxfun(@le,frame',FrameIndex(:,2))); % find FrameIndices above and below each trigger
inds = all(inds,3);         % find specific trial that encompasses each trigger
[inds,~] = find(inds);      % determine index of each trial

% Create header for each file
info.frame = [];
info.line = [];
info.event_id = [];
info = repmat({info},N,1);
for t = 1:N
    info{t}.frame = frame(inds==t)-FrameIndex(t,1); % offset due to chunking (doesn't take into account depth)
    info{t}.line = line(inds==t);
    info{t}.event_id = event_id(inds==t);
end


%% Save trials to individual files
Files = cell(N,1);
for t = 1:N
    imgs = load2P(File,'Frames',FrameIndex(t,1):FrameIndex(t,2),'Depth',Depth); % load data
    fn = sprintf('%s_%03d.sbx',SaveBase,t-1);                                   % create filename
    Files{t} = save2P(fn,imgs,'Header',info{t});                                % save data to file
end


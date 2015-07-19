function info = fixTimeStamps(InfoFile)
% assumes 6 second trials with:
% start bar moving pulse
% 1 frame later end bar moving pulse
% 22 frames later start bar moving pulse
% 1 frame later end bar moving pulse
% 70 frames later, repeat
correctDiff = [1,22,1,70]';
saveOut = false;

%% Check input arguments
narginchk(0, 1)
if ~exist('InfoFile', 'var') || isempty(InfoFile)
    directory = CanalSettings('DataDirectory');
    [InfoFile, p] = uigetfile({'*.mat'}, 'Choose Info file for sbx file', directory);
    if isnumeric(InfoFile)
        return
    end
    InfoFile = fullfile(p, InfoFile);
    load(InfoFile, 'info');
elseif isstruct(InfoFile)
    info = InfoFile;
elseif isdir(InfoFile)
    [InfoFile, p] = uigetfile({'*.mat'}, 'Choose Info file for sbx file', InfoFile);
    if isnumeric(InfoFile)
        return
    end
    InfoFile = fullfile(p,InfoFile);
    load(InfoFile, 'info');
elseif ischar(InfoFile)
    load(InfoFile, 'info');
end

%% Save original info
if saveOut
    vars = whos(matfile(InfoFile));
    if ~any(strcmp('origInfo', {vars.name})) % first time doing it
        origInfo = info;
        save(InfoFile, 'origInfo', '-append');
    else % load original info for fixing time stamps
        load(InfoFile, 'origInfo', '-mat');
        info = origInfo;
    end
end

%% Throw out all stamps not on first port
bad = info.event_id~=1;
info.event_id(bad) = [];
info.line(bad) = [];
info.frame(bad) = [];

%% Move to first valid start trigger
framesBetween = diff(info.frame);
while sum(abs(framesBetween(1:2)-correctDiff(1:2))) > 4 %assumes first trigger is start of trial
    framesBetween(1) = [];
    info.frame(1) = [];
    info.event_id(1) = [];
    info.line(1) = [];
end

%% Fix time stamps
index = 1;
go = false;
while ~go
    fprintf('\n%d\t%d', framesBetween(index), correctDiff(rem(index-1,4)+1));
    if abs(framesBetween(index)-correctDiff(rem(index-1,4)+1)) > 5 %more than 5 frames apart
        fprintf('\t fixed');     
        if info.frame(index+1)-info.frame(index) > 5*max(correctDiff)
            info.frame(index+1) = []; % remove it
            info.line(index+1) = [];
            info.event_id(index+1) = [];
            index = index - 1;
        else
            info.frame(index+2:end+1) = info.frame(index+1:end); % shift back one
            info.frame(index+1) = info.frame(index)+correctDiff(rem(index-1,4)+1); % replace
            info.event_id(end+1) = 1;
            info.line(index+2:end+1) = info.line(index+1:end); % shift back one
            info.line(index+1) = NaN;
        end
        framesBetween = diff(info.frame);
    end
    index = index + 1;
    if index > numel(framesBetween)
        go = true;
    end
end

%% Add frames to end
while rem(numel(info.frame),4) ~= 0
    info.frame(index+1) = info.frame(index)+correctDiff(rem(index-1,4)+1); % replace
    info.event_id(end+1) = 1;
    info.line(index+1) = NaN;
    index = index + 1;
end

%% Save to file
if saveOut
    save(InfoFile, 'info', '-append');
end

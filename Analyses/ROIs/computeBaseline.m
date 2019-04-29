function Baseline = computeBaseline(Data,AnalysisInfo,varargin)
% Data is numFrames x numROIs
% Stimulus is vector of length numFrames

type = 'catch'; % 'catch' or 'movprctile'
numFrames = [];

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case 'numFrames'
                numFrames = varargin{index+1};
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

if isempty(numFrames)
    switch type
        case 'catch'
            numFrames = 4;
        case 'movprctile'
            numFrames = 3001;
    end
end

%%
switch type
    case 'catch'
        %% catch
        
        Index = AnalysisInfo.StimID==0;
        Index = [0;Index(1:end-1)];
        frames = AnalysisInfo.ExpStimFrames(Index,1)-1;
        frames = [max(frames-numFrames,1),frames];
        
        
        numStim = size(frames,1);
        Frames = nan(numFrames,numROIs,numStim);
        for s = 1:numStim
            Frames(:,:,s) = Data(frames(s,1):frames(s,2),:);
        end
        Baseline = squeeze(mean(Frames));
        
    case 'movprctile'
        %% movprctile
        Baseline = movprctile(Data,5,numFrames,1); % SLOW
        
end




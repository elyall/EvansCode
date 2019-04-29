function Baseline = computeBaseline(Data,AnalysisInfo,varargin)
% Data is numFrames x numROIs
% Stimulus is vector of length numFrames
% Baseline is 1 x numROIs (catch) or numFrames x numROIs (movprctile)

type = 'catch'; % 'catch' or 'movprctile'
numFrames = [];
prc = 5;

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
        
        Index = AnalysisInfo.StimID==0; % index of catch stimuli
        Index = logical([0;Index(1:end-1)]); % index of stimuli after catch
        frames = AnalysisInfo.ExpStimFrames(Index,1)-1;
        frames = [max(frames-numFrames+1,1),frames];
        
        
        numStim = size(frames,1);
        numROIs = size(Data,2);
        Frames = nan(numFrames,numROIs,numStim);
        for s = 1:numStim
            Frames(:,:,s) = Data(frames(s,1):frames(s,2),:); % gather last second after catch before next stim
        end
        Baseline = squeeze(nanmean(nanmean(Frames,1),3)); % average over time and over trials
        
        
    case 'movprctile'
        %% movprctile
        Baseline = movprctile(Data,prc,numFrames,1); % SLOW
        
end




function plotIndividualTraces(ROIFile, ROIid)

saveToPDF = false;

% Data to display
datatype = 'dFoF'; % 'dFoF' or 'raw'
sorttype = 'stim'; % 'stim' or 'trial' (will not show averages)
showdata = true;
showavgs = true;

% Computation parameters
minrunspeed = 0; % (will only run if 'RunningSpeed' found in ROIFile)
numSTDsOutlier = 6; % Inf means no trials are thrown out
spikeThreshold = 0.01; % set to "0" to not overlay spikes

% Display parameters
zscoreby = 'all'; %'all' or 'trial' or 'none'
desiredOffset = 1;
avgPixelHeight = 6; % number of pixels of average plots height

% Constant parameters
controlindex = 1; %index of control stimulus ID relative to other StimIDs (ID is minimum value => 1)


%% Check input arguments
narginchk(0, 2);
if ~exist('ROIFile','var') || isempty(ROIFile)
    directory = CanalSettings('DataDirectory');
    [ROIFile, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIFile)
        return
    end
    ROIFile = fullfile(p,ROIFile);
end
if ~exist('ROIid','var') || isempty(ROIid)
    ROIid = 'all';
end

if ~exist('Fs','var') || isempty(Fs)
    Fs=1;
    xlab='Frames';
else
    xlab='Time (s)';
end
if ~exist('framelimit', 'var')
    framelimit='all';
end

%% Load in data
load(ROIFile, 'ROIdata', 'AnalysisInfo', 'RunningSpeed', '-mat');
if ~exist('RunningSpeed', 'var') && minrunspeed ~= 0
    minrunspeed = 0;
end

 if spikeThreshold > 0 && ~isfield(ROIdata.rois, 'spikes');
     spikeThreshold = 0;
 end
 
%% Plot rasters
numROIs = numel(ROIdata.rois);
for rindex = 1:numROIs %cycle through rois
    if (ischar(ROIid) && strcmp(ROIid, 'all')) || ismember(rindex,ROIid) %only display for ROIs selected
        
        % Sort stimuli
        [StimIDs, ~, TrialIndex] = unique(AnalysisInfo.StimID);
        numStimuli = numel(StimIDs);
        numTrials = numel(TrialIndex);
        switch sorttype
            case 'stim'
                TrialIndex(TrialIndex == controlindex) = -Inf; %doesn't matter if controlindex is 1 (but does if controlindex is higher)
                [TrialIndex, displayIndex] = sort(TrialIndex);
                TrialIndex(TrialIndex == -Inf) = controlindex;
                blockLength = [];
                for sindex = 1:numStimuli
                    blockLength = cat(1,blockLength,sum(TrialIndex==sindex));
                end
                labels = [{'control'}; cellstr(num2str((1:numStimuli-1)'))];
            case 'trial'
                displayIndex = 1:numTrials;
                showavgs = false;
        end
        
        % Sort data
        switch datatype
            case 'dFoF'
                data = ROIdata.rois(rindex).dFoF(displayIndex, :)-1;
            case 'raw'
                data = ROIdata.rois(rindex).data(displayIndex, :);
        end
        if spikeThreshold > 0
            spikes = ROIdata.rois(rindex).spikes(displayIndex,:) > spikeThreshold;
        end
        StimulusFrames = cat(2, AnalysisInfo.TrialStimFrames(displayIndex, :), zeros(numTrials,1));
        StimulusFrames(TrialIndex==controlindex,1:2) = nan; % remove stim frames for control trials
        
        % Removing non-running trials
        if minrunspeed
            runspeed = RunningSpeed(displayIndex, :);
            StimPeriodAvgSpeed = mean(runspeed(:, StimulusFrames(:,1):StimulusFrames(:,2)),2);
            NotRunning = StimPeriodAvgSpeed <= minrunspeed;
            data(NotRunning, :) = [];
            StimulusFrames(NotRunning, :) = [];
            TrialIndex(NotRunning, :) = [];
            if spikeThreshold > 0
                spikes(NotRunning, :) = [];
            end
        end
        
        % Remove outliers for each stimulus
        if strcmp(sorttype, 'stim') && strcmp(zscoreby, 'all')
            for sindex = 1:numStimuli
                n = inf;
                currentData = data(TrialIndex==sindex,:);
                [H,W] = size(currentData);
                currentData = reshape(zscore(currentData(:)), H, W);
                while n ~= blockLength(sindex)
                    n = blockLength(sindex);
                    bad = find(max(abs(currentData),[],2) > numSTDsOutlier); %throw out outlier trials
                    if ~isempty(bad)
                        currentData(bad, :) = [];
                        blockLength(sindex) = blockLength(sindex) - numel(bad);
                        bad = bad + sum(blockLength(1:sindex-1));
                        data(bad, :) = [];
                        StimulusFrames(bad, :) = [];
                        TrialIndex(bad, :) = [];
                        if spikeThreshold > 0
                            spikes(bad, :) = [];
                        end
                    end
                end
            end
        end
        
        % Compute average stimuli
        AvgStim = zeros(numStimuli, size(data, 2));
        AvgStimuliFrames = [];
        for sindex = 1:numStimuli
            AvgStim(sindex, :) = mean(data(TrialIndex==sindex, :), 1);
            AvgStimuliFrames = cat(1, AvgStimuliFrames, [mode(StimulusFrames(TrialIndex==sindex, :), 1),0]);
        end
        
        % Z-score data
        switch zscoreby
            case 'trial'
                data = zscore(data,0,2); %zscore each row individually
                AvgStim = zscore(AvgStim,0,2);
                ylab = 'Trials: z-score (per trial)';
                data = bsxfun(@rdivide, data, max(abs(data), [], 2)); % Normalize each trace to between 0 to 1
                AvgStim = bsxfun(@rdivide, AvgStim, max(abs(data), [], 2))*avgPixelHeight;
            case 'all'
                [H,W] = size(data);
                data = reshape(zscore(data(:)),H,W);
                [H,W] = size(AvgStim);
                AvgStim = reshape(zscore(AvgStim(:)),H,W);
                ylab = 'Trials: z-score (all)';
                data = data/max(abs(data(:))); % Normalize data to between 0 to 1
                AvgStim = AvgStim/max(abs(AvgStim(:)))*avgPixelHeight;
            case 'none'
                % do nothing
                data = bsxfun(@rdivide, data, max(abs(data), [], 2)); % Normalize each trace to between 0 to 1
                AvgStim = bsxfun(@rdivide, AvgStim, max(abs(AvgStim), [], 2))*avgPixelHeight;
                switch datatype
                    case 'raw'
                        ylab = 'Trials: Fluorescence (AU)';
                    case 'dFoF'
                        ylab = 'Trials: dF/F';
                end
        end
        x = repmat(0:1/Fs:(size(data,2)-1)/Fs, size(data,1), 1);
        
        % Crop data
        if isnumeric(framelimit) && numel(framelimit)==2
            xshift = -framelimit(1)+1;
            data = data(:,framelimit(1):framelimit(2));
            AvgStim = AvgStim(:,framelimit(1):framelimit(2));
        else
            xshift = 0;
        end
        
        % Offset traces
        offset = desiredOffset*(size(data,1)-1):-desiredOffset:0;
        data = bsxfun(@plus, data, offset');
        if showavgs
            offset = avgPixelHeight*(numStimuli-1):-avgPixelHeight:0;
            AvgStim = bsxfun(@plus, AvgStim, offset');
            AvgStim_x = repmat(0:1/Fs:(size(AvgStim,2)-1)/Fs, numStimuli, 1);
            data = data + avgPixelHeight*numStimuli;
        end
        
        % Add average stimuli to matrix to display
        if showavgs && showdata
            blockLength = cat(1, blockLength, repmat(avgPixelHeight, numStimuli,1));
            labels = cat(1, labels, [{'control avg'}; cellstr([num2str((1:numStimuli-1)'),repmat(' avg',numStimuli-1,1)])]);
        elseif showavgs
            blockLength = repmat(avgPixelHeight, numStimuli,1);
            labels = [{'control avg'}; cellstr([num2str((1:numStimuli-1)'),repmat(' avg',numStimuli-1,1)])];
        end
        
        % Display Data
        hf = figure();
        hold on
        if showdata
            plot(x',data');
        end
        if showavgs
            plot(AvgStim_x', AvgStim');
        end
        set(gca,'TickDir','out')
        xlabel(xlab);
        ylabel(ylab);
        axis tight
        
        % Plot Stimulation Shading
%         if exist('stimindex','var') && ~isempty(stimindex)
%             maxy=get(gca,'YLim');
%             % find first and last frame for each stimulus
%             stims=[0;stimindex]-[stimindex;0];
%             firsts=find(stims<0)/Fs-1/Fs;
%             lasts=(find(stims>0)-1)/Fs-1/Fs;
%             % plot stimulus marks
%             hold on
%             t=size(d,2);
%             for i=1:length(firsts);
%                 if firsts(i)>t
%                     return
%                 elseif lasts(i)>t
%                     lasts(i)=t;
%                 end
%                 h=area([firsts(i),lasts(i)],[maxy(2),maxy(2)],-.5);
%                 child=get(h,'Children');
%                 set(child,'FaceAlpha',0.5,'EdgeColor','none');
%             end
%         end
        
        % Overlay spikes
        if spikeThreshold > 0
            spikex = x(spikes);
            spikey = data(spikes);
            plot(spikex, spikey, 'r*');
        end
        
        % Plot block lines
        if strcmp(sorttype, 'stim')
            blockLines = cumsum(flipud(blockLength))-1;
            blockLines = blockLines(1:end-1);
            if showavgs
                blockLines(1:numStimuli-1) = [];
            end
            plot(repmat([0,size(data,2)-1],length(blockLines),1)',repmat(blockLines+desiredOffset/2,1,2)','k-','LineWidth',2); % trial dividing lines
        end
        
        % Plot stim lines
%         plot(repmat([StimulusFrames(:,1)-.5, StimulusFrames(:,2)+.5]+xshift,2,1), repmat((0.5:1:size(StimulusFrames,1))',2,2), 'k--','LineWidth',2);
        
        % Label y-axis
        ticks = cumsum(flipud(blockLength))-diff(cumsum([0;flipud(blockLength)]))/2 + 0.5;
        set(gca,'YTick',ticks, 'YTickLabel',flipud(labels));
        ylabel('Stimulus');
        
        % Label x-axis
%         set(gca,'XTick',2*Fs:2*Fs:framelimit(2), 'XTickLabel',num2str((2:2:framelimit(2)/Fs)'))
%         xlabel('Frame');
              
        % Save plot to PDF
        if saveToPDF
            % Fix label
            if isempty(ROIdata.rois(rindex).label)
                ROIdata.rois(rindex).label = {'unlabeled'};
            end
            title(sprintf('ROI: %d, Label: %s', rindex, ROIdata.rois(rindex).label{1}));
            export_fig(hf, 'test.pdf','-append'); close(hf);
        end
        
    end %good roi flag
end %roi cycle


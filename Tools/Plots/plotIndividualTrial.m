function colors = plotIndividualTrial(ROIdata, ROIid, TrialID, StimID, varargin)

type = 'traces'; % 'raster' or 'traces' or 'spikes'
datatype = 'dFoF'; % 'dFoF' or 'raw'
zscoreby = 'none'; % 'all' or 'individual' or 'none'

desiredOffset = 4; % for 'traces' plotting
FrameRate = 15.45;

showSpikes = false;
spikeThresh = .18;

% spikes-only
colors = jet(10);
colors = colors(randperm(10),:);

hA = []; %axes handle
Title = '';

%% Check input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    directory = CanalSettings('DataDirectory');
    [ROIFile, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIFile)
        return
    end
    ROIFile = fullfile(p,ROIFile);
    load(ROIFile, 'ROIdata');
elseif ischar(ROIdata)
    load(ROIdata, 'ROIdata');
end

if ~exist('ROIid','var') || isempty(ROIid)
    ROIid = 'all';
end

if ~exist('TrialID','var') || isempty(TrialID)
    TrialID = 1;
elseif ischar(TrialID)
    TrialID = 1:numel(ROIdata.DataInfo.TrialIndex);
end

if ~exist('StimID','var')
    StimID = [];
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'datatype'
                datatype = varargin{index+1};
                index = index + 2;
            case 'zscoreby'
                zscoreby = varargin{index+1};
                index = index + 2;
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case 'desiredOffset'
                desiredOffset = varargin{index+1};
                index = index + 2;
            case 'FrameRate'
                FrameRate = varargin{index+1};
                index = index + 2;
            case 'showData'
                showData = varargin{index+1};
                index = index + 2;
            case 'showSpikes'
                showSpikes = varargin{index+1};
                index = index + 2;
            case 'spikeThresh'
                spikeThresh = varargin{index+1};
                index = index + 2;
            case 'colors'
                colors = varargin{index+1};
                index = index + 2;
            case 'Title'
                Title = varargin{index+1};
                index = index + 2;
            case 'axes'
                hA = varargin{index+1};
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

if strcmp(type, 'spikes')
    showSpikes = true;
end

%% Determine ROIs
if ischar(ROIid) && strcmp(ROIid, 'all')
    ROIid = 1:numel(ROIdata.rois);
end
numROIs = numel(ROIid);

%% Create axes
if isempty(hA)
    for tindex = 1:numel(TrialID)
        figure();
        hA(tindex) = axes();
    end
end
nColors = size(colors,1);

%% Display each trial
numFrames = ROIdata.DataInfo.numFramesBefore + ROIdata.DataInfo.numFramesAfter + 1;
for tindex = 1:numel(TrialID)
    
    % Determine trial info
    if isempty(StimID)
        TrialIndex = TrialID(tindex);
    else
        TrialIndex = find(ROIdata.DataInfo.StimID == StimID);
        TrialIndex = TrialIndex(TrialID(tindex));
    end
    StimIndex = ROIdata.DataInfo.StimID(TrialIndex);
    StimFrames = [ROIdata.DataInfo.numFramesBefore + 1, ROIdata.DataInfo.numFramesBefore + ROIdata.DataInfo.numStimFrames(TrialIndex)];
    
    % Extract data
    if any(strcmp(type, {'raster', 'traces'}))
        data = zeros(numROIs, numFrames);
    end
    if showSpikes
        spikes = zeros(numROIs, numFrames);
    end
    for rindex = 1:numROIs
        if any(strcmp(type, {'raster', 'traces'}))
            switch datatype
                case 'dFoF'
                    data(rindex,:) = ROIdata.rois(ROIid(rindex)).dFoF(TrialIndex,:);
                case 'raw'
                    data(rindex,:) = ROIdata.rois(ROIid(rindex)).data(TrialIndex,:);
            end
        end
        if showSpikes
            spikes(rindex, :) = ROIdata.rois(ROIid(rindex)).spikes(TrialIndex,:);
        end
    end
    
    % Z-score data
    if any(strcmp(type, {'raster', 'traces'}))
        switch zscoreby
            case 'all'
                data = reshape(zscore(data(:)),numROIs,numFrames);
                switch datatype
                    case 'raw'
                        colorbarTitle = 'z-scored fluorescence';
                    case 'dFoF'
                        colorbarTitle = 'z-scored dF/F';
                end
            case 'individual'
                data = zscore(data, [], 2);
                switch datatype
                    case 'raw'
                        colorbarTitle = 'z-scored individual fluorescence';
                    case 'dFoF'
                        colorbarTitle = 'z-scored individual dF/F';
                end
            case 'none'
                switch datatype
                    case 'raw'
                        colorbarTitle = 'Fluorescence (AU)';
                    case 'dFoF'
                        colorbarTitle = 'dF/F';
                end
        end
    end
    
    % Display data
    axes(hA(tindex));
    switch type
        case 'raster'
            imagesc(data);
            temp = get(gca,'CLim');
            cmap = b2r(temp(1),temp(2));
            colormap(cmap);
            hold on
            
            % spikes
            if showSpikes
                [Rindex, Tindex] = find(spikes>=spikeThresh);
                plot(Tindex, Rindex, 'g*');
            end
            
            % stim lines
            StimFrames(1) = StimFrames(1)-0.5;
            StimFrames(2) = StimFrames(2)+0.5;
            
            % axes
            set(gca, 'YTick', 1:10:numROIs, 'YTickLabel', cellstr(num2str((1:10:numROIs)')));
            set(gca, 'XTick', 0.5:FrameRate:numFrames+0.5, 'XTickLabel', cellstr(num2str((0:numFrames/FrameRate)')));
            
            % colorbar
            hc = colorbar;
            ylabel(hc, colorbarTitle);
            colors = colormap;
            
        case 'traces'
            x = repmat(0:1/FrameRate:(numFrames-1)/FrameRate, numROIs, 1);
            offset = desiredOffset*(numROIs-1):-desiredOffset:0;
            data = bsxfun(@plus, data, offset');
            plot(x',data');
            axis tight
            hold on
            
            % spikes
            if showSpikes
                [Rindex, Tindex] = find(spikes>=spikeThresh);
                Rindex = sub2ind([numROIs, numFrames], Rindex, Tindex);
                plot(x(1,Tindex), data(Rindex), 'k.');
            end
            
            % stim lines
            StimFrames(1) = (StimFrames(1)-0.5)/FrameRate;
            StimFrames(2) = (StimFrames(2)+0.5)/FrameRate;
            
            % axes
            if strcmp(zscoreby, 'none') && strcmp(datatype, 'dFoF')
                set(gca, 'YTick', 1:desiredOffset:desiredOffset*(numROIs-1)+1, 'YTickLabel', cellstr(num2str((numROIs:-1:1)')));
            else
                set(gca, 'YTick', 0:desiredOffset:desiredOffset*(numROIs-1), 'YTickLabel', cellstr(num2str((numROIs:-1:1)')));
            end
            set(gca, 'XTick', 0:numFrames/FrameRate, 'XTickLabel', cellstr(num2str((0:numFrames/FrameRate)')));
            colors = get(gca, 'ColorOrder');
            
        case 'spikes'
            
            % stim lines
            plot(repmat(StimFrames,2,1), repmat([0.5,numROIs+0.5],2,1)', 'k--','LineWidth',2);
            
            % raster
            [Rindex, Tindex] = find(spikes>=spikeThresh);
            Rindex = numROIs-Rindex+1;
            YLocation = [Rindex-0.4, Rindex+0.4];
            Tindex = repmat(Tindex, 1, 2);
            hold on
            for rindex = 1:numROIs
                index = Rindex==rindex;
                plot(Tindex(index,:)', YLocation(index, :)','-','LineWidth',2,'Color',colors(rem(numROIs-rindex,nColors)+1,:));
            end
            ylim([0.5,numROIs+0.5]);
            
            % axes
            set(gca, 'YTick', 1:numROIs, 'YTickLabel', cellstr(num2str((numROIs:-1:1)')));
            set(gca, 'XTick', 0:FrameRate:numFrames, 'XTickLabel', cellstr(num2str((0:numFrames/FrameRate)')));
    end
    set(gca,'TickDir','out')
    
    % stim lines
    if any(strcmp(type, {'raster', 'traces'}))
        YLim = get(gca, 'YLim');
        plot(repmat(StimFrames,2,1), repmat(YLim,2,1)', 'k--','LineWidth',2);
    end
    
    % labels
    if isequal(Title, true) || isempty(Title)
        title(sprintf('Trial: %d (Stimulus: %d)', TrialIndex, StimIndex));
    elseif ischar(Title)
        title(Title);
    end
    ylabel('ROI');
    xlabel('Time (s)');
      
end
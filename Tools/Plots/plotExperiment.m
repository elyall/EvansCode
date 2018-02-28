function [hA,Data,NormalizeBy,SubtractOff] = plotExperiment(Data, Neuropil, varargin)

% Data
ROIindex = [1 inf];
FrameIndex = [1 inf];
ExpStimFrames = [];

% Computations
NeuropilWeight = false;
baselinePrctile = 0;
adjustment = 'normalize'; % '', 'zscore', 'normalize'
NormalizeBy = [];
SubtractOff = [];
smoothType = ''; % '', 'moving', 'lowess', 'loess', 'sgolay', 'rlowess', or 'rloess'

% Display
Type = '2Dlines'; % '2Dlines', '3Dlines', or 'image'
LineWidth = 2;
spacing = 1;
XLabel = 'Time (sec)';
frameRate = 15.45/4;
Colors = [];

hA = [];
directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Type'
                Type = varargin{index+1};
                index = index + 2;
            case 'LineWidth'
                LineWidth = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
                index = index + 2;
            case 'ExpStimFrames'
                ExpStimFrames = varargin{index+1};
                index = index + 2;
            case 'NeuropilWeight'
                NeuropilWeight = varargin{index+1};
                index = index + 2;
            case {'Colors','Color'}
                Colors = varargin{index+1};
                index = index + 2;
            case 'adjustment'
                adjustment = varargin{index+1};
                index = index + 2;
            case 'smoothType'
                smoothType = varargin{index+1};
                index = index + 2;
            case 'NormalizeBy'
                NormalizeBy = varargin{index+1};
                index = index + 2;
            case 'SubtractOff'
                SubtractOff = varargin{index+1};
                index = index + 2;
            case 'hA'
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

if ~exist('Data','var') || isempty(Data)
    [Data, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory);
    if isnumeric(Data)
        return
    end
    Data = fullfile(p,Data);
end


%% Load ROIs
if ischar(Data)
    ROIFile = Data;
    load(ROIFile, 'ROIdata', '-mat');
    Data = ROIdata;
end
if isstruct(Data)
    totalROIs = numel(Data.rois);
    totalFrames = numel(Data.rois(1).rawdata);
else
    [totalFrames,totalROIs] = size(Data);
end


%% Determine data to analyze

% Determine ROIs to display
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(1:end-1)+1:totalROIs];
end
numROIs = numel(ROIindex);

% Convert ROIdata to matrices
if isstruct(Data)
    Neuropil = reshape([Data.rois(ROIindex).rawneuropil], totalFrames, numROIs);
    Data = reshape([Data.rois(ROIindex).rawdata], totalFrames, numROIs);
else
    Data = Data(:,ROIindex);
    if exist('Neuropil','var') && ~isempty(Neuropil)
        Neuropil = Neuropil(:,ROIindex);
    end
end

% Determine frames to pull from
if FrameIndex(end) == inf
    FrameIndex = [FrameIndex(1:end-1), FrameIndex(1:end-1)+1:totalFrames];
end

% Format experiment data
if ~isempty(ExpStimFrames)
    if isvector(ExpStimFrames)
        Stim = ExpStimFrames;
    else
        Stim = zeros(1,totalFrames);
        for tindex = 1:size(ExpStimFrames,1)
            Stim(ExpStimFrames(tindex,1):ExpStimFrames(tindex,2)) = 1;
        end
    end
    Stim = Stim(FrameIndex);
% ExpStimFrames = ExpStimFrames - FrameIndex(1) + 1;
% ExpStimFrames(all(ExpStimFrames<0,2),:) = [];
% ExpStimFrames(all(ExpStimFrames>FrameIndex(end),2),:) = [];
% if ExpStimFrames(1)<0
%     ExpStimFrames(1) = 0;
% end
% if ExpStimFrames(end) > FrameIndex(end)
%     ExpStimFrames(end) = FrameIndex(end);
% end
end

% Determine colors of lines
if ismember(Type,{'2Dlines','3Dlines'}) && ~isempty(Colors) && size(Colors,1) ~= numROIs
    Colors = repmat(Colors,ceil(numROIs/size(Colors,1)),1);
end


%% Subtract off neuropil
if NeuropilWeight
    if isequal(NeuropilWeight,true)
        NeuropilWeight = determineNeuropilWeight(Data, Neuropil);
    end
    Data = Data - bsxfun(@times,NeuropilWeight',Neuropil);
end


%% Compute dFoF
if baselinePrctile
    baseline = prctile(Data, baselinePrctile);
    Data = bsxfun(@rdivide, bsxfun(@minus, Data, baseline), baseline);
end


%% Smooth data
if ~isempty(smoothType)
    Data = smooth(Data, smoothType);
end


%% Adjust data
Data = Data(FrameIndex,:);
if ~isempty(adjustment)
    switch adjustment
        case 'zscore'
            Data = zscore(Data);
        case 'normalize'
            if isempty(SubtractOff)
                SubtractOff = min(Data); % subtract off minimum
            end
            Data = bsxfun(@minus, Data, SubtractOff);
            if isempty(NormalizeBy)
                NormalizeBy = max(Data); % normalize by maximum
            end
            Data = bsxfun(@rdivide, Data, NormalizeBy);
    end
end


%% Plot activity
if isempty(hA)
    figure;
    hA = axes();
else
    axes(hA); hold on;
end

switch Type
    case 'image'
        imagesc(Data');
        set(hA, 'XTick', 0.5:frameRate:range(FrameIndex)+0.5, 'XTickLabel', cellstr(num2str((0:range(FrameIndex)/frameRate)')));
        xlabel(XLabel);
        
    case '2Dlines'
        
        % Shift data
        x = Data;
        if numROIs > 1
%             x = bsxfun(@minus, x, min(x));   % subtract off minimum
%             x = bsxfun(@rdivide, x, max(x)); % normalize by maximum
            x = bsxfun(@plus, x, 0.5:spacing:spacing*(numROIs-.5)); % space out lines
        end
                
        % Plot stimuli
        if ~isempty(ExpStimFrames)
            ExpStimFrames = bsxfun(@plus,ExpStimFrames,[-.5,.5]);
            hold on;
            YLim = [min(x(:)),max(x(:))]+[-.02,.02]*range(x(:));
            for tindex = 1:size(ExpStimFrames,1)
                patch(ExpStimFrames(tindex,[1,1,2,2])/frameRate,YLim([1,2,2,1]),[.9,.9,.9],'EdgeColor',[.9,.9,.9]);
            end
        end
        
        % Plot data
        h = plot(repmat((0:1/frameRate:range(FrameIndex)/frameRate)',1,numROIs),x,'LineWidth',LineWidth);
        if ~isempty(Colors)
            for rindex = 1:numROIs
                h(rindex).Color = Colors(rindex,:);
            end
        end
        xlabel(XLabel);
        ylabel('ROI');
        axis tight
        
    case '3Dlines'
        
        % Plot data
        [Y,X] = meshgrid(1:spacing:spacing*numROIs,0:1/frameRate:range(FrameIndex)/frameRate);
        h=plot3(X,Y,Data);
        hold on;
        if ~isempty(Colors)
            for rindex = 1:numROIs
                h(rindex).Color = Colors(rindex,:);
            end
        end
        xlabel(XLabel);
        ylabel('ROI');
        
        % Plot area underneath curves
        ZLim = get(gca,'ZLim');
        for rindex = 1:numROIs
            fill3(X([1:end,end,1],rindex),Y([1:end,end,1],rindex),[Data(:,rindex);ZLim(1);ZLim(1)],[1,1,1],'EdgeAlpha',0);
        end
        
        % Plot stimuli
        if exist('Stim','var')
            hold on;
            Stim = Stim-[0;Stim(1:end-1)];
            temp = [find(Stim==1),find(Stim==-1)];
            for sindex = 1:size(temp,1)
                fill3([temp(sindex,:),fliplr(temp(sindex,:))],repmat((numROIs+1)*spacing,1,2),[ZLim,fliplr(ZLim)],[.9,.9,.9]);
            end
        end
        
        axis tight
end



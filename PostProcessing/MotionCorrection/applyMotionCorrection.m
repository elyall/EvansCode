function [Images, MCdata, loadObj] = applyMotionCorrection(Images, MCdata, loadObj, varargin)


MCdataIndex = [];       % N x 2 matrix specifying the images file and depth for each MCdata object input
FileID = [];            % numFrames x 2 matrix specifying which MCdata object each frame corresponds to and what frame number within that object
FrameIndex = [1 inf];   % vector specifying which frames in Images to motion correct
DepthIndex = [1 inf];   % vector specifying which depths in Images to motion correct


directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'MCdataIndex'
                MCdataIndex = varargin{index+1};
                index = index + 2;
            case 'FileID'
                FileID = varargin{index+1};
                index = index + 2;
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
                index = index + 2;
            case 'DepthIndex'
                DepthIndex = varargin{index+1};
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

if ~exist('Images','var') || isempty(Images) % prompt for file selection
    [Images,directory] = uigetfile({'*.imgs;*.sbx;*.tif'},'Select images:',directory,'MultiSelect','on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(directory,Images);
end
if ischar(Images) || iscellstr(Images)
    [Images, loadObj] = load2P(Images,'Frames',FrameIndex); % load in images
    FrameIndex = [1 inf]; % update frame index since only frames requested were loaded
end
[H, W, numDepths, numChannels, numFrames] = size(Images); % determine image dimensions

if ~exist('MCdata','var') || isempty(MCdata)
    [MCdata,p] = uigetfile({'*.mat'},'Choose Experiment file',directory,'MultiSelect','on');
    if isnumeric(MCdata)
        return
    end
    MCdata = fullfile(p,MCdata);
end
if ischar(MCdata)
    MCdata = {MCdata};
end

if ~exist('loadObj','var') || isempty(loadObj)
    loadObj.FrameIndex = [ones(numFrames,1),(1:numFrames)',reshape(1:numFrames*numDepths,numDepths,numFrames)']; % assumes frames are from only 1 file input
end


%% Determine data to motion correct
if FrameIndex(end)==inf
    FrameIndex = [FrameIndex(1:end-1), FrameIndex(end-1)+1:numFrames];
end
numFrames = numel(FrameIndex);
if DepthIndex(end)==inf
    DepthIndex = [DepthIndex(1:end-1), DepthIndex(end-1)+1:numDepths];
end
numDepths = numel(DepthIndex);


%% Create MCdata object dictionary
if iscellstr(MCdata) % load in MCdata from files
    for findex = 1:numel(MCdata)
        temp = load(MCdata{findex},'MCdata','-mat');
        if ~isfield(temp,'MCdata')
            error('File %s does not contain an MCdata variable',MCdata{findex});
        end
        MCdata{findex} = temp.MCdata;
    end
end
if iscell(MCdata) % convert to struct
    MCdata = cat(2,MCdata{:});
end
if isempty(MCdataIndex) % build dictionary
    numMCFiles = numel(MCdata);
    numDataFiles = numel(unique(loadObj.FrameIndex(FrameIndex,1)));
    if numMCFiles ~= numDataFiles*numDepths
        error('Specify MCdata object dictionary');
    end
    MCdataIndex = [repelem(1:numDataFiles,numDepths)',repmat(DepthIndex',numDataFiles,1)];
end
    
    
%% Determine indices for each frame
if isempty(FileID)
    FileID = loadObj.FrameIndex(:,1:2);
end


%% Apply motion correction
switch MCdata(1).type
    
    case 'doLucasKanade'
        for iF = 1:numFrames
            for iC = 1:numChannels
                for iD = 1:numDepths
                    
                    % to specify depth of image passed in
                    if ~exist('depthIndex', 'var')
                        Di = iD;
                    else
                        Di = depthIndex;
                    end
                    
                    if iF == 1 % initialize variables
                        [Images(:,:,iD,iC,iF), B, xi, yi] = applyDoLucasKanade(...
                            Images(:,:,iD,iC,iF),...
                            MCdata(FileIndex(iF,1)).dpx(:,FileIndex(iF,2),Di),...
                            MCdata(FileIndex(iF,1)).dpy(:,FileIndex(iF,2),Di));
                    else % use previously created variables
                        Images(:,:,iD,iC,iF) = applyDoLucasKanade(...
                            Images(:,:,iD,iC,iF),...
                            MCdata(FileIndex(iF,1)).dpx(:,FileIndex(iF,2),Di),...
                            MCdata(FileIndex(iF,1)).dpy(:,FileIndex(iF,2),Di), B, xi, yi);
                    end
                end
            end
        end
        
    case 'Translation'
        
        % Determine limits of regions to fill with NaNs
        for index = 1:numel(MCdata)
            T = MCdata(index).T;
            T = repmat(T,1,2);
            T(T(:,1)>0,1) = 1; % positive offset: min is front edge and max is offset
            T(T(:,1)<0,[1,3]) = [H+T(T(:,1)<0,1)+1,repmat(H,nnz(T(:,1)<0),1)]; % negative offest: min is far edge minus offset and max is far edge
            T(T(:,2)>0,2) = 1; % positive offset: min is front edge and max is offset
            T(T(:,2)<0,[2,4]) = [W+T(T(:,2)<0,2)+1,repmat(W,nnz(T(:,2)<0),1)]; % negative offest: min is far edge minus offset and max is far edge
            MCdata(index).bounds = T;
        end
        
        for iF = 1:numFrames
            cF = FrameIndex(iF);
            for iD = 1:numDepths
                try 
                    cD = DepthIndex(iD);
                    mind = ismember(MCdataIndex,[FileID(iF,1),cD],'rows');
                    if any(MCdata(mind).T(FileID(cF,2),:))
                        Images(:,:,cD,:,cF) = circshift(Images(:,:,cD,:,cF), MCdata(mind).T(FileID(cF,2),1:2));
                        if MCdata(mind).T(FileID(cF,2),1)
                            Images(MCdata(mind).bounds(FileID(cF,2),1):MCdata(mind).bounds(FileID(cF,2),3),:,cD,:,iF) = nan; % fill top-bottom non-sampled region with nans
                        end
                        if MCdata(mind).T(FileID(cF,2),2)
                            Images(:,MCdata(mind).bounds(FileID(cF,2),2):MCdata(mind).bounds(FileID(cF,2),4),cD,:,iF) = nan; % fill left-right non-sampled region with nans
                        end
                        % Images(:,:,iD,:,iF) = imtranslate(Images(:,:,1,:,iF), MCdata(MCindex(iF,1)).T(MCindex(iF,2),[2,1]), 'FillValues',nan); % very slow
                    end
                catch ME
                    if FileID(cF,2) ~= size(MCdata(mind).T,1)+1 % bug where code doesn't motion correct last image -> leave last image as is
                        rethrow(ME)
                    end
                end
            end
        end
end
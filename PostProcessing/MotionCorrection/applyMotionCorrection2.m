function Images = applyMotionCorrection(Images, MCdata, loadObj, varargin)


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

if ~exist('Images','var') || isempty(Images) % Prompt for file selection
    [ImageFile,directory] = uigetfile({'*.imgs;*.sbx;*.tif'},'Select images:',directory,'MultiSelect','on');
    if isnumeric(ImageFile)
        return
    end
    ImageFile = fullfile(directory, ImageFile);
    [Images, loadObj] = load2P(ImageFile,'Frames',FrameIndex);
    FrameIndex = [1 inf];
elseif ischar(Images) || iscellstr(Images)
    [Images, loadObj] = load2P(Images,'Frames',FrameIndex);
    FrameIndex = [1 inf];
end
[H, W, numDepths, numChannels, numFrames] = size(Images);

if ~exist('MCdata','var') || isempty(MCdata)
    [ExperimentFile, p] = uigetfile({'*.mat'},'Choose Experiment file',directory,'MultiSelect','on');
    if isnumeric(ExperimentFile)
        return
    end
    ExperimentFile = fullfile(p,ExperimentFile);
    load(ExperimentFile, 'MCdata');
elseif ischar(MCdata)
    load(MCdata, 'MCdata','-mat');
end

if ~exist('loadObj','var') || isempty(loadObj)
    loadObj.FrameIndex = [ones(numFrames,1),(1:numFrames)',reshape(1:numFrames,numDepths,numFrames)']; % assumes frames are from only 1 file input
end


%% Determine data to motion correct
if FrameIndex(end)==inf
    FrameIndex = [FrameIndex(1:end-1), FrameIndex(end-1)+1:numFrames];
end
if DepthIndex(end)==inf
    DepthIndex = [DepthIndex(1:end-1), DepthIndex(end-1)+1:numDepths];
end


%% Create MCdata object dictionary
if isempty(MCdataIndex)
    numMCFiles = numel(MCdata);
    numDataFiles = numel(unique(loadObj.FileIndex(FrameIndex,1)));
    if numMCFiles ~= numDataFiles*numDepths
        error('Specify MCdata object dictionary');
    end
    MCdataIndex = [repelem(1:numDataFiles,numDepths)',repmat((1:numDepths)',numDataFiles,1)];
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
        for index = 1:numel(MCdata)
            T = MCdata(index).T;
            T = repmat(T,1,2);
            T(T(:,1)>0,1) = 1;
            T(T(:,2)>0,2) = 1;
            T(T(:,1)<0,[1,3]) = [H+T(T(:,1)<0,1)+1,repmat(H,nnz(T(:,1)<0),1)];
            T(T(:,2)<0,[2,4]) = [W+T(T(:,2)<0,2)+1,repmat(W,nnz(T(:,2)<0),1)];
            MCdata(index).T = T;
        end
        
        for iF = 1:numFrames
            try
                Images(:,:,1,:,iF) = circshift(Images(:,:,1,:,iF), MCdata(FileIndex(iF,1)).T(FileIndex(iF,2),1:2));
                if MCdata(FileIndex(iF,1)).T(FileIndex(iF,2),1)
                    Images(MCdata(FileIndex(iF,1)).T(FileIndex(iF,2),1):MCdata(FileIndex(iF,1)).T(FileIndex(iF,2),3),:,1,:,iF) = nan; % fill top-bottom non-sampled region with nans
                end
                if MCdata(FileIndex(iF,1)).T(FileIndex(iF,2),2)
                    Images(:,MCdata(FileIndex(iF,1)).T(FileIndex(iF,2),2):MCdata(FileIndex(iF,1)).T(FileIndex(iF,2),4),1,:,iF) = nan; % fill left-right non-sampled region with nans
                end
                % Images(:,:,1,iC,iF) = imtranslate(Images(:,:,1,iC,iF), MCdata(MCindex(iF,1)).T(MCindex(iF,2),[2,1]), 'FillValues',nan); % very slow
            catch ME
                if FileIndex(iF,2) ~= size(MCdata(FileIndex(iF,1)).T,1)+1 % old bug where code doesn't motion correct last image
                    rethrow(ME)
                end
            end
        end
end
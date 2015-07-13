function [Images, Config, InfoFile] = readSbx(SbxFile, InfoFile, varargin)
% Loads 'Frames' of single .sbx file ('SbxFile'). Requires
% corresponding information file ('InfoFile').

LoadType = 'Direct'; % 'MemMap' or 'Direct'
Frames = 1:20; % indices of frames to load in 'Direct' mode, or 'all'
Channels = 1;
SaveDirect = false;

% Defaults
invert = true;
xavg = 4; %scanbox version 1 only

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Type','type'}
                LoadType = varargin{index+1}; %'Direct' or 'MemMap'
                index = index + 2;
            case {'Frames','frames'} % indices of frames to load in 'Direct' mode
                Frames = varargin{index+1};
                index = index + 2;
            case {'Channels','channels'}
                Channels = varargin{index+1};
                index = index + 2;
            case {'Invert', 'invert'}
                invert = varargin{index+1};
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

if ~exist('SbxFile', 'var') || isempty(SbxFile)
    directory = CanalSettings('DataDirectory');
    [f,p] = uigetfile({'*.sbx'}, 'Choose scanbox file to load', directory);
    if isnumeric(f)
        Images = []; return
    end
    SbxFile = fullfile(p,f);
end

if ~exist('InfoFile', 'var') || isempty(InfoFile)
    InfoFile = sbxIdentifyFiles(SbxFile);
    InfoFile = InfoFile{1};
end

%% Load In Acquisition Information
load(InfoFile, 'info'); %load in 'info' variable

% Update header
if any(~isfield(info, {'Height', 'numChannels', 'numFrames', 'Width', 'scanbox_version'}))
    info = updateInfoFile(SbxFile, InfoFile);
end
% Parse acquisition information
Config = parseSbxHeader(InfoFile);
Config.loadType = LoadType;

%% Load In Images
switch LoadType
    
    case {'MemMap', 'memmap'}
        
        Images = MappedTensor(SbxFile, Config.size, 'Class', Config.Precision);
        Images = permute(Images, [3 2 5 1 4]);
        % Images = intmax('uint16') - Images; % computation doesn't carry
        % over
        Config.DimensionOrder = Config.DimensionOrder([3 2 5 1 4]);
        Config.size = size(Images);
        % Images = memmapfile(SbxFile,...
        %     'Format', {'uint16', [info.numChannels info.Width info.Height info.numFrames], 'Frames'}); %format data as best as possible into individual frames
        % 'Offset', nsamples*2*info.nchan,... %skip first frame as it's incomplete
        % Frames = permute(intmax('uint16')-x.Data.Frames,[3,2,5,1,4]);
        
    case {'Direct', 'direct'}
        
        % Determine frames to load
        if ischar(Frames) || (numel(Frames)==1 && Frames == inf)
            Frames = 1:info.numFrames;
        elseif Frames(end) == inf
            Frames = [Frames(1:end-2),Frames(end-1):info.numFrames];
        end
        numFramesToLoad = numel(Frames);
        Frames = Frames - 1; % offset b/c of "0" indexing
        
        % Determine channels to load
        if ischar(Channels) || (numel(Channels)==1 && Channels == inf)
            Channels = 1:info.numChannels;
        elseif Channels(end) == inf
            Channels = [Channels(1:end-2),Channels(end-1):info.numChannels];
        end
        numChannelsBeingRead = numel(Channels);
        
        % Preallocate output
        if info.scanbox_version == 1
            Images = zeros(numChannelsBeingRead, info.Width/xavg, info.Height, numFramesToLoad, 'uint16');
        elseif info.scanbox_version == 2
            Images = zeros(numChannelsBeingRead, info.Width, info.Height, numFramesToLoad, 'uint16');
        end

        
        % Determine number and location of seek operations due to non-contiguous frame requests (jumps)
        seekoperations = find(diff(Frames)~=1); %find any jumps within the frames to read out (jumps requiring seeking)
        if isempty(seekoperations) %no jumps
            numframesperread = numFramesToLoad; %all frames will be read in one read
            seekoperations = 1; %only one seek operation to first frame of FrameIndex
        else
            numframesperread = diff([0,seekoperations,numFramesToLoad]); %multiple reads required with various numbers of frames per read
            seekoperations = [1,seekoperations+1]; %indexes the first frame of each read within FrameIndex
        end
        
        % Open File
        info.fid = fopen(SbxFile);
        if(info.fid ~= -1)
            
            % Load Images
            fprintf('Loading\t%d\tframe(s) from\t%s...', numel(Frames), SbxFile);
            for index = 1:length(seekoperations)
                if(fseek(info.fid, info.Height * info.Width * 2 * Frames(seekoperations(index)) * info.numChannels, 'bof')==0) % "2" b/c assumes uint16 => 2 bytes per record
                    temp = fread(info.fid, info.Height * info.Width * info.numChannels * numframesperread(index), 'uint16=>uint16');
                    if info.scanbox_version == 1 % downsample/averaging
                        temp = mean(reshape(temp,[info.numChannels xavg info.Width/xavg info.Height numframesperread(index)]), 2);
                        info.Width = info.Width/xavg;
                        temp = reshape(temp,[info.numChannels info.Width info.Height numframesperread(index)]);
                    elseif info.scanbox_version == 2
                        temp = reshape(temp,[info.numChannels info.Width info.Height numframesperread(index)]);
                    end
                    Images(:,:,:,seekoperations(index):seekoperations(index)+numframesperread(index)-1) = temp(Channels,:,:,:); % save only requested channels
                else
                    warning('fseek error...');
                    Images = [];
                end
            end
            
            % Reorder dimensions
            Images = permute(Images, [3 2 5 1 4]); % flip colormap and reorder (original [1,3,2,4] => rotate images)
            
            % Flip colormap
            if invert
                Images = intmax('uint16') - Images;
            end
            
            % Correct for nonuniform spatial sampling
            if info.scanbox_version == 1
                S = sparseint;
                info.Width = size(S,2);
                good = zeros(info.Height, info.Width, 1, numChannelsBeingRead, info.numFrames);
                for ii = 1:numChannelsBeingRead
                    for iii = 1:1
                        for iiii = 1:info.numFrames
                            good(:,:,iii,ii,iiii) = Images(:,:,iii,ii,iiii) * S; % correct for non-uniform sampling
                        end
                    end
                end
                Images = good;
            end

        else
            warning('unable to open file: %s', SbxFile);
            Images = [];
        end
        fclose(info.fid);
        
        Config.DimensionOrder = Config.DimensionOrder([3 2 5 1 4]);
        if info.scanbox_version == 1
            Config.Width = info.Width;
        end
        Config.size = size(Images);
        
        fprintf('\tComplete\n');
        if SaveDirect
            fprintf('\tSaving frames to mat file');
            savefn = [SbxFile(1:end-4),'_EL.mat'];
            n=1;
            while exist(savefn,'file')
                savefn = [SbxFile(1:end-4),sprintf('%d_EL.mat',n)];
                n = n + 1;
            end
            save(savefn, 'Images', 'info', '-v7.3');
        end
end
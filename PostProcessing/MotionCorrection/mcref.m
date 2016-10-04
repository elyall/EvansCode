function refframe = mcref(files,refChannel,runChannel,auto,nFsame)

threshspeed=8;

%% Initialization
if ~exist('auto','var')
    auto=true;
end
if ~exist('same','var')
    nFsame=1;
elseif ~nFsame
    nFsame=1:length(files);
end
try %in case file(s) not saved with scim header
    [~,nFrames,ChannelIndex]=parseScim(files(nFsame));
catch
    nFrames=zeros(size(nFsame));
    for i=nFsame
        nFrames(i)=numel(imfinfo(files{i}));
    end
    ChannelIndex=[1,0,0,0]; % assume only first channel if scim header doesn't exit
end
if ~exist('refChannel','var') || isempty(refChannel)
    if ChannelIndex(2) % assumes secondary fluorophore is second channel
        refChannel=2;
    else
        refChannel=1;
    end
end
if ~exist('runChannel','var') || isempty(runChannel)
    if ChannelIndex(3) % assumes run speed is third channel
        runChannel=3;
    else
        runChannel=false;
    end
end

%% Code
if ispc && ~auto
    refframe=[]; % reference frame
    i = 1; % file index
    while ~isinteger(refframe)
        % update file index if beyond either edge
        if i < 1
            i = length(files);
        elseif i > length(files)
            i = 1;
        end
        % determine frames to load
        if strcmp(Frames,'all')
            imgindex = Channel:nChannels:numel(imfinfo(files{i}));
        else
            imgindex = Frames;
        end
        % load frames
        imgs = [];
        for j = imgindex
            imgs = cat(3,imgs,imread(files{i},'tif',j));
        end
        stackdisplay(imgs); % display frames
        refframe = input('What frame (#) should be the baseline? (F=next file, B=last file): ','s');
        close(gcf);
        % update file index or specify reference frame to load
        if strcmp(refframe,'F')
            i = i + 1;
        elseif strcmp(refframe,'B')
            i = i - 1;
        else
            refframe = str2double(refframe);
            refframe = uint16(refframe); % isinteger() doesn't like double
        end
    end
    refframe = imread(files{i},'tif',refframe*nChannels-(nChannels-Channel));
    refframe = double(refframe);
else % auto find frame to reference off of
    m=round(length(files)/2);
    if runChannel % find frame with running little or no running
        go=0;
        while ~go
            runspeed=brunspeed(files(m));
            f=find(runspeed<threshspeed,1);
            if ~isempty(f)
                go=1;
            else
                m=m+1;
                if m>length(files)
                    m=1;
                end
            end
        end
    else
        if numel(nFrames)>1
            nFrames=nFrames(m);
        end
        f=round(nFrames/2);
    end
    fprintf('\nloading reference frame from middle of stack: %s \t frame %d',files{m},f);
    refframe = load2P(files(m),refChannel,f);
end
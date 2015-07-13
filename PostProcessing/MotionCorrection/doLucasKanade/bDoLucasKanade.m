function OutFiles = bDoLucasKanade(ImgsFiles, ItemToSave)

% User settings
Nbasis = 16;
niter = 25;
damping = 1;
deltacorr = .0005;
Channel2AlignFrom = 1;

%% Check input arguments
narginchk(0,2);
if ~exist('ImgsFiles','var') || isempty(ImgsFiles) % Prompt for file selection
    directory = loadCanalSettings;
    [ImgsFiles,p] = uigetfile({'*.imgs'}, 'Choose ''imgs'' file to motion correct:', directory, 'MultiSelect', 'on');
    if isnumeric(ImgsFiles)
        return
    elseif iscellstr(ImgsFiles)
        for index = 1:length(ImgsFiles)
            ImgsFiles{index} = fullfile(p, ImgsFiles{index});
        end
    else
        ImgsFiles = {fullfile(p,ImgsFiles)};
    end
elseif ischar(ImgsFiles) % File is specified
    ImgsFiles = {ImgsFiles};
end
if ~exist('ItemToSave', 'var') || isempty(ItemToSave)
    ItemToSave = 'Frames'; % 'Map' or 'Frames' or 'Both'
end


%% Determine what to save
switch ItemToSave
    case 'Frames'
        saveimgs = true;
        savemap = false;
    case 'Map'
        saveimgs = false;
        savemap = true;
    case 'Both'
        saveimgs = true;
        savemap = true;
end
    

%% Determine number of files
nFiles = length(ImgsFiles);
if saveimgs
    OutImgsFiles = cell(nFiles,1);
end
if savemap
    AlignFiles = cell(nFiles,1);
end

%% Cycle through files aligning each frame of each file
for F = 1:nFiles
  
    % Determine # of frames to read in at a time
    Images = load2P(ImgsFiles{F}, 'Type', 'Direct', 'Frames', 2);
    mem = memory;
    sizeFrame = whos('Images');
    sizeFrame = sizeFrame.bytes;
    nFramesPerLoad = max(1, floor(0.1*mem.MaxPossibleArrayBytes/sizeFrame));
    
    % Load Header
    config = load2PConfig(ImgsFiles{F});
    if Channel2AlignFrom > config.Channels;
        Channel2AlignFrom = 1;
    end
    
    % initialize imgs file
    if saveimgs
        OutImgsFiles{F} = [strtok(ImgsFiles{F}, '.'), '_dLK.imgs'];
        header = createImgsHeader(config,...
            'type',     'imgs',... 
            'mcType', 'doLucasKanade',...
            'mcReferenceChannel', Channel2AlignFrom,...
            'mcNbasis', Nbasis,...
            'mcNiter', niter,...
            'mcDamping', damping,...
            'mcDeltaCorr', deltacorr...
            );
        fid = writeImgsHeader(header, OutImgsFiles{F});
        Channels = 1:config.Channels;
    end
    
    % initialize map file
    if savemap
        AlignFiles{F} = [strtok(ImgsFiles{F}, '.'), '.align'];
        MC.type = 'doLucasKanade';
        MC.Nbasis = Nbasis;
        MC.niter = niter;
        MC.damping = damping;
        MC.deltacorr = deltacorr;
        MC.Channel2AlignFrom = Channel2AlignFrom;
    end
    MC.dpx = zeros(Nbasis + 1, config.Frames);
    MC.dpy = zeros(Nbasis + 1, config.Frames);
    
    % Compute template
    Frames = load2P(ImgsFiles{F}, 'Type', 'Direct', 'Frames', 1:500);
    template = double(mean(Frames(:,:,:,Channel2AlignFrom,:),5));
    
    % For each frame determine the necessary transformation
    wb=waitbar(0, sprintf('Aligning %d of %d: %s', F, nFiles, ImgsFiles{F}));
    
    % Load Batch
    for f = 1:nFramesPerLoad:config.Frames
        lastframe = min(f+nFramesPerLoad-1, config.Frames);
        Frames = load2P(ImgsFiles{F}, 'Type', 'Direct', 'Frames', f:lastframe);
        
        % Compute Transformation
        for currentFrame = 1:lastframe-f+1
            FrameIndex = f+currentFrame-1;
            [AlignedFrame, MC.dpx(:,FrameIndex), MC.dpy(:,FrameIndex)] = doLucasKanade(template, double(Frames(:,:,1,Channel2AlignFrom,currentFrame)));
            %         if f > 15 % repeat alignment to last 15 frames
            %            [MC.frames(:,:,1,Channel2AlignFrom,f),MC.dpx(:,f),MC.dpy(:,f)] = doLucasKanade(double(mean(Frames(:,:,1,Channel2AlignFrom,f-15:f-1),5)), double(Frames(:,:,1,Channel2AlignFrom,f)), MC.dpx(:,f), MC.dpy(:,f));
            %         end
            
            % Save motion corrected images
            if saveimgs
                % Save frame from each channel in order
                for c = Channels
                    if c == Channel2AlignFrom
                        % Save already computed image
                        img = reshape(AlignedFrame,config.Height*config.Width,1); % reshape back to a vector for saving
                        fwrite(fid,img,config.Precision);
                    else
                        % Align frame from current channel
                        if f == 1
                            [currentframe, B, xi, yi] = applyDoLucasKanade(double(Frames(:,:,1,c,currentFrame)), MC.dpx(:,FrameIndex), MC.dpy(:,FrameIndex), Nbasis);
                        else
                            currentframe = applyDoLucasKanade(double(Frames(:,:,1,c,f)), MC.dpx(:,FrameIndex), MC.dpy(:,FrameIndex), Nbasis, B, xi, yi);
                        end
                        img = reshape(currentframe,config.Height*config.Width,1); % reshape back to a vector for saving
                        fwrite(fid,img,config.Precision);
                    end
                end
            end
            
            waitbar(FrameIndex/config.Frames, wb);
        end
    end
    close(wb);
    
    % Save imgs: close file that was created
    if saveimgs
        fclose(fid);
    end
    
    % Save map
    if savemap
        save(AlignFiles{F}, 'MC', '-mat', '-v7.3');
    end
    
    % Format output
    if saveimgs && savemap
        OutFiles = [OutImgsFiles, AlignFiles];
    elseif saveimgs
        OutFiles = OutImgsFiles;
    elseif savemap
        OutFiles = AlignFiles;
    end
    
    clear template Frames % clear MemMap variables
    clear Images
end
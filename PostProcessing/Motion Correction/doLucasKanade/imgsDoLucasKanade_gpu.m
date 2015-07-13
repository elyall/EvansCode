function [OutImgsFile, AlignFile] = imgsDoLucasKanade_gpu(ImgsFile, ItemToSave, Channel2AlignFrom, Nbasis, niter, damping, deltacorr)


%% Check input arguments
narginchk(0,7);
if ~exist('ImgsFile','var') || isempty(ImgsFile) % Prompt for file selection
    directory = CanalSettings('DataDirectory');
    [ImgsFile,p] = uigetfile({'*.imgs'}, 'Choose ''imgs'' file to motion correct:', directory, 'MultiSelect', 'on');
    if isnumeric(ImgsFile)
        return
    end
    ImgsFile = fullfile(p,ImgsFile);
end

if ~exist('ItemToSave', 'var') || isempty(ItemToSave)
    ItemToSave = 'Frames'; % 'Map' or 'Frames' or 'Both'
end

if ~exist('Channel2AlignFrom', 'var') || isempty(Channel2AlignFrom)
    Channel2AlignFrom = 1;
end

if ~exist('Nbasis', 'var') || isempty(Nbasis)
    Nbasis = 16;
end

if ~exist('niter', 'var') || isempty(niter)
    niter = 25;
end

if ~exist('damping', 'var') || isempty(damping)
    damping = 1;
end

if ~exist('deltacorr)', 'var') || isempty(deltacorr)
    deltacorr = .0005;
end


%% Determine what to save
OutImgsFile = [];
AlignFile = [];
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


%% Cycle through files aligning each frame of each file
try
    
    % load images and header
    [Images, config] = loadImgs(ImgsFile);
    Frames = Images.Data.Frames;
    if Channel2AlignFrom > config.Channels;
        Channel2AlignFrom = 1;
    end
    
    % initialize imgs file
    if saveimgs
        OutImgsFile = [strtok(ImgsFile, '.'), '_dLK.imgs'];
        header = createImgsHeader(config,...
            'mcType', 'doLucasKanade',...
            'mcReferenceChannel', Channel2AlignFrom,...
            'mcNbasis', Nbasis,...
            'mcNiter', niter,...
            'mcDamping', damping,...
            'mcDeltaCorr', deltacorr...
            );
        fid = writeImgsHeader(header, OutImgsFile);
        Channels = 1:config.Channels;
    end
    
    % initialize map file
    if savemap
        AlignFile = [strtok(ImgsFile, '.'), '.align'];
        MC.type = 'doLucasKanade';
        MC.Nbasis = Nbasis;
        MC.niter = niter;
        MC.damping = damping;
        MC.deltacorr = deltacorr;
        MC.Channel2AlignFrom = Channel2AlignFrom;
    end
    MC.dpx = zeros(Nbasis + 1, config.Frames);
    MC.dpy = zeros(Nbasis + 1, config.Frames);
    
    % For each frame determine the necessary transformation
    template = double(mean(Frames(:,:,1,Channel2AlignFrom,1:min(500,end)),5));
    wb=waitbar(0, sprintf('Aligning: %s', ImgsFile));
    for f = 1:config.Frames
        [AlignedFrame, MC.dpx(:,f), MC.dpy(:,f)] = doLucasKanade_gpu(template, double(Frames(:,:,1,Channel2AlignFrom,f)), [], [], Nbasis, niter, damping, deltacorr);
        %         if f > 15 % repeat alignment to last 15 frames
        %            [MC.frames(:,:,1,Channel2AlignFrom,f),MC.dpx(:,f),MC.dpy(:,f)] = doLucasKanade_gpu(double(mean(Frames(:,:,1,Channel2AlignFrom,f-15:f-1),5)), double(Frames(:,:,1,Channel2AlignFrom,f)), MC.dpx(:,f), MC.dpy(:,f));
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
                        [currentframe, B, xi, yi] = applyDoLucasKanade_gpu(double(Frames(:,:,1,c,f)), MC.dpx(:,f), MC.dpy(:,f), Nbasis);
                    else
                        currentframe = applyDoLucasKanade_gpu(double(Frames(:,:,1,c,f)), MC.dpx(:,f), MC.dpy(:,f), Nbasis, B, xi, yi);
                    end
                    img = reshape(currentframe,config.Height*config.Width,1); % reshape back to a vector for saving
                    fwrite(fid,img,config.Precision);
                end
            end
        end
        
        waitbar(f/config.Frames, wb);
    end
    close(wb);
    
    % Save imgs: close file that was created
    if saveimgs
        fclose(fid);
    end
    
    % Save map
    if savemap
        save(AlignFile, 'MC', '-mat', '-v7.3');
    end
    
    % Format output
    if saveimgs && savemap
        OutFiles = [OutImgsFile, AlignFile];
    elseif saveimgs
        OutFiles = OutImgsFile;
    elseif savemap
        OutFiles = AlignFile;
    end
    
catch exception
    try close(wb); end
    warning('Failed to align file %d: %s\n%s', F, ImgsFile, exception.message);
end

clear template Frames % clear MemMap variables
clear Images
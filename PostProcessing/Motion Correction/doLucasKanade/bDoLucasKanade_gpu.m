function OutFiles = bDoLucasKanade_gpu(ImgsFiles, ItemToSave)

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
    
    try
        
        % load images and header
        [Images, config] = loadImgs(ImgsFiles{F});
        Frames = Images.Data.Frames;
        if Channel2AlignFrom > config.Channels;
            Channel2AlignFrom = 1;
        end
        
        % initialize imgs file
        if saveimgs
            OutImgsFiles{F} = [strtok(ImgsFiles{F}, '.'), '_dLK.imgs'];
            header = createImgsHeader(config,...
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
        
        % For each frame determine the necessary transformation
        template = double(mean(Frames(:,:,1,Channel2AlignFrom,1:min(500,end)),5));
        wb=waitbar(0, sprintf('Aligning %d of %d: %s', F, nFiles, ImgsFiles{F}));
        for f = 1:config.Frames
            [AlignedFrame, MC.dpx(:,f), MC.dpy(:,f)] = doLucasKanade_gpu(template, double(Frames(:,:,1,Channel2AlignFrom,f)));
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
        
    catch exception
        try close(wb); end
        warning('Failed to align file %d: %s\n%s', F, ImgsFiles{F}, exception.message);
    end
    
    clear template Frames % clear MemMap variables
    clear Images
end
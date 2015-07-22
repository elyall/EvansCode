function saveFile = Vid(saveFile,ImageFiles,frameRate,stimindex,speed,varargin)
% FN - filename to save video to
% FILES - cell array of strings of filenames, or stack of images with
% dimensions [X,Y,t]
% FS - framerate of the video (can be left blank)
% STIMINDEX - logical vector of length t specifying whether a stimbox
% should be located in that frame. Or can be location of stimulation file
% to calculate index off of, or can be 'true' to prompt for stim file
% SPEED - vector of length t specifying the running speed of the animal in
% each frame. Or can be location of stimulation file to load speed off of, 
% or can be 'true' to prompt for stim file. If stim file doesn't have
% speed, then speed will be calculated off original files
% ...,'Channel',S,... - S is a scalar specifying which Channel to load the
% frames from in the case that FILES is a cell array of strings


directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    switch varargin{index}
        case 'Channel'
            Channel = varargin{index+1};
            index = index + 2;
        case 'SpeedUp'
            SpeedUp = varargin{index+1};
            index = index + 2;
        otherwise
            warning('Argument ''%s'' not recognized',varargin{index});
            index = index + 1;
    end
end

if ~exist('ImageFiles','var') || isempty(ImageFiles)
    ImageFiles=listfns();
    if isempty(ImageFiles)
        return
    end
end
if ~exist('saveFile','var') || isempty(saveFile)
    if iscellstr(ImageFiles)
        [filepath,saveFile,~]=fileparts(ImageFiles{1});
    else
        filepath=cd;
        saveFile='';
    end
    [saveFile,p]=uiputfile({'*.avi'},'Save as',fullfile(filepath,strcat(saveFile,'.avi')));
    if isempty(saveFile)
        return
    end
    saveFile=fullfile(p, saveFile);
end

if ~exist('Channel','var')
    Channel=1;
end

if ~exist('SpeedUp','var')
    SpeedUp=3;
end

if ~exist('frameRate','var') || isempty(frameRate)
    frameRate = 15.45;
end

% if ~exist('stimindex','var') || isempty(stimindex) || (isscalar(stimindex) && stimindex==false)
    stimtoggle=0;
% else
%     stimtoggle=1;
%     if (isscalar(stimindex) && stimindex==true) || ischar(stimindex) % if filename given rather than indices of stim
%         if ischar(stimindex) && ~isdir(stimindex)
%             StimFile = stimindex;
%         else
%             if ischar(stimindex) && isdir(stimindex)
%                 prompt = stimindex; % change prompt directory
%             else
%                 prompt = [];
%             end
%             if iscellstr(ImageFiles)
%                 sites = uniquesites(ImageFiles);
%             else
%                 sites = [];
%             end
%             StimFile = stimid(sites,prompt); % prompt to select stimfile
%         end
%         stimindex = load(StimFile{1},'Experiment');
%         stimindex = stimindex.Experiment.StimulusIndex;
%         
%         % legacy Nov 2013
% %         if ~ischar(stimindex)
% %             % prompt for file name
% %             if ~exist('filepath','var')
% %                 filepath=cd;
% %             end
% %             [f,p]=uigetfile({'*.mat'},'Choose stim file',filepath); %.xls=legacy Nov2013
% %             stimfn=[p,f];
% %         end
% %         % create stimindex from file
% %         if iscellstr(files)
% %             currentfiles=files;
% %         else
% %             currentfiles=listfns();
% %         end
% %         [stimindex,~,~] = StimIndex(currentfiles,stimfn);
% %         stimindex=logical(stimindex);
%     end
% end
% if ~exist('speed','var') || isempty(speed) || (isscalar(speed) && speed==false)
    speedtoggle=0;
% else
%     speedtoggle=1;
%     if (isscalar(speed) && speed==true) || ischar(speed) % if filename given rather than indices of stim
%         if ~exist('StimFile','var')
%             if ischar(speed) && ~isdir(speed)
%                 StimFile = speed;
%             else
%                 if ischar(speed) && isdir(speed)
%                     prompt = speed; % change prompt directory
%                 else
%                     prompt = [];
%                 end
%                 if iscellstr(ImageFiles)
%                     sites = uniquesites(ImageFiles);
%                 else
%                     sites = [];
%                 end
%                 StimFile = stimid(sites,prompt); % prompt to select stimfile
%             end
%         end
%         speed = load(StimFile{1},'Experiment');
%         speed = speed.Experiment.RunningSpeed;
%         
%         % legacy Nov2013
% %         % get stim filename
% %         if ~ischar(speed)
% %             if ~exist('stimfn','var')
% %                 if ~exist('filepath','var')
% %                     filepath=cd;
% %                 end
% %                 [f,p]=uigetfile({'*.xls'},'Choose stim file',filepath);
% %                 stimfn=[p,f];
% %             end
% %         else
% %             stimfn=speed;
% %         end
% %         % create speed from file
% %         try
% %             speed=xlsread(stimfn,'runspeed');
% %         catch % in case run speed wasn't saved to stim file
% %             if iscellstr(files)
% %                 currentfiles=files;
% %             else
% %                 currentfiles=listfns();
% %             end
% %             speed=brunspeed(currentfiles);
% %         end
%     end
% end


%% Code
if iscellstr(ImageFiles)
   ImageFiles = load2P(ImageFiles,'Type','Direct','Channel',Channel); 
end
fprintf('\nframerate = %f, saving channel %d',frameRate,Channel);

video=VideoWriter(saveFile,'Motion JPEG AVI');
set(video,'FrameRate',frameRate);
open(video)

highval=prctile(ImageFiles(:),99.9);
[y,x,f]=size(ImageFiles);
maxdim=max(x,y);
for i=1:f
    img=imresize(ImageFiles(:,:,i),[maxdim,maxdim]); %make image square
    img=uint8(img/(highval/(double(intmax('uint8'))+1))) - 1; %convert image to uint8 (or double would work too)
%     if stimtoggle
%         if stimindex(i)
%             img(round(.9*maxdim):maxdim,round(.9*maxdim):maxdim)=255;
%         end
%     end
%     if speedtoggle
%         s=round(speed(i));
%         if s>0
%             img(maxdim-s:maxdim,1:max(round(0.01*maxdim),1))=255;
%         end
%     end
    writeVideo(video,img);
end
video.close;
fprintf('\nvideo saved: %s',saveFile);

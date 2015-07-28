function plotStimSep(d,Fs,ROIid,baseline,stimulus,framelimit,minrunspeed)
% need to include stimulus names

%% Intitialization
if ~exist('d','var')
    d = ROIorganize();
elseif iscellstr(d)
    d = ROIorganize(d);
end
if ~exist('Fs','var') || isempty(Fs)
    Fs=1;
    xlab='Frames';
else
    xlab='Time (s)';
end
if ~exist('ROIid','var') || isempty(ROIid)
    ROIid = 'all';
end
if ~exist('baseline','var') || isempty(baseline)
    baseline=5; % baseline frames
end
if ~exist('stimulus','var') || isempty(stimulus)
    stimulus=0; % stim period (in seconds if Fs exists, otherwise in frames)
end
if ~exist('framelimit', 'var')
    framelimit=[2,40];
end
if ~exist('minrunspeed','var')
    minrunspeed = 50;
end


%% Code
numfiles=length(d);
for f=1:numfiles
    for r=1:d(f).nrois
        if ischar(ROIid) || ismember(r,ROIid) %only display for ROIs selected
            data={};
            yMax=zeros(d(f).nstim,1);
            yMin=zeros(d(f).nstim,1);
            baselineframes = zeros(d(f).nstim,2);
            stimulusframes = zeros(d(f).nstim,2);
            avg={};
            for s=1:d(f).nstim
                
                if length(baseline) == 1
                    if ~isempty(d(f).stims(s).firstframe)
                        baselineframes(s,1) = d(f).stims(s).firstframe - baseline;
                        baselineframes(s,2) = d(f).stims(s).firstframe - 1;
                    else %control trial
                        baselineframes(s,1) = 2;
                        baselineframes(s,2) = 6;
                    end
                elseif length(baseline) == 2
                    baselineframes(s,1) = baseline(1);
                    baselineframes(s,2) = baseline(2);
                end
                if length(stimulus) == 1
                    if ~isempty(d(f).stims(s).firstframe)
                        stimulusframes(s,1) = d(f).stims(s).firstframe;
                        stimulusframes(s,2) = d(f).stims(s).lastframe;
                    else %control trial
                        stimulusframes(s,1) = 7;
                        stimulusframes(s,2) = 12;
                    end
                elseif length(stimulus) == 2
                    stimulusframes(s,1) = stimulus(1);
                    stimulusframes(s,2) = stimulus(2);
                end
                
                data{s} = d(f).rois(r).data(d(f).rois(r).stimindex==d(f).stimids(s),:); % pull out raw data
                data{s}(isnan(data{s})) = 0;
                speed = d(f).rois(r).speed(d(f).rois(r).stimindex==d(f).stimids(s),:);
                speed(isnan(speed)) = 0;
                avgspeed = mean(speed(:,12:24),2); %detemine average running speed during stimulus period
                data{s}(avgspeed<minrunspeed,:) = []; %remove trials where the mouse wasn't running
                baselineF = mean(data{s}(:,baselineframes(s,1):baselineframes(s,2)),2); % calculate baseline for each trial remaining (average fluorescence during pre-stimulus period)
                data{s} = bsxfun(@rdivide,bsxfun(@minus,data{s},baselineF),baselineF); % calculate DFoF
                avg{s}=mean(data{s},1);
                yMax(s)=max(data{s}(:));
                yMin(s)=min(data{s}(:));
            end
            yMax = max(yMax);
            yMin = min(yMin);
            figure('NumberTitle','off','Name',sprintf('ROI: %d, File: %d',r,f));
            for s=1:d(f).nstim
                subplot(1,d(f).nstim,s);
                plot(framelimit(1)/Fs:1/Fs:framelimit(2)/Fs,data{s}(:,framelimit(1):framelimit(2))','k','LineWidth',1); % plot raw data points      
                hold on
                if strcmp(d(f).rois(r).label, 'PV+')
                    plot(framelimit(1)/Fs:1/Fs:framelimit(2)/Fs,avg{s}(:,framelimit(1):framelimit(2)),'r','LineWidth',4);
                else
                    plot(framelimit(1)/Fs:1/Fs:framelimit(2)/Fs,avg{s}(:,framelimit(1):framelimit(2)),'b','LineWidth',4);
                end 
                % title(sprintf('%s',num2str(d(f).stimids(s))));
                if d(f).stimids(s) == 0
                    title('control');
                    ylabel('dF/F');
                    xlabel(sprintf('%s',xlab));
                else
                    title(sprintf('%d',s-1));
                    set(gca, 'YTickLabel', [])
                    xshift = -framelimit(1)+1;
                    h=area([(stimulusframes(s,1)-.5+xshift)/Fs,(stimulusframes(s,2)+.5+xshift)/Fs],[yMax,yMax],yMin);
                    child=get(h,'Children');
                    set(child,'FaceAlpha',0.2,'EdgeColor','none');
                    axis off
                end
                ylim([yMin,yMax]);
                xlim([0,d(f).minframespertrial/Fs]);                              
            end
        end
    end
end

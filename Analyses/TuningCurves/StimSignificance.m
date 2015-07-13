function ROIdata = StimSignificance(ROIdata,ROIid,baseline,stimulus,minrunspeed,neuropilweight,maxpvalue)

FitTuningCurves = true;

saveToROIFile = false;
saveToPDF = false;

ControlIndex = false;

%% Check input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    directory = CanalSettings('DataDirectory');
    [ROIdata, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('ROIid','var') || isempty(ROIid)
    ROIid = 'all';
end

if ~exist('baseline','var') || isempty(baseline)
    baseline=10; % baseline frames
end
if ~exist('stimulus','var') || isempty(stimulus)
    stimulus=[];
end
if ~exist('minrunspeed','var') || isempty(minrunspeed)
    minrunspeed=100;
end
if ~exist('neuropilweight', 'var') || isempty(neuropilweight)
    neuropilweight = 0.65;
end
if ~exist('maxpvalue','var') || isempty(maxpvalue)
    maxpvalue=0.05;
end
if ~exist('outlierweight', 'var') || isempty(outlierweight)
    outlierweight = 3;
end

%% Load data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
end
numROIs = numel(ROIdata.rois);
StimIDs = unique(ROIdata.DataInfo.StimID);
numStimuli = numel(StimIDs);

if neuropilweight
    if ~isfield(ROIdata.rois, 'neuropil')
        warning('Neuropil doesn''t exist for current file. Continuing anyway...');
        neuropilweight = 0;
    end
end

if minrunspeed>0
    if ~isfield(ROIdata, 'RunningSpeed')
        minrunspeed = 0;
    end
end

%% Determine Stimulus Frames
% StimFrames = zeros(numStimuli, 2);
% BaselineFrames = zeros(numStimuli, 2);
% for sindex = 1:numStimuli
%     temp = AnalysisInfo(AnalysisInfo.StimID==StimIDs(sindex), :);
%     StimFrames(sindex,:) = temp(1,:).TrialStimFrames;
%     BaselineFrames(sindex,:) = temp(1,:).TrialBaselineFrames;
% end

%% Calculate average response for each stimulus
for rindex = 1:numROIs
    ROIdata.rois(rindex).curve = zeros(1, numStimuli);
    ROIdata.rois(rindex).StdError = zeros(1, numStimuli);
    ROIdata.rois(rindex).nTrials = zeros(1, numStimuli);
    ROIdata.rois(rindex).Raw = cell(numStimuli, 1);
    if ControlIndex
        ROIdata.rois(rindex).SigPos = zeros(1, numStimuli);
        ROIdata.rois(rindex).PValue = zeros(1, numStimuli);
    end
    
    controlindex = find(StimIDs==ControlIndex); %locate control 'stimulus' in structure
    stimuliindex = find(StimIDs~=ControlIndex); %locate the rest of the stimuli
    indexorder = cat(1,controlindex,stimuliindex); %run through the stimuli analyzing the control trials first
    
    for s=indexorder'
        
        if length(baseline) == 2
            baselinefirstframe = baseline(1);
            baselinelastframe = baseline(2);
        else
            baselinefirstframe = BaselineFrames(s,1);
            baselinelastframe = BaselineFrames(s,2);
        end
        if length(stimulus) == 2
            stimfirstframe = stimulus(1);
            stimlastframe = stimulus(2);
        else
            stimfirstframe = StimFrames(s,1);
            stimlastframe = StimFrames(s,2);
        end
        
        CaTraces = ROIdata.rois(rindex).data(ROIdata.DataInfo.StimID==StimIDs(s),:); % pull out all trials for current stimulus
        CaTraces(isnan(CaTraces)) = 0;
        if neuropilweight % adjust for neuropil
            neuropil = ROIdata.rois(rindex).neuropil(ROIdata.DataInfo.StimID==StimIDs(s),:);
            CaTraces = CaTraces-neuropilweight*neuropil;
        end
        if minrunspeed > 0 % remove trials where mouse wasn't running
            speed = ROIdata.RunningSpeed(ROIdata.DataInfo.StimID==StimIDs(s),:);
            speed(isnan(speed)) = 0;
            avgspeed = mean(speed(:,stimfirstframe:stimlastframe),2);
            CaTraces = CaTraces(avgspeed>=minrunspeed,:); %only analyze trials where the mouse is running
        end
        
        BaselineData = CaTraces(:,baselinefirstframe:baselinelastframe);
        StimulusData = CaTraces(:,stimfirstframe:stimlastframe);
        BaselineAvg = mean(BaselineData,2);
        StimulusDFoF = bsxfun(@rdivide,bsxfun(@minus,StimulusData,BaselineAvg),BaselineAvg);
        StimulusDFoF = mean(StimulusDFoF,2);
        
        % Remove Outliers
        n = inf;
        while n ~= length(StimulusDFoF)
            n = length(StimulusDFoF);
            StimulusDFoF(abs(StimulusDFoF-mean(StimulusDFoF))>outlierweight*std(StimulusDFoF)) = [];
        end
        
        if s == ControlIndex %control trials
            ControlDFoF = StimulusDFoF;
        end
        
        % Save tuning curves
        ROIdata.rois(rindex).curve(s) = mean(StimulusDFoF); % record average dF/F for stimulus
        ROIdata.rois(rindex).StdError(s) = std(StimulusDFoF)/sqrt(length(StimulusDFoF)); % record standard error for stimulus
        ROIdata.rois(rindex).Raw{s} = StimulusDFoF;
        ROIdata.rois(rindex).nTrials(s) = length(StimulusDFoF);
        
        % Perform t-test
        if ControlIndex
            [ROIdata.rois(rindex).SigPos(s), ROIdata.rois(rindex).PValue(s)] = ttest2(...
                ControlDFoF,....
                StimulusDFoF,...
                'Alpha',maxpvalue);
        end
        
    end %cycle stimuli
    
    % Fit tuning curves
    if numStimuli > 1 && FitTuningCurves
        ROIdata.rois(rindex).min = min(ROIdata.rois(rindex).curve(2:end));
        ROIdata.rois(rindex).max = max(ROIdata.rois(rindex).curve(2:end));
        if ControlIndex
            excitatorycurve = ROIdata.rois(rindex).curve(2:end)-ROIdata.rois(rindex).min;
            inhibitorycurve = ROIdata.rois(rindex).curve(2:end)-ROIdata.rois(rindex).max;
        else
            excitatorycurve = ROIdata.rois(rindex).curve-ROIdata.rois(rindex).min;
            inhibitorycurve = ROIdata.rois(rindex).curve-ROIdata.rois(rindex).max;
        end
        [ExicFit, ExicGoFit] = FitFull(excitatorycurve);
        [InhFit, InhGoFit] = FitFull(inhibitorycurve);
        if ExicGoFit.rsquare >= InhGoFit.rsquare
            ROIdata.rois(rindex).FitDirection = 'Excitatory';
            ROIdata.rois(rindex).Fit = ExicFit;
            ROIdata.rois(rindex).GoFit = ExicGoFit;
            ROIdata.rois(rindex).offset = ROIdata.rois(rindex).min;
        else
            ROIdata.rois(rindex).FitDirection = 'Inhibitory';
            ROIdata.rois(rindex).Fit = InhFit;
            ROIdata.rois(rindex).GoFit = InhGoFit;
            ROIdata.rois(rindex).offset = ROIdata.rois(rindex).max;
        end
        ROIdata.rois(rindex).Coeff = coeffvalues(ROIdata.rois(rindex).Fit);
        ROIdata.rois(rindex).ConfIntervals = confint(ROIdata.rois(rindex).Fit);
        ROIdata.rois(rindex).rsquare = ROIdata.rois(rindex).GoFit.rsquare;
    end
    

    
    % Plot figures
    if (ischar(ROIid) && (strcmp(ROIid, 'all') || strcmp(ROIid, ROIdata.rois(rindex).label))) || (isnumeric(ROIid) && ismember(rindex,ROIid)) %only display for ROIs selected
        
        % Fix Label
        if isempty(ROIdata.rois(rindex).label)
            ROIdata.rois(rindex).label = {'none'};
        end
        
        hf = figure(); % create figure 'NumberTitle','off','Name',sprintf('ROI: %ROIdata, Label: %s',r,ROIdata.rois(r).label{1})
        hold on
        
        % Plot each trial's average dF/F for each stimulus
        plot(ones(ROIdata.rois(rindex).nTrials(1),1), ROIdata.rois(rindex).Raw{1}, 'k.') %plot raw data points for control position
        for s = 2:numStimuli
            plot((s)*ones(ROIdata.rois(rindex).nTrials(s),1), ROIdata.rois(rindex).Raw{s}, 'k.') %plot raw data points for all stimuli
        end
        
        % Plot tuning curves
        if (~isempty(ROIdata.rois(rindex).label) && ~strcmp(ROIdata.rois(rindex).label{1},'none')) && ~strcmp(ROIdata.rois(rindex).label, 'none') % plot in red
            if ControlIndex
                errorbar(1,ROIdata.rois(rindex).curve(1),ROIdata.rois(rindex).StdError(1),'r','LineWidth',2,'Marker','.','MarkerSize',20);  %plot control position
                errorbar(2:numStimuli,ROIdata.rois(rindex).curve(2:end),ROIdata.rois(rindex).StdError(2:end),'r-','LineWidth',2,'Marker','.','MarkerSize',20); %plot curve
            else
                errorbar(1:numStimuli,ROIdata.rois(rindex).curve,ROIdata.rois(rindex).StdError,'r-','LineWidth',2,'Marker','.','MarkerSize',20); %plot curve=
            end
        else % plot in blue
            if ControlIndex
                errorbar(1,ROIdata.rois(rindex).curve(1),ROIdata.rois(rindex).StdError(1),'bo','LineWidth',2,'Marker','.','MarkerSize',20);  %plot control position
                errorbar(2:numStimuli,ROIdata.rois(rindex).curve(2:end),ROIdata.rois(rindex).StdError(2:end),'b-','LineWidth',2,'Marker','.','MarkerSize',20); %plot curve
            else
                errorbar(1:numStimuli,ROIdata.rois(rindex).curve,ROIdata.rois(rindex).StdError,'b-','LineWidth',2,'Marker','.','MarkerSize',20); %plot curve=
            end
        end
        
        % Plot fit
        if numStimuli > 1
            if ControlIndex
                yfit = feval(ROIdata.rois(rindex).Fit,1:0.001:numStimuli-1);
                plot(2:0.001:numStimuli, yfit + ROIdata.rois(rindex).offset, 'g-','LineWidth',2);
            else
                yfit = feval(ROIdata.rois(rindex).Fit,1:0.001:numStimuli);
                plot(1:0.001:numStimuli, yfit + ROIdata.rois(rindex).offset, 'g-','LineWidth',2);
            end
            Ydim = get(gca,'YLim');
            Xdim = get(gca,'XLim');
            text(Xdim(2),Ydim(2),sprintf('r^2=%.2f\nFWHM=%.2f', ROIdata.rois(rindex).rsquare, ROIdata.rois(rindex).Coeff(3)),'Color','g','FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top');
        end
        
        % Set axes labels
%         set(gca, 'XTick', 1:numStimuli, 'XTickLabel', [{'control'}; cellstr(num2str((30:15:150)'))]);
        set(gca, 'XTick', 1:numStimuli, 'XTickLabel', cellstr(num2str(StimIDs)));
        xlabel('Angle');
%         set(gca,'XTick',1:numStimuli,'XTickLabel',[{'control'}; cellstr(num2str((1:7)'))]);
%         xlabel('Position');
        ylabel('Average Stimulus-Evoked dF/F');
        xlim([0,numStimuli+1]);
        
        % Plot stars for stimuli that evoked a significant response
%         Ydim = get(gca,'YLim');
%         sigindices = find(ROIdata.rois(r).SigPos); %locate any responses that are significant
%         if ~isempty(sigindices)
%             text(sigindices,(Ydim(2)-(Ydim(2)-Ydim(1))/10)*ones(length(sigindices),1),'*','Color',[0,1,1],'FontSize',15,'HorizontalAlignment','center'); %display significance star
%         end
        
        % Display number of trials per average
%         Ydim = get(gca,'YLim');
%         for s = 1:numStimuli
%             text(s,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%ROIdata',ROIdata.rois(r).nTrials(s)),'HorizontalAlignment','center');
%         end
        
        % Display p-values
%         Ydim = get(gca,'YLim');
%         for s = 1:numStimuli
%             text(s,Ydim(1)+(Ydim(2)-Ydim(1))/20,sprintf('p=%.3f',ROIdata.rois(r).PValue(s)),'HorizontalAlignment','center');
%         end
        
        hold off      
        
        % Save plot to PDF
        if saveToPDF
            if ischar(ROIdata)
                title(sprintf('ROI: %d, Label: %s, File: %s',rindex,ROIdata.rois(rindex).label{1},ROIdata));
                filename = strcat(ROIdata(1:end-3), 'pdf');
            else
                title(sprintf('ROI: %d, Label: %s', rindex, ROIdata.rois(rindex).label{1}));
                filename = 'TuningCurves.pdf';
            end
            export_fig(hf, filename,'-append'); close(hf);
        end
        
    end % plot figure
    
end %cycle ROIs


if saveToROIFile && ischar(ROIdata)
    for rindex = 1:numROIs
        if ~isempty(ROIdata.rois(rindex).label) && strcmp(ROIdata.rois(rindex).label{1}, 'none')
            ROIdata.rois(rindex).label = {};
        end
    end
    save(ROIdata, 'ROIdata', '-append', '-mat');
end

function [LinDiff,CI] = LinDiff(Raw,StimLog,varargin)

N = 10000; % # of bootstraps
alpha = 0.05;
verbose = false;
labels = {};

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'N'
                N = varargin{index+1};
                index = index + 2;
            case 'alpha'
                alpha = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = ~verbose;
                index = index + 1;
            case {'Labels','labels'}
                labels = varargin{index+1};
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


%% Compute linear difference
numROIs = size(Raw,1);
numW = sum(StimLog,2);
StimLog = cellfun(@find, mat2cell(StimLog,ones(size(StimLog,1),1),size(StimLog,2)),'UniformOutput',false);
Index = find(numW>1);
numMult = sum(numW>1);
LinDiff = nan(numROIs,numMult);
CI = nan(numROIs,2,numMult);
if verbose; p = parfor_progress(numMult*numROIs); end
if N~=0
    parfor s = 1:numMult
        if numW(Index(s))==2
            func = @(x,y,z) z-(x+y);
        elseif numW(Index(s))==3
            func = @(w,x,y,z) z-(w+x+y);
        elseif numW(Index(s))==4
            func = @(v,w,x,y,z) z-(v+w+x+y);
        elseif numW(Index(s))==5
            func = @(u,v,w,x,y,z) z-(u+v+w+x+y);
        end
        for r = 1:numROIs
            [LinDiff(r,s),CI(r,:,s)] = BootStrap(N,func,[Raw(r,StimLog{Index(s)}+1),Raw(r,Index(s))],'alpha',alpha);
            %             if verbose; parfor_progress(p); end
        end
    end
else
    Curves = cellfun(@nanmean, Raw);
    parfor s = 1:numMult
        LinDiff(:,s) = Curves(:,Index(s)) - sum(Curves(:,StimLog{Index(s)}+1),2);
    end
end
if verbose; parfor_progress(p,0); end


%% Plot results
if verbose
    figure; 
    plot([.5,numMult+.5],[0,0],'k--');
    hold on;
    for r = 1:numROIs
        errorbar(1:numMult,LinDiff(r,:),CI(r,:,1),CI(r,:,2));
    end
    xlim([.5,numMult+.5]);
    if ~isempty(labels)
        set(gca,'XTick',1:numMult,'XTickLabel',labels,'XTickLabelRotation',90);
    end
    ylabel('Difference from Linear Sum (dF/F)')
    xlabel('Stimulus');
end


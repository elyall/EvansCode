function [val,ci,bootdist] = BootStrap(N,bootfunc,data,varargin)

alpha = .05;
verbose = false;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'alpha'
                alpha = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = ~verbose;
                index = index + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end


%% Compute value
Mean = cellfun(@mean,data,'UniformOutput',false);
val = feval(bootfunc,Mean{:});

%% Create bootstrap distribution
numVars = numel(data);
num = cellfun(@numel,data);
bootdist = nan(N,1);
for n = 1:N
    current = arrayfun(@(x) datasample(data{x},num(x)),1:numVars,'UniformOutput',false);
    Mean = cellfun(@mean,current,'UniformOutput',false);
    bootdist(n) = feval(bootfunc,Mean{:});
end

%% Compute confidence intervals
pct1 = 100*alpha/2;
pct2 = 100-pct1;
lower = prctile(bootdist,pct1,1); 
upper = prctile(bootdist,pct2,1);
ci =[lower;upper]; % return

%% Plot output
if verbose
    figure; histogram(bootdist,'Normalization','probability');
    ylabel('Fraction of Boot Distribution'); xlabel('Value');
    YLim = ylim(gca);
    hold on;
    h1 = plot([val,val],YLim,'r');
    h2 = plot([ci(1),ci(1)],YLim,'k');
    plot([ci(2),ci(2)],YLim,'k');
    legend([h1,h2],'observed',sprintf('%d%% C.I.',round((1-alpha)*100)));
end
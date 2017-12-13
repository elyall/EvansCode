function [val,ci,bootdist] = BootStrap(N,bootfunc,data,varargin)
% BOOTSTRAP runs bootstrp on multiple random variables that are not sampled
% at the same point in time.
%
% VAL = BOOTSTRAP(N,BOOTFUNC,DATA) computes the mean of each cell of DATA
% and then plugs those means into BOOTFUNC to produce the output: VAL. DATA
% is a vector cell array with each cell containing a vector of raw data.
% BOOTFUNC is an anonymous function whose number of input variables is
% equal to numel(DATA).
%
% [VAL,CI,BOOTDIST] = BOOTSTRAP(N,BOOTFUNC,DATA) bootstraps the sample mean
% of each cell of DATA to generate N  means for each dataset, and then
% plugs those sample means into BOOTFUNC to generate BOOTDIST, a
% bootstrapped distribution of length N of the output value VAL. The 95th
% percent confidence intervals of the bootstrap distribution are returned
% as CI.
%
% [VAL,CI,BOOTDIST] = BOOTSTRAP(...,'alpha',ALPHA) computes the
% 100*(1-ALPHA) percent bootstrap confidence interval of the statistic
% defined by the function BOOTFUN. ALPHA is a scalar between 0 and 1. The
% default value of ALPHA is 0.05.
%
% [VAL,CI,BOOTDIST] = BOOTSTRAP(...,'verbose') plots a histogram of the
% bootstrap distribution and overlays VAL and CI.


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

if nargout > 1
    %% Create bootstrap distribution
    numVars = numel(data);
    num = cellfun(@numel,data);
    bootdist = nan(N,1);
    parfor n = 1:N
        current = arrayfun(@(x) datasample(data{x},num(x)),1:numVars,'UniformOutput',false);
        Mean = cellfun(@mean,current,'UniformOutput',false);
        if ~isequal(bootfunc,false)
            bootdist(n) = feval(bootfunc,Mean{:});
        else
            bootdist(n) = Mean{1};
        end
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
end
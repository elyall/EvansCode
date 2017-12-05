function [result,Var,FuncCov,Cov] = linearFuncOnRV(data,func,varargin)
% LINEARFUNCONRV - computes linear function(s) on a set of random variables
% and propagates error

% data is n x m with n observations and m variables; or a cell array 1xm
% where each cell contains an arbitrary # of observations

% func is k x n matrix of combination coefficients (e.g. if subtracting one
% variable from another variable then it would be: [1,-1]). k is # of
% linear functions to compute.

type = 'uncorrelated'; % 'correlated' or 'uncorrelated' specifying whether the observations are correlated or not

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'type'
                type = varargin{index+1};
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


%% Compute output and propagate error
if iscell(data) % assumes data is uncorrelated since it's in a cell
    N = length(data);                    % # of variables
    Mean = cellfun(@mean,data);          % compute means of variables
    Cov = zeros(N,N);                    % create variance-covariance matrix
    Cov(1:N+1:end) = cellfun(@var,data); % set diagonal to variances (covariances are all 0)
    if isrow(Mean)
        Mean = Mean';                    % ensure column vector for computing result
    end
else
    N = size(data,2);                    % # of variables
    Mean = mean(data)';                  % compute means of variables
    switch type
        case 'correlated'
            Cov = cov(data);             % compute variance-covariance matrix
        case 'uncorrelated'
            Cov = zeros(N,N);            % create variance-covariance matrix
            Cov(1:N+1:end) = var(data);  % set diagonal to variances (covariances are all 0)
    end
end

result = func*Mean;       % compute outputs of linear functions
FuncCov = func*Cov*func'; % compute functions' variance-covariance matrix
Var = diag(FuncCov);      % pull out variance of each function


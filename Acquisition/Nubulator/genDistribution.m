function out = genDistribution(N,varargin)

Distribution = 'specific'; % uniform, normal, specific, constant
Weights = {[.1,.6,.9,1,1,.9,.6,.1],[1,1]}; % values between 0 to 1, numel>=2

%% Check input arguments
if ~exist('N','var') || isempty(N)
    N = 100;
end

if isempty(varargin)
    varargin = {[0,1]};
end

if ~iscell(Weights)
    Weights = {Weights};
end


%% Generate points
numDims = numel(varargin);
switch Distribution
    case 'uniform'
        out = rand(N,numDims);
    case 'normal'
        out = randn(N,numDims);
    case 'specific'
        if numel(Weights)==1 && numDims>1
            Weights = repmat(Weights,1,numDims);
        end
        tempDist = cell(numDims,1);
        tempWeights = cell(numDims,1);
        for dindex = 1:numDims
            if numel(varargin{dindex})>2
                tempDist{dindex} = 0:1/(varargin{dindex}(3)-1):1;
            else
                tempDist{dindex} = 0:1/10000:1;
            end
            Xw = 0:1/(numel(Weights{dindex})-1):1;
            tempWeights{dindex} = interp1(Xw, Weights{dindex}, tempDist{dindex});
        end
        if numDims == 1
            out = tempDist{1};
            Weights = tempWeights{1};
        else
            [X,Y] = meshgrid(tempDist{1},tempDist{2});
            out = [X(:),Y(:)];
            [X,Y] = meshgrid(tempWeights{1},tempWeights{2});
            Weights = prod([X(:),Y(:)],2);
        end
        if N > size(out,1)
            out = datasample(out, N, 1, 'Replace', true, 'Weights', Weights);
        else
            out = datasample(out, N, 1, 'Replace', false, 'Weights', Weights);
        end
    case 'constant'
        if numDims == 1
            out = linspace(0,1,N)';
        else %if numDims > 1
            temp = repmat(N^(1/numDims), 1, numDims);
            for dindex = 1:numDims
                if numel(varargin{dindex})>2
                    temp(dindex) = varargin{dindex}(3);
                end
            end
            if numDims == 2
                [X,Y] = meshgrid(linspace(0,1,temp(1)),linspace(0,1,temp(2)));
                out = [X(:),Y(:)];
            elseif numDims == 3
                [X,Y,Z] = meshgrid(linspace(0,1,temp(1)),linspace(0,1,temp(2)),linspace(0,1,temp(3)));
                out = [X(:),Y(:),Z(:)];
            end
        end
end

% % Scale to be between 0 and 1 (guarantees a point on each extreme)
% out = bsxfun(@minus, out, min(out));
% out = bsxfun(@rdivide, out, max(out));

% Scale to requested region
for dindex = 1:numDims
    out(:,dindex) = out(:,dindex)*(varargin{dindex}(2)-varargin{dindex}(1)) + varargin{dindex}(1);
end

% Fill in to produce the desired number of outputs
n = size(out,1);
if n < N
    T = floor(N/n);
    out = cat(1, repmat(out,T,1), datasample(out, rem(N,n), 1, 'Replace', true));
end


%% Check result
nbins = 20;
figure;
if numDims == 1
    subplot(1,2,1); plot(out,'k.'); axis tight; ylabel('Value'); xlim([varargin{1}(1),varargin{1}(2)]);
    subplot(1,2,2); histogram(out,nbins); axis tight; xlabel('X Position'); ylabel('count'); xlim([varargin{1}(1),varargin{1}(2)]);
elseif numDims == 2
    subplot(1,3,1); plot(out(:,1),out(:,2),'k.'); axis tight; xlabel('X'); ylabel('Y');
    xlim([varargin{1}(1),varargin{1}(2)]); ylim([varargin{2}(1),varargin{2}(2)]);
    subplot(1,3,2); histogram(out(:,1),nbins); axis tight; xlabel('X Position'); ylabel('count'); xlim([varargin{1}(1),varargin{1}(2)]);
    subplot(1,3,3); histogram(out(:,2),nbins); axis tight; xlabel('Y Position'); ylabel('count'); xlim([varargin{2}(1),varargin{2}(2)]);
end

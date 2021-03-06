function [p,nuldist,actual] = permutationTest(Data, GroupID, varargin)

numPerms = 1000;
func = @(x,y) mean(x)-mean(y); % operates along 1st dimension

saveOut = false;
saveFile = '';


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numPerms'
                numPerms = varargin{index+1};
                index = index + 2;
            case 'func'
                func = varargin{index+1};
                index = index + 2;
            case {'save','Save'}
                saveOut = true;
                index = index + 1;
            case 'saveFile'
                saveFile = varargin{index+1};
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

if iscell(Data)
    if ~exist('GroupID','var') || isempty(GroupID)
        GroupID = repelem(1:numel(Data),cellfun(@numel,Data))';
    end
    try
        Data = cat(1,Data{:});
    catch
        Data = cat(2,Data{:})';
    end
elseif isrow(Data)
    Data = Data';
end
if islogical(GroupID) && size(GroupID,2)>1
    temp = GroupID;
    GroupID = zeros(size(temp,1),1);
    for ind = 1:size(temp,2)
        GroupID(temp(:,ind)) = ind;
    end
end
N = numel(Data);
[Groups,~,GroupID] = unique(GroupID); % ensure GroupID is 1:N
if numel(Groups)~=2
    error('Can only do permutation test with 2 groups!');
end


%% Calculate actual value
actual = func(Data(GroupID==1),Data(GroupID==2));


%% Calculate null distribution

% try
    
    % Create matrix of different permutations
    Perms = zeros(N,numPerms);
    for n = 1:numPerms
        Perms(:,n) = GroupID(randperm(N));
    end
    
    % Compute null distribution based off permutations
    Data = repmat(Data,1,numPerms);
    nuldist = func(reshape(Data(Perms==1),sum(GroupID==1),numPerms),reshape(Data(Perms==2),sum(GroupID==2),numPerms));
    
% catch % likely matrix is too big? (too many permutations or data points)
%     warning('Failed at fast method, using slower method. (likely too many permutations or data points to fit whole matrix in memory)');
%     
%     % Initialize output
%     nuldist = nan(numPerms,1);
%     
%     % Cycle through permutations
%     % fprintf('Computing %d permutations...\n',numPerms);
%     % pfH = parfor_progress(numPerms);
%     parfor pindex = 1:numPerms
%         
%         % Shuffle data
%         currentID = GroupID(randperm(N));
%         
%         % Compute value
%         nuldist(pindex) = func(Data(currentID==1),Data(currentID==2));
%         
%         %     parfor_progress(pfH); % update status
%     end
%     % parfor_progress(pfH,0);
    
% end


%% Calculate p-value
Mean = mean(nuldist); % compute mean of null distribution
tempnul = nuldist-Mean; % mean-center null distribution
tempactual = actual-Mean; % rectify actual with distribution
p = sum(abs(tempnul)>abs(tempactual))/numPerms; % compute two-sided p value(s)


%% Save outputs
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-v7.3');
    else
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-append');
    end
    fprintf('Saved permutation results to: %s\n',saveFile);
end


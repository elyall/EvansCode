function out = movprctile(data,prc,sz,dim)

if ~exist('dim','var')
    dim=1;
end

if dim==2
    data = data';
end

% Determine parameters
if numel(sz)==1
    if ~mod(sz,2)
        error('size must be odd or 1x2');
    end
    window = repmat(floor(sz/2),[1,2]);
elseif numel(sz)==2
    window = abs(sz);
end
inds = nan(sum(window+1),1);
buffer = nan(sum(window)+1,size(data,2));

% Compute moving prctile
out = nan(size(data));
parfor r = 1:size(data,1)
    out(r,:) = prctile(data(max(1,r-window(1)):min(end,r+window(2)),:),prc);
end

% buffer(1:window(2),:) = data(1:window(2),:);
% inds = (sum(window)+1)*prc/100;
% for r = 1:size(data,1)
%   if r>sum(window)+1
%       bad = inds<r-window(1);
%       buffer(bad,:) = [];
%       inds(bad) = [];
% end

if dim==2
    out = out';
end
function [p,Label] = calcPValues(Data,GroupID,Names,test)
% Unpaired two-sample significance testing
%
% Data - vector of data points, or cell array of data points where each
% cell is a unique group
%
% GroupID - (ignored if Data is cell array) vector of length equal to Data
% vector specifying which group each data point is a part of, or logical
% matrix of height equal to length of Data vector and width equal to number
% of groups specifying which data point is a part of which group.
% (GroupID==0 is ignored)
%
% Names - cell array of strings equal in length to the number of groups
%
% test - 'ranksum' or 'ttest2' specifying which test to use

MultipleComparisons_correction = true;

if ~exist('test','var') || isempty(test)
    test = 'ranksum';
end

if iscell(Data)
    GroupID = repelem(1:numel(Data),cellfun(@numel,Data));
    try 
        Data = cat(1,Data{:});
    catch
        Data = cat(2,Data{:});
    end
elseif islogical(GroupID) && size(GroupID,2)>1
    temp = GroupID;
    GroupID = zeros(size(temp,1),1);
    for ind = 1:size(temp,2)
        GroupID(temp(:,ind)) = ind;
    end
end

Groups = unique(GroupID);
Groups(Groups==0) = [];
combs = nchoosek(Groups,2);

if ~exist('Names','var')
    Names = {};
end
if isempty(Names) && numel(Groups)>2
    Names = strcat('G',cellstr(num2str((1:max(Groups))')));
end

p = nan(1,size(combs,1));
Label = cell(size(combs,1),1);
for c = 1:size(combs,1)
    switch test
        case {'ranksum'}
            p(c) = ranksum(Data(GroupID==combs(c,1)),Data(GroupID==combs(c,2)));
        case {'t','ttest','ttest2','t-test'}
            p(c) = ttest2(Data(GroupID==combs(c,1)),Data(GroupID==combs(c,2)));
        case {'ks','kurskalwallis'}
            p(c) = kruskalwallis(Data,GroupID);
    end
end

if numel(p)>1 && MultipleComparisons_correction==true
   [~,~,p] = fdr_bh(p);
end

for c = 1:length(p)
    if p(c)==0
        L = 'p= 0';
    elseif p(c)<.01
        L = sprintf('p=%.1e',p(c));
    else
        L = sprintf('p=%.3f',p(c));
    end
    if isempty(Names)
        Label{c} = L;
    else
        Label{c} = sprintf('%s:%s %s',Names{combs(c,1)},Names{combs(c,2)},L);
    end
end

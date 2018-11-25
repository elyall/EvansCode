function [p,Label] = calcPValues(Data,GroupID,Names)

test = 'Wilcoxon';

if islogical(GroupID) && size(GroupID,2)>1
    temp = GroupID;
    GroupID = zeros(size(temp,1),1);
    for ind = 1:size(temp,2)
        GroupID(temp(:,ind)) = ind;
    end
end

Groups = unique(GroupID);
Groups(Groups==0) = [];
combs = nchoosek(Groups,2);

if ~exist('Names','var') || isempty(Names)
    Names = strcat('G',cellstr(num2str(Groups)));
end

p = nan(1,size(combs,1));
Label = cell(size(combs,1),1);
for c = 1:size(combs,1)
    switch test
        case {'Wilcoxon','ranksum'}
            p(c) = ranksum(Data(GroupID==combs(c,1)),Data(GroupID==combs(c,2)));
        case {'t','ttest','t-test'}
            p(c) = ttest2(Data(GroupID==combs(c,1)),Data(GroupID==combs(c,2)));
    end
    if p(c)<.01
        Label{c} = sprintf('%s:%s p=%.1e',Names{combs(c,1)},Names{combs(c,2)},p(c));
    else
        Label{c} = sprintf('%s:%s p=%.2f',Names{combs(c,1)},Names{combs(c,2)},p(c));
    end
end

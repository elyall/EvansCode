function [Index,Actual,Sum,Indiv] = LinearityIndex(Curves,Stim)
% Curves - ROIs x Stim
% Stim - Stim x Cond

verbose = false;

CondIndex = find(sum(Stim,2)==1);     % determine row location of single stimuli
[~,X] = find(Stim(sum(Stim,2)==1,:)); % determine column location of single stimuli
[~,order] = sort(X,'ascend');         % determine column-wise order of stimuli
CondIndex = CondIndex(order);         % order row indices of single stimuli by which column it's in (important for indexing later)

MultIndices = find(sum(Stim,2)>1); % determine stimuli which have multiple conditions

N = numel(MultIndices);
numROIs = size(Curves,1);
Indiv = cell(numROIs,N);
Sum = zeros(numROIs,N);
Actual = zeros(numROIs,N);
parfor index = 1:numROIs
    for c = 1:N
        Indiv{index,c} = Curves(index,CondIndex(Stim(MultIndices(c),:)));
        Sum(index,c) = sum(Curves(index,CondIndex(Stim(MultIndices(c),:))));
        Actual(index,c) = Curves(index,MultIndices(c));
    end
end
Index = (Actual-Sum)./Sum;

if verbose
    figure;
    scatter(Sum(:),Actual(:),'.');
    xlabel('Sum (dF/F)');
    ylabel('Actual (dF/F)');
    XLim = get(gca,'XLim');
    YLim = get(gca,'YLim');
    Lim = [min(XLim(1),YLim(1)),max(XLim(2),YLim(2))];
    hold on;
    plot(Lim,Lim,'--');
end
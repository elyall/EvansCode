function h = plotmatrix_mult(Data,GroupID,Names,Colors,GroupNames)

numBins = 20;
Size = 20;
showP = true;
AxisEqual = true;
UnityLine = true;
XLim = [0,1];
YLim = [0,1];
Corr = false;
Median = false;

if ~exist('GroupID','var') || isempty(GroupID)
    plotmatrix(Data);
    return
end

N = size(Data,2);
if ~exist('Names','var') || isempty(Names)
    Names = strcat('Var',cellstr(num2str((1:N)')));
end

[IDs,~,GroupID] = unique(GroupID);
numGroups = numel(IDs);
if ~exist('Colors','var') || isempty(Colors)
    Colors = jet(numGroups);
end


%%

if isrow(GroupID)
    GroupID = GroupID';
end
Index = GroupID==1:numGroups;

edges = [min(Data);max(Data)];
edges = arrayfun(@(x) linspace(edges(1,x),edges(2,x),numBins+1), 1:N, 'UniformOutput',false);
edges = cat(1,edges{:});

h = nan(N);
figure;
for x = 1:N
    for y = 1:x
        if showP && numGroups>1
            ind = sub2ind([N,N+1],x,y);
            h(y,x) = subplot(N,N+1,ind+y); hold on;
        else
            ind = sub2ind([N,N],x,y);
            h(y,x) = subplot(N,N,ind); hold on;
        end
        if Corr
            rho = nan(numGroups,1);
            p = nan(numGroups,1);
        end
        for g = 1:size(Index,2)
            if x==y
                histogram(Data(Index(:,g),x),edges(x,:),'EdgeColor',Colors(g,:),'DisplayStyle','stairs','Normalization','probability');
                xlabel(Names{x});
                ylabel('Fraction');
                
            else
                scatter(Data(Index(:,g),x),Data(Index(:,g),y),Size,Colors(g,:),'.',...
                    'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
                if ~isempty(YLim)
                    ylim(YLim);
                end
                if ~isempty(XLim)
                    xlim(XLim);
                end
                if UnityLine
                    X = get(gca,'XLim');
                    Y = get(gca,'YLim');
                    Lim = [min(X(1),Y(1)),max(X(2),Y(2))];
                    hold on;
                    plot(Lim,Lim,'r--');
                end
                if AxisEqual
                    axis equal;
                end
                if Corr
                    [rho(g),p(g)] = corr(Data(Index(:,g),x),Data(Index(:,g),y));
                end
                xlabel(Names{x});
                ylabel(Names{y});
            end
        end
        
        if x==y
            
            % Show legend in last plot
            if x==N && exist('GroupNames','var') && ~isempty(GroupNames)
                legend(GroupNames,'Location','best'); legend boxoff;
            end
            
            if Median
                Y = get(gca,'YLim');
                hold on;
                for g = 1:numGroups
                    plot(repmat(median(Data(Index(:,g),x)),1,2),Y,'--','Color',Colors(g,:));
                end
            end
            
            if showP && numGroups>1
                subplot(N,N+1,ind+y-1);
                combs = nchoosek(1:numGroups,2);
                p = nan(1,size(combs,1));
                temp = cell(size(combs,1),1);
                for c = 1:size(combs,1)
                    p(c) = ranksum(Data(Index(:,combs(c,1)),x),Data(Index(:,combs(c,2)),x));
                    if p(c)<.01
                        temp{c} = sprintf('%s:%s p=%.1e',GroupNames{combs(c,1)},GroupNames{combs(c,2)},p(c));
                    else
                        temp{c} = sprintf('%s:%s p=%.2f',GroupNames{combs(c,1)},GroupNames{combs(c,2)},p(c));
                    end
                end
                text(.2,.5,temp,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','middle');
                axis off;
            end
   
        elseif x~=y
            
            if Corr
                temp = cell(numGroups,1);
                for g = 1:numGroups
                    if p(g)<.01
                        temp{g} = sprintf('%s rho=%.2f (p=%.1e)',GroupNames{g},rho(g),p(g));
                    else
                        temp{g} = sprintf('%s rho=%.2f (p=%.2f)',GroupNames{g},rho(g),p(g));
                    end
                end
                text(0,1,temp,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top');
            end            
        end
            
    end
end


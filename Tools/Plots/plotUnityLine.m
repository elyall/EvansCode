function h = plotUnityLine(hA)

if ~exist('h','var') || isempty(hA)
    hA = gca;
end

XLim = get(gca,'XLim');
YLim = get(gca,'YLim');
Lim = [min(XLim(1),YLim(1)),max(XLim(2),YLim(2))];

axes(hA);
Hold = ishold(hA);
hold on;
h = plot(Lim,Lim,'k--');
if ~Hold
    hold off
end
temp = get(hA,'Children');
set(hA,'Children',temp([2:end,1]));
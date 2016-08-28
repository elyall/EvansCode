function [PreCoM, PosCoM] = statSignificance(PreROIdata, PreTrialIndex, PosROIdata, PosTrialIndex, RealCoM, RealSel, ROIindex)

% verbose = true;

%% SPLIT
% is change in stats during either half significantly different from that
% observerd with trimming

[PreCoM, PreSel] = permuteData(PreROIdata,'TrialIndex',PreTrialIndex,'type','split','ROIindex',ROIindex,'DistBtwn',2.5);
[PosCoM, PosSel] = permuteData(PosROIdata,'TrialIndex',PosTrialIndex,'type','split','ROIindex',ROIindex,'DistBtwn',2.5);

%% CoM
p = nan(4,1);

figure('Name','Center of Mass','Position',[50, 50, 1400, 800]);

% Difference in first half vs second half of full pad
subplot(2,2,1);
histogram(diff(PreCoM,[],2));
p(1)=signrank(diff(PreCoM,[],2));
title('Full Pad: difference in halves');
xlabel('Change in CoM');
ylabel('# of ROIs');

% Difference in first half vs second half of single whisker
subplot(2,2,2);
histogram(diff(PosCoM,[],2));
p(2)=signrank(diff(PosCoM,[],2));
title('Single Whisker: difference in halves');
xlabel('Change in CoM');
ylabel('# of ROIs');

subplot(2,2,3);
h=histogram(diff(PreCoM,[],2));
hold on;
histogram(diff(RealCoM,[],2),'BinWidth',h.BinWidth);
p(3)=signrank(diff(PreCoM,[],2),diff(RealCoM,[],2));
title('Full Pad halves vs. Actual');
xlabel('Change in CoM');
ylabel('# of ROIs');
legend({'Full Pad halves','Trimming'},'Location','Best');

subplot(2,2,4);
h=histogram(diff(PosCoM,[],2));
hold on;
histogram(diff(RealCoM,[],2),'BinWidth',h.BinWidth);
p(4)=signrank(diff(PosCoM,[],2),diff(RealCoM,[],2));
title('Single Whisker halves vs. Actual');
xlabel('Change in CoM');
ylabel('# of ROIs');
legend({'Single Whisker halves','Trimming'},'Location','Best');

XLim = [inf,-inf];
YLim = [inf,-inf];
for index = 1:4
    subplot(2,2,index);
    temp = get(gca,'XLim');
    XLim = [min(XLim(1),temp(1)),max(XLim(2),temp(2))];
    temp = get(gca,'YLim');
    YLim = [min(YLim(1),temp(1)),max(YLim(2),temp(2))];
end
for index = 1:4
    subplot(2,2,index);
    xlim(XLim);
    ylim(YLim);
    text(XLim(1)+range(XLim)/20, YLim(2)-range(YLim)/20, sprintf('wilcoxon\np=%.1e', p(index)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end


%% Sel
p = nan(4,1);

figure('Name','Selectivity','Position',[50, 50, 1400, 800]);

% Difference in first half vs second half of full pad
subplot(2,2,1);
histogram(diff(PreSel,[],2));
p(1)=signrank(diff(PreSel,[],2));
title('Full Pad: difference in halves');
xlabel('Change in Selectivity');
ylabel('# of ROIs');

% Difference in first half vs second half of single whisker
subplot(2,2,2);
histogram(diff(PosSel,[],2));
p(2)=signrank(diff(PosSel,[],2));
title('Single Whisker: difference in halves');
xlabel('Change in Selectivity');
ylabel('# of ROIs');

subplot(2,2,3);
h=histogram(diff(PreSel,[],2));
hold on;
histogram(diff(RealSel,[],2),'BinWidth',h.BinWidth);
p(3)=signrank(diff(PreSel,[],2),diff(RealSel,[],2));
title('Full Pad halves vs. Actual');
xlabel('Change in Selectivity');
ylabel('# of ROIs');
legend({'Full Pad halves','Trimming'},'Location','NorthEast');

subplot(2,2,4);
h=histogram(diff(PosSel,[],2));
hold on;
histogram(diff(RealSel,[],2),'BinWidth',h.BinWidth,'Normalization','probability');
p(4)=signrank(diff(PosSel,[],2),diff(RealSel,[],2));
title('Single Whisker halves vs. Actual');
xlabel('Change in Selectivity');
ylabel('# of ROIs');
legend({'Single Whisker halves','Trimming'},'Location','NorthEast');

XLim = [inf,-inf];
YLim = [inf,-inf];
for index = 1:4
    subplot(2,2,index);
    temp = get(gca,'XLim');
    XLim = [min(XLim(1),temp(1)),max(XLim(2),temp(2))];
    temp = get(gca,'YLim');
    YLim = [min(YLim(1),temp(1)),max(YLim(2),temp(2))];
end
for index = 1:4
    subplot(2,2,index);
    xlim(XLim);
    ylim(YLim);
    text(XLim(1)+range(XLim)/20, YLim(2)-range(YLim)/20, sprintf('wilcoxon\np=%.1e', p(index)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end


%% PERMUTATION
% Where do stats for each half fall within permuted distribution of
% respected half as well as other half

% PercCoM = nan(numROIs,4);
% PercSel = nan(numROIs,4);
function removeBlankMasks(SegmentFiles)

if ischar(SegmentFiles)
    SegmentFiles = {SegmentFiles};
end
numFiles = numel(SegmentFiles);

for findex = 1:numFiles
    load(SegmentFiles{findex},'mask','-mat');
    badIndex = sum(mask)==0;
    if any(badIndex)
        mask(:,badIndex) = [];
        save(SegmentFiles{findex},'mask','-mat','-append');
        fprintf('fixed\n');
    end
end
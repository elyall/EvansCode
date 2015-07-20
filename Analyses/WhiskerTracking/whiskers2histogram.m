function whiksers2histogram(ImageFiles, TrialIndex, StimIndex, FrameIndex, TrialsToAnalyze)
%WHISKERS2HISTOGRAM Convert high speed videos of whiskers to 2d histograms.
%   Have the user enter details of the experiment, and specify where the video
%   files are located, and which .dat file corresponds to those videos. This
%   script will then measure the whisker density for running trials during the
%   specified analysis period.
%
%   Saves a 3-d matrix in the local directory. The third dimension corresponds
%   to each trial type. To visualize the density map you can use imagesc like
%   so: figure; imagesc(whisker_density(:,:,1), [0 0.2]), colormap hot
%
%   Dependancies: Norpix2MATLABopenSingleFrame, Norpix2MATLAB,
%   classify_run_trials, MakeGaussWindow, fwhm.
%
%   UC Berkeley
%   Adesnik Lab
%   G. Telian
%   20150624

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory = cd;
saveOut = false;
saveFile = '';
ROI = true;

if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles, p] = uigetfile('*.tif;*.avi;*.seq;*.raw', 'Select files to analyze', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles)
        return
    elseif iscell(ImageFiles)
        for findex = 1:numel(ImageFiles)
            ImageFiles{findex} = fullfile(p, ImageFiles{findex});
        end
    else
        ImageFiles = {fullfile(p, ImageFiles)};
    end
end
numFiles = numel(ImageFiles);

if ~exist('TrialIndex', 'var') || isempty(TrialIndex)
    TrialIndex = ones(numFiles, 1);
end
numTrials = numel(unique(TrialIndex));

if ~exist('StimIndex', 'var') || isempty(StimIndex)
    StimIndex = ones(numTrials, 1);
end
StimIDs = unique(StimIndex);
numStimuli = numel(StimIDs);

if ~exist('FrameIndex', 'var') || isempty(FrameIndex)
    FrameIndex = repmat([1,inf], numStimuli, 1);
end

if ~exist('TrialsToAnalyze', 'var') || isempty(TrialsToAnalyze)
    TrialsToAnalyze = 1:numTrials;
end

if saveOut && isempty(saveFile)
    index = 0;
    saveFile = [ImageFiles{1}(1:end-4), '.den'];
    while exist(saveFile, 'file')
        index = index + 1;
        saveFile = [ImageFiles{1}(1:end-4), num2str(index), '.den'];
    end
    fprintf('Will save whisker density info to %s\n', saveFile);
end

%% Make whisker density matrices for each stimulus condition

% Load in one movie frame to get image size for array pre-allocation
img = readFlea3(ImageFiles(1));
[H, W, ~] = size(img);
whisker_density = zeros(H, W, numStimuli);

% Select ROI to analyze
if ROI
    hF = figure('Name', 'Choose ROI to analyze');
    imagesc(img);
    mask = roipoly;
    close(hF);
end

% Calculate whisker density for each stimulus
for sindex = 1:numStimuli
    goodTrials = intersect(find(StimIndex == StimIDs(sindex)), TrialsToAnalyze);

    % Iterate through all movies of a particular trial type
    frame_counter = 0;
    for tindex = goodTrials

        % load trial
        img = [];
        currentFiles = ImageFiles(TrialIndex==tindex);
        img = readFlea3(currentFiles);
%         if ~strcmp(class(img), 'uint8')
%             img = uint8(img);
%         end
        numFrames = size(img, 3);
        
        % Remove data outside ROI
        if ROI
            img(repmat(~mask, 1, 1, numFrames)) = 0;
        end

        % Compute threshold
        % used to binarize image and separate reflective whisker from
        % background and noise.
        if ~ROI
            threshold = mean(img(:)) + 4*std(img(:));
        else
            threshold = mean(img(mask)) + 4*std(img(mask));
        end

        % Determine period of trial to analyze
        if FrameIndex(sindex, 2) == inf
            lastframe = numFrames;
        else
            lastframe = FrameIndex(sindex, 2);
        end
        numFrames = lastframe - FrameIndex(sindex, 1) + 1;
        
        % Compute whisker counts for density calculation
        img_thresh = img(:,:,FrameIndex(sindex, 1):lastframe)>threshold;
        whisker_density(:, :, sindex) = whisker_density(:, :, sindex) + sum(img_thresh, 3);
        frame_counter = frame_counter + numFrames;
        
        % visualization
        figure('Position', [50, 50, 1400, 800]);
        for findex = 1:numFrames
            subplot(1,2,1);
            imagesc(img(:,:,findex));
            colormap(gray);
            axis off
            subplot(1,2,2);
            imagesc(img_thresh(:,:,findex));
            axis off
            drawnow
            pause(0.05);
        end
    end
    
    % End Compute whisker counts for density calculation
    whisker_density(:, :, sindex) = whisker_density(:, :, sindex)/frame_counter;

end

if saveOut
    save(saveFile, 'ImageFiles', 'whisker_density', '-v7.3');
end








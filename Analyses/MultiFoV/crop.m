function [Images, Maps] = crop(Images, rect, Maps)


%% Parse input arguments
if iscell(Images)
    output = 'cell';
elseif ~isempty(Images)
    Images = {Images};
    output = 'numeric';
else
    output = 'empty';
end

if ~exist('rect', 'var') || isempty(rect)
    rect = true;
%     rect = repmat([32.51, 0, 729.98, 512], numFiles, 1);
end

if ~exist('Maps', 'var') || isempty(Maps)
    Maps = [];
end


%% Crop images
if ~isempty(Images) 
    numFiles = numel(Images);
    
    % Fix input
    if ~islogical(rect) && size(rect, 1) == 1 && numFiles>1
        rect = repmat(rect, numFiles, 1);
    end
    
    % Crop images
    for findex = 1:numFiles
        
        % UI-crop
        if islogical(rect) && rect == true
            [~, rect(findex,:)] = imcrop(Images{findex}(:,:,1,1,1));
        end
        
        % Make rect integers
        rect(findex,[1,2]) = floor(rect(findex,[1,2]));
        rect(findex,[3,4]) = ceil(rect(findex,[3,4]));
        
        % Crop image
        Images{findex} = Images{findex}(rect(findex,2)+1:rect(findex,2)+rect(findex,4), rect(findex,1)+1:rect(findex,1)+rect(findex,3),:,:,:);
        
    end
    
end

% Close figure (UI-crop only)
if islogical(rect) && rect == true
    close gcf;
end

% Fix output
if strcmp(output, 'numeric')
    Images = Images{1};
end


%% Crop maps
if ~isempty(Maps)
    numFiles = numel(Maps);
    
    % Fix input
    if ~islogical(rect) && size(rect, 1) == 1 && numFiles>1
        rect = repmat(rect, numFiles, 1);
    end
    
    for findex = 1:numFiles
        
        % Make rect integers
        rect(findex,[1,2]) = floor(rect(findex,[1,2]));
        rect(findex,[3,4]) = ceil(rect(findex,[3,4]));
        
        % Crop map
        Maps(findex).XWorldLimits(1) = Maps(findex).XWorldLimits(1)+rect(findex,1);
        Maps(findex).XWorldLimits(2) = Maps(findex).XWorldLimits(1)+rect(findex,3);
        Maps(findex).YWorldLimits(1) = Maps(findex).YWorldLimits(1)+rect(findex,2);
        Maps(findex).YWorldLimits(2) = Maps(findex).YWorldLimits(1)+rect(findex,4);
        Maps(findex).ImageSize = rect(findex,[4,3]);
    end
end




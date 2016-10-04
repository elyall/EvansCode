%% FIX MappedTensor Subtraction FUCK UP

% Files = {'/home/elyall/Documents/Data/1955/150227/1955_110_001.sbx',...
%          '/home/elyall/Documents/Data/1955/150227/1955_110_009.sbx',...
%          '/home/elyall/Documents/Data/1955/150227/1955_110_011.sbx',...
%          '/home/elyall/Documents/Data/1955/150227/1955_110_019.sbx',...
%          '/home/elyall/Documents/Data/1955/150227/1955_110_101.sbx',...
%          '/home/elyall/Documents/Data/1955/150227/1955_110_109.sbx',...
%          '/home/elyall/Documents/Data/1955/150227/1955_110_111.sbx',...
%          '/home/elyall/Documents/Data/1955/150227/1955_110_119.sbx'};

Files = {...
    '/home/elyall/Documents/Data/1956/150304/1956_180_002.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_009.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_008.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_011.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_019.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_101.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_102.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_109.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_111.sbx',...
    '/home/elyall/Documents/Data/1956/150304/1956_180_119.sbx'};


SaveFiles = Files;
for findex = 1:numel(SaveFiles)
    SaveFiles{findex} = [SaveFiles{findex}(1:end-4),'_2.sbx'];
end

numFramesPerLoad = 760;


%% Fix data
for findex = 1:numel(Files)
    if ~exist(SaveFiles{findex}, 'file')
        
        Config = load2PConfig(Files{findex});
        totalFrames = Config.Frames;
        
        for bindex = 1:numFramesPerLoad:totalFrames % load frames in batches
            lastframe = min(bindex+numFramesPerLoad-1, totalFrames);
            currentFrames = bindex:lastframe;
            
            % Load frames
            [Images,~,Config] = load2P(Files{findex}, 'Type', 'Direct', 'Frames', currentFrames);
            
            % Compute mean
            FM = mean(reshape(Images(:,:,1,1,:), Config.size(1)*Config.size(2), numel(currentFrames)));
            FM = FM > 30000;
            
            % Fix subtraction
            Images(:,:,:,:,FM) = intmax('uint16') - Images(:,:,:,:,FM);
            
            % Save
            if bindex == 1
                load([Files{findex}(1:end-3),'mat'], 'info');
                save2P(SaveFiles{findex}, Images, 'Append', false, 'Header', info);
            else
                save2P(SaveFiles{findex}, Images, 'Append', true);
            end
            
        end
    end
end


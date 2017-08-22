
% % Load images
% FrameIndex = 1:7000; %1:2000;
% [imgs, loadObj, Config] = load2P('/media/elyall/Data/7120/170731/7120_250_003.sbx','Frames',FrameIndex,'Depth',2,'IndexType','relative','verbose');

% % Motion correct
% load('/media/elyall/Data/7120/170731/7120_250_003_depth2.align','MCdata','-mat');
% imgs = applyMotionCorrection(imgs,MCdata,loadObj);
% figure; imagesc(mean(imgs,5));

% % Crop bad regions
% Bad = [20,0,10,0]; % top,bottom,left,right
% imgs(end-Bad(2)+1:end,:,:,:,:) = [];
% imgs(1:Bad(1),:,:,:,:) = [];
% imgs(:,end-Bad(4)+1:end,:,:,:) = [];
% imgs(:,1:Bad(3),:,:,:) = [];

% % Save to tiff
% directory = '/media/elyall/Data/7120/170731/depth2/';
% str = 'image';
% filename = cell(size(imgs,5),1);
% for index = 1:size(imgs,5)
%     filename{index} = strcat(directory, str, num2str(FrameIndex(index)-1,'%05d' ), '.tiff');
%     save2P(filename{index},imgs(:,:,1,1,index));
% end

generateH5File('/home/elyall/Documents/Code/MATLAB/OutsideCode/HNCcorr/H5files/image2', '/home/elyall/Data/7120/170731/depth2', [size(imgs,1),size(imgs,2)], size(imgs,5), 31, 10);
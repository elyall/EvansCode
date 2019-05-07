function Savio_compAvgResponse(N)

addpath(genpath('/global/home/users/elyall/Code/Matlab/'))
parpool('local', 20);


% User defined inputs
str = '_catch';
secsBefore = 1;
secsAfter = 3;
% SaveDir = '/media/elyall/Data3/dFoFnew/';
SaveDir = '/global/scratch/elyall/dFoFnew/';

% set file to analyze
File = {};
% L2/3
for d = 1:4
    File = [File,{'/global/scratch/elyall/7142/170710/7142_220_002'}];
    File = [File,{'/global/scratch/elyall/6994/170725/6994_210_000'}];
    File = [File,{'/global/scratch/elyall/7120/170731/7120_250_003'}];
    File = [File,{'/global/scratch/elyall/7197/170807/7197_160_001'}];
end % 16
% L4
File = [File,{'/global/scratch/elyall/7734/180112/7734_338_000'}]; % 17
File = [File,{'/global/scratch/elyall/7734/180112/7734_308_001'}];
File = [File,{'/global/scratch/elyall/7736/180117/7736_300_000'}];
File = [File,{'/global/scratch/elyall/7736/180117/7736_265_001'}];
File = [File,{'/global/scratch/elyall/7737/180118/7737_291_000'}];
File = [File,{'/global/scratch/elyall/7737/180118/7737_326_001'}]; % 22
% L2/3 anesthetized
File = [File,{'/global/scratch/elyall/9445/181015/9445_180_005'}]; % 23
for d = 1:4
    File = [File,{'/global/scratch/elyall/9019/181021/9019_165_000'}];
    File = [File,{'/global/scratch/elyall/9025/181021/9025_180_002'}];
end % 31

% for d = 1:4
%     File = [File,{['/media/elyall/Data/7142/170710/7142_220_002_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data/6994/170725/6994_210_000_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data/7120/170731/7120_250_003_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data/7197/170807/7197_160_001_depth',num2str(d)]}];
% end
% File = [File,{'/media/elyall/Data/7734/180112/7734_338_000'}]; % 1st depth
% File = [File,{'/media/elyall/Data/7734/180112/7734_308_001'}]; % 2nd depth
% File = [File,{'/media/elyall/Data/7736/180117/7736_300_000'}]; % 1st depth
% File = [File,{'/media/elyall/Data/7736/180117/7736_265_001'}]; % 2nd depth
% File = [File,{'/media/elyall/Data/7737/180118/7737_291_000'}]; % 1st depth
% File = [File,{'/media/elyall/Data/7737/180118/7737_326_001'}]; % 2nd depth
% 
% File = [File,{'/media/elyall/Data2/9445/181015/9445_180_005'}];
% for d = 1:4
%     File = [File,{['/media/elyall/Data2/9019/181021/9019_165_000_depth',num2str(d)]}];
%     File = [File,{['/media/elyall/Data2/9025/181021/9025_180_002_depth',num2str(d)]}];
% end

% Set save name
[~,fn,~] = fileparts(File{N});
SaveFile = fullfile(SaveDir,[fn,str,'.avg']);

% Deal with multiple depths vs 1 depth
Depth = [repelem((1:4)',4);ones(7,1);repelem((1:4)',2)];
if N<17 || N>23
    AlignFile = sprintf([File{N},'_depth%d.align'],Depth);
else
    AlignFile = [File{N},'.align'];
end

% Compute average response per stimulus
load([File{N},'.exp'],'TrialIndex','-mat'); % load running index
computeAverageStimResponse([File{N},'.sbx'], [File{N},'.exp'], AlignFile,...
    'Depth',        Depth(N),...
    'TrialIndex',   TrialIndex,...
    'timeBefore',   secsBefore,...
    'timeAfter',    secsAfter,...
    'Save', 'SaveFile', SaveFile);



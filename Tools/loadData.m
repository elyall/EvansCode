function [data,Vars] = loadData(Files, Vars)
% Files is string or cell array of strings of filenames
% Vars is string or cell array of strings of variable names

if ischar(Files)
    Files = {Files};
end
numFiles = numel(Files);

if ~exist('Vars','var') || isempty(Vars)
    Vars = whos(matfile(Files{1}));
    Vars = {Vars(:).name};
elseif ischar(Vars)
    Vars = {Vars};
end
numVars = numel(Vars);

data = cell(numVars,numFiles);
str = sprintf('''%s'',',Vars{:});
for f = 1:numFiles
    eval(sprintf('temp=load(''%s'',%s''-mat'');',Files{f},str));
    data(:,f) = struct2cell(temp);
end

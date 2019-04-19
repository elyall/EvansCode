function saveEPS(fn,hF)
% SAVEEPS saves figure to .eps file
%   fn (string) - filename to save figure to (default: 'temp.eps')
%   hF (figure handle) - handle of figure to save (default: current figure)'

if ~exist('fn','var') || isempty(fn)
    fn = 'temp.eps';
end
if ~exist('hF','var') || isempty(hF)
    hF = gcf;
end

print(hF,fn,'-painters','-depsc','-tiff','-r300');

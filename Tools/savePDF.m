function savePDF(fn,hF)
% SAVEPDF saves figure to .pdf file
%   fn (string) - filename to save figure to (default: 'temp.pdf')
%   hF (figure handle) - handle of figure to save (default: current figure)

if ~exist('fn','var') || isempty(fn)
    fn = 'temp.pdf';
end
if ~exist('hF','var') || isempty(hF)
    hF = gcf;
end

print(hF,fn,'-painters','-dpdf','-r300');
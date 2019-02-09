function saveEPS(fn,hF)

if ~exist('fn','var') || isempty(fn)
    fn = 'temp.eps';
end
if ~exist('hF','var') || isempty(hF)
    hF = gcf;
end

print(hF,fn,'-depsc','-tiff','-r300','-painters');

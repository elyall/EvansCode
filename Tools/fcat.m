function out = fcat(in,dim)

if dim == 1
    mdim = 2;
elseif dim == 2
    mdim = 1;
end

sz = cellfun(@(x) size(x), in, 'UniformOutput', false);
sz = cat(1,sz{:});
out_size = max(sz(:,mdim));
toadd = out_size-sz(:,mdim);
for n = 1:numel(in)
    if mdim == 1
        in{n} = cat(mdim,in{n},nan(toadd(n),sz(n,2)));
    elseif mdim == 2
        in{n} = cat(mdim,in{n},nan(sz(n,1),toadd(n)));
    end
end
out = cat(dim,in{:});
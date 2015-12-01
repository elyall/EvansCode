function Filename=save2P(Images,Filename,Header,prompt,Channel)


%% Initialization
if ~exist('Images','var') || isempty(Images)
    if ~exist('Channel','var') % only matters if 'imgs' var is left out or blank
        Channel=1;
    end
    [Images,m]=load2P([],Channel);
    if isempty(Images)
        return
    end
end

if ~exist('Filenames','var') || isempty(Filename) % No filename input
    if exist('m','var')
        ids=uniquesites(m);
        Filename=strcat(ids{1},'.tif');
    else
        Filename='*.tif';
    end
    prompt=true; % not a full filename
elseif isdir(Filename) % Directory input
    Filename=strcat(Filename,'*.tif');
    prompt=true; % not a full filename
end

if ~exist('Header','var') || isempty(Header)
    Header='';
end

if ~exist('prompt','var')
    prompt=false;
end

%% Code
if prompt
    [Filename,p]=uiputfile({'*.tif'},'Save As',Filename); %only works in GUI environment (non-command line call)
    if isempty(Filename)
        return
    end
    Filename=fullfile(p,Filename);
end

Images=uint16(Images);
n=size(Images,3);
t=Tiff(Filename,'w');
wb=waitbar(0,sprintf('Saving: %s...',Filename));
for i=1:n
    if i~=1
        t.writeDirectory();
    end
    t.setTag('ImageLength',size(Images,1));
    t.setTag('ImageWidth', size(Images,2));
    t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
    t.setTag('BitsPerSample', 16);
    t.setTag('SamplesPerPixel', 1);
    t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    t.setTag('Software', 'MATLAB');
    if exist('header','var')
        t.setTag('ImageDescription',Header);
    end
    t.write(Images(:,:,i));
    waitbar(i/n,wb);
end
t.close();
close(wb);


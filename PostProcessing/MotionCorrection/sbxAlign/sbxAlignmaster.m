function sbxAlignmaster(fname,Depth,rect,Frames)

% rect - [y_top, y_bottom, x_left, x_right]

    global info
    
    if ~exist('rect','var')
        rect = [];
    end
    if isequal(rect,true)
        [~, ~, rect] = crop(sbxreadpacked(fname,0,1), rect);
        rect = round([rect(3),rect(3)+rect(4),rect(1),rect(1)+rect(2)]);
    end

    fprintf('Starting sbxAlignmaster: %s\n',[fname,'.sbx']);

    %% Edit to load and analyze single depth
    if ~exist('Depth','var') || isempty(Depth)
        Depth = 1;
    end
    Config = load2PConfig([fname,'.sbx']);
    numDepths = Config.Depth;
    if numDepths>1
        str = sprintf('_depth%d',Depth);
    else
        str = '';
    end
    
    if ~exist('Frames','var') || isempty(Frames)
        Frames = idDepth([fname,'.sbx'],[],'Depth',Depth)';
    else
        Frames = idDepth([fname,'.sbx'],[],'Depth',Depth,'Frames',Frames,'IndexType','absolute')';
    end
    Frames(Frames==1) = []; % throw out very first frame as it's wrong (blank row in T added in at end)
    numFrames = numel(Frames);
    Frames = [Frames;1:numFrames];

    %%

    z = sbxreadpacked(fname,0,1);

    szz = size(z);

    

    %Computing first order stats

    fprintf('Getting first-order stats\n');

    

    ms = 0;

    vs = 0;

    

    X = [ones(numFrames,1),linspace(-1,1,numFrames)'];

    X = bsxfun(@times,X,1./sqrt(sum(X.^2)));

    

    parfor jj = 1:numFrames

        z = double(sbxreadpacked(fname,Frames(1,jj)-1,1));


        ms = ms + z(:)*X(jj,:);

        vs = vs + z(:).^2;

    end

    

    s = reshape(sqrt(1/numFrames*(vs - sum(ms.^2,2))),szz);
    s = real(s); % occurs when constant 0 is acquired at some pixel on the frame; imaginary number caused by sqrt of negative number above (added by Evan)
    
    thestd = medfilt2(s,[31,31],'symmetric'); % 

    

    gl = X(:,2);

    l  = reshape(ms(:,2),szz);



    %%

    

    fprintf('Alignment first pass\n');

    

    [m,~,T] = sbxAlignpar(fname,thestd,gl,l,Frames,numDepths,rect); %Takes about 2.5 minutes


    if isempty(rect)
        rgx = (1:size(m,2))+45;
        rgy = 32 + (1:size(m,1));
        mask = [];
    else
        rgy = rect(1):rect(2);
        rgx = rect(3):rect(4);
        mask = false(szz);
        mask(rgy,rgx) = true;
    end

    T0 = T;

    

    for nn = 1:10

        fprintf('Refining alignment... pass %d\n',nn);

        [m,~,T] = sbxAligniterative(fname,m,rgy,rgx,thestd(rgy,rgx),gl,l,Frames,rect);

        dT = sqrt(mean(sum((T0-T).^2,2)));

        T0 = T;

        if dT < .25

            break;

        end

        fprintf('delta: %.3f\n',dT);

    end



    fprintf('Getting aligned first-order stats\n');

    

    ms = 0;

    vs = 0;

    

    m2 = 0;

    v2 = 0;

    

    X = [ones(numFrames,1),linspace(-1,1,numFrames)'];

    X = bsxfun(@times,X,1./sqrt(sum(X.^2)));

    

    g = exp(-(-5:5).^2/2/1.6^2);

    tic;

    parfor jj = 1:numFrames

        z = single(sbxreadpacked(fname,Frames(1,jj)-1,1));
        if ~isempty(mask) %Evan
            z(~mask) = 0;
        end
        
        z = z./thestd;
        z(isinf(z)) = 0; %Evan

        z = circshift(z,T(jj,:));

        z = double(z);



        ms = ms + z(:)*X(jj,:);

        vs = vs + z(:).^2;

        

        z = conv2(g,g,z,'same');

        m2 = m2 + z(:)*X(jj,:);

        v2 = v2 + z(:).^2;

    end

    toc;

    

    ss = sqrt(1/numFrames*(vs - sum(ms.^2,2)));

    m = reshape(ms(:,1)*X(1,1),size(l));

    v = reshape(ss.^2,size(l));

    
    ms(isnan(ms)) = 0;%Evan
    m2(isnan(m2)) = 0;%Evan
    v2(isnan(v2)) = 0;%Evan
    

    [Q,~] = qr(ms,0);

    [Q2,~] = qr(m2,0);

    s2 = reshape(sqrt(1/numFrames*(v2 - sum(m2.^2,2))),size(l));

    m2 = reshape(m2(:,1)*X(1,1),size(l));


    

    %In two passes, compute m, stdx, stdy

    %{

    m0 = 0;

    v0 = 0;

    

    m1 = 0;

    v1 = 0;

    

    g = exp(-(-5:5).^2/2/1.6^2);

    

    fprintf('Computing simple stats... pass %d\n',1);

    parfor jj = 1:info.max_idx

        z = double(sbxreadpacked(fname,jj-1,1));

        z = z - gl(jj)*l;

        z = circshift(z,T(jj,:));

        z = z./thestd;

        

        m0 = m0 + z;

        v0 = v0 + z.^2;



        z = conv2(g,g,z,'same');

        

        m1 = m1 + z;

        v1 = v1 + z.^2;

    end

    

    m = m0/info.max_idx;

    v = v0/info.max_idx;



    m2 = m1/info.max_idx;

    v2 = v1/info.max_idx;



    s2 = sqrt(v2 - m2.^2);

    %}

    
    s2(s2==0) = inf; %Evan
    
    k = 0;

    fprintf('Computing simple stats... pass %d\n',2);

    parfor jj = 1:numFrames % parfor

        z = double(sbxreadpacked(fname,Frames(1,jj)-1,1));
        if ~isempty(mask) %Evan
            z(~mask) = 0;
        end
        
        z = z - gl(jj)*l;

        z = circshift(z,T(jj,:));

        z = z./thestd;
        z(isinf(z)) = 0; %Evan
        

        z = conv2(g,g,z,'same');

        z = reshape(z(:) - (Q2*(Q2'*z(:))),size(z));

        k  =  k + (z./s2).^4;

    end

    s2(s2==inf) = 0; %Evan

    sm = s2./m2;
    sm(isnan(sm)) = 0; %Evan
    sm(isinf(sm)) = 0; %Evan
    
    k = k/numFrames - 3;

    

    if Depth==1
        T = cat(1,zeros(1,2),T); % add back values for first frame that was ignored
        Frames = [[1;0],Frames];
    end
    
    try

        save([fname,str,'.align'],'m','v','thestd','sm','k','T','Q','Frames','-append');

    catch

        save([fname,str,'.align'],'m','v','thestd','sm','k','T','Q');

    end

    

    tic;

    sbxComputeci(fname,Depth,Frames(1,:)); %Takes about 10 minutes, eats up a ton of RAM

    toc;

end
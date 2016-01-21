function sbxalignmaster(ImageFile)


Config = load2PConfig(ImageFile);
szz = Config.Size(1:2);



%Computing first order stats

fprintf('Getting first-order stats\n');



ms = 0;

vs = 0;



X = [ones(Config.Frames,1),linspace(-1,1,Config.Frames)'];

X = bsxfun(@times,X,1./sqrt(sum(X.^2)));



parfor jj = 1:Config.Frames
    
    z = load2P(ImageFile,'Frames',jj,'Double');
    
    
    ms = ms + z(:)*X(jj,:);
    
    vs = vs + z(:).^2;
    
end



s = reshape(sqrt(1/Config.Frames*(vs - sum(ms.^2,2))),szz);
try
    thestd = medfilt2(s,[31,31],'symmetric');
catch
    thestd = medfilt2(real(s),[31,31],'symmetric');
end



gl = X(:,2);

l  = reshape(ms(:,2),szz);



%%



fprintf('Alignment first pass\n');



[m,~,T] = sbxalignpar(ImageFile,thestd,gl,l); %Takes about 2.5 minutes



rgx = (1:size(m,2))+45;

rgy = 32 + (1:size(m,1));

T0 = T;



for nn = 1:10
    
    fprintf('Refining alignment... pass %d\n',nn);
    
    [m,~,T] = sbxaligniterative(ImageFile,m,rgy,rgx,thestd(rgy,rgx),gl,l);
    
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



X = [ones(Config.Frames,1),linspace(-1,1,Config.Frames)'];

X = bsxfun(@times,X,1./sqrt(sum(X.^2)));



g = exp(-(-5:5).^2/2/1.6^2);

tic;

parfor jj = 1:Config.Frames
    
    z = single(load2P(ImageFile,'Frames',jj));
    
    z = z./thestd;
    
    z = circshift(z,T(jj,:));
    
    z = double(z);
    
    
    
    ms = ms + z(:)*X(jj,:);
    
    vs = vs + z(:).^2;
    
    
    
    z = conv2(g,g,z,'same');
    
    m2 = m2 + z(:)*X(jj,:);
    
    v2 = v2 + z(:).^2;
    
end

toc;



ss = sqrt(1/Config.Frames*(vs - sum(ms.^2,2)));

m = reshape(ms(:,1)*X(1,1),size(l));

v = reshape(ss.^2,size(l));





[Q,~] = qr(ms,0);

[Q2,~] = qr(m2,0);

s2 = reshape(sqrt(1/Config.Frames*(v2 - sum(m2.^2,2))),size(l));

m2 = reshape(m2(:,1)*X(1,1),size(l));






k = 0;

fprintf('Computing simple stats... pass %d\n',2);

parfor jj = 1:Config.Frames
    
    z = load2P(ImageFile,'Frames',jj,'Double');
    
    z = z - gl(jj)*l;
    
    z = circshift(z,T(jj,:));
    
    z = z./thestd;
    
    
    
    z = conv2(g,g,z,'same');
    
    z = reshape(z(:) - (Q2*(Q2'*z(:))),size(z));
    
    k  =  k + (z./s2).^4;
    
end



sm = s2./m2;

k = k/Config.Frames - 3;


[p,f,~] = fileparts(ImageFile);
fname = fullfile(p,f);

try
    
    save([fname '.align'],'m','v','thestd','sm','k','T','Q','-mat','-append');
    
catch
    
    save([fname '.align'],'m','v','thestd','sm','k','T','Q','-mat','-v7.3');
    
end



tic;

sbxcomputeci(ImageFile); %Takes about 10 minutes, eats up a ton of RAM

toc;

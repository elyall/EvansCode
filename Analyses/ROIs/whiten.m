function [whData, invT, whT, Mean] = whiten(Data)


%% Center data
Mean = mean(Data);                  % calculate mean
Data = bsxfun(@minus, Data, Mean);  % subtract off mean

%% Perform Singular Value Decomposition (SVD)

% Compute SVD
A = Data'*Data;
[V,D,~] = svd(A);

% Ensure all elements on diagonal are nonzero
ind = find(diag(D)==0); % find zero elements on diagonal
D(ind,ind) = 0.00001;   % set all zero elements on diagonal to be nonzero

%% Compute whitening transform (and inverse)
whT = sqrt(size(Data,1)-1)*V*sqrtm(inv(D))*V'; % determine whitening transform
invT = pinv(whT);                      % determine inverse of transform

%% Whiten data
whData = Data*whT;


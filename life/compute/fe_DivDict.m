function [fe] = fe_DivDict(varargin)
fe = varargin{1};
n = varargin{2};
lambda = varargin{3};

% feFitPhi() function that given a fix vector of weights w optimize the tensor Phi
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

[nAtoms] = feGet(fe,'natoms');
[nTheta] = size(fe.life.M.DictSig,1);
[nVoxels] = size(fe.life.M.Phi,2);

w = fe.life.fit.weights;

ind_vox = fe.life.M.ind_vox{n};
nVoxels_red = length(ind_vox);

% restrict signal to voxels
dSig = reshape(feGet(fe,'dsigdemeaned'),[nTheta,nVoxels]);
dSig = dSig(:,ind_vox);

% restrict Phi to voxels
Phi = fe.life.M.Phi;
Phi = Phi(:,ind_vox,:);

B = ttv(Phi,w,3);
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels_red);

% Free space for extra dictionary
nDict = size(fe.life.M.Dictionaries,2);
for m=nDict+1:-1:n+2
    fe.life.M.Dictionaries{m} = fe.life.M.Dictionaries{m-1};
    fe.life.M.ind_vox{m} = fe.life.M.ind_vox{m-1};
end

%% Divide Dictionary n into (n, n+1)
Dict = fe.life.M.Dictionaries{n};
e = sum((dSig - Dict*B).^2,1)./sum(dSig.^2,1);
epsilon = median(e);
ind_voxA = e < epsilon; % indices to voxels with low error

%fe.life.M.Dictionaries{n} = Dict;
fe.life.M.ind_vox{n} = ind_vox(ind_voxA);

ind_voxB = e >= epsilon; %indices to voxels with high error need to be fitted with a new dictionary (next level)
fe.life.M.Dictionaries{n+1} = Dict;
fe.life.M.ind_vox{n+1} = ind_vox(ind_voxB);

end


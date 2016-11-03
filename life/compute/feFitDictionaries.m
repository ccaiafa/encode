function [fe] = feFitDictionaries(varargin)
fe = varargin{1};
lambda = varargin{2};


% feFitPhi() function that given a fix vector of weights w optimize the tensor Phi
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

%load variables.mat

%[nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
[nAtoms] = feGet(fe,'natoms');
[nTheta] = size(fe.life.M.DictSig,1);
[nVoxels] = feGet(fe,'nvoxels');

w = fe.life.fit.weights;

dSig = reshape(feGet(fe,'dsigdemeaned'),[nTheta,nVoxels]);

B = ttv(fe.life.M.Phi,w,3);
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);

nDict = size(fe.life.M.Dictionaries,2);
for n=1:nDict
disp(['level=',num2str(n)])    
[fe.life.M.Dictionaries{n}, B(:,fe.life.M.ind_vox{n})] = update_DandB(dSig(:, fe.life.M.ind_vox{n}), fe.life.M.Dictionaries{n}, B(:,fe.life.M.ind_vox{n}), lambda);    
end

%% Compute Phi compatible with B(a,v)
[sub, ~] = find(fe.life.M.Phi);
ind = sub2ind(size(B), sub(:,1), sub(:,2));
b = B(:);

A = sptensor(sub, ones(size(sub,1),1), size(fe.life.M.Phi));
A = ttv(A,w.^2,3);
[subA, val] = find(A);
A = sparse(subA(:,1),subA(:,2),val,nAtoms,nVoxels);
a = A(:);
div = a(ind);
[ind1,~,val] = find(div);

newval = full(w(sub(:,3)).*b(ind));
newval(ind1) = newval(ind1)./val;

fe.life.M.Phi = sptensor(sub, newval, size(fe.life.M.Phi));

end

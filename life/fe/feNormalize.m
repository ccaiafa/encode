function [ fe ] = feNormalize( fe )
% Adjust weights according to the mean diffusion signal per fascicle
nTheta  = feGet(fe,'nbvecs');
nVoxels = feGet(fe,'nvoxels');
A = ttv(fe.life.M.Phi,ones(size(fe.life.M.Phi,1),1),1); % sum over atoms
A = ttv(A, ones(size(A,1),1),1); % sum over voxels
a = double(A);
a = full(a)/(nTheta*nVoxels);
fe.life.M.Phi = sptensor(fe.life.M.Phi.subs, fe.life.M.Phi.vals./a(fe.life.M.Phi.subs(:,3)),size(fe.life.M.Phi)); % Normalize tensor Phi

fe.life.fit.weights = fe.life.fit.weights.*a; % Normalize vector of weights

end


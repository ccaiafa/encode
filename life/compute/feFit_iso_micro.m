function [fit] = feFit_iso_micro(M, w, bval, bvecs, S0, dSig)
% 
%   
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F.  Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

%% The following part allows fitting isotropic weight (w0) ad microstructure (Qa)

% Compute signal contributed by fascicles (with mean)
%Sf = M_times_w_with_mean(M,w);
nTheta = size(M.Dict,1);
nVoxels = size(M.Phi,2);
Sf = M_times_w(M.Phi.subs(:,1),M.Phi.subs(:,2),M.Phi.subs(:,3),M.Phi.vals,M.Dict,w,nTheta,nVoxels);

Sf = reshape(Sf,[nTheta,nVoxels]);

% Extract unique bval
bval = bval(1);

% Compute sum of weights per voxel
nFibers = size(M.Phi,3);
vox_fib = [M.Phi.subs(:,2),M.Phi.subs(:,3)];
vox_fib = unique(vox_fib,'rows');

A = sparse(vox_fib(:,1),vox_fib(:,2),ones(size(vox_fib,1),1),nVoxels,nFibers);
sumwf = A*w;
sumwf(sumwf>1) = 1;

D0 = 3; %[Alexander, 2001] water free isotropic diffusion

% Compute isotropic weight (w0) and microstructure tensor (Qa) per voxel
w0 = zeros(nVoxels,1);
wa = zeros(nVoxels,1);
Qa = cell(nVoxels,1);
parfor v=1:nVoxels
    rangewf = linspace(1e-9,1-sumwf(v)-1e-9,20);
    [w0(v), Qa{v}] = fit_iso_micro(S0(v), bval, D0, bvecs', dSig(:,v), Sf(:,v), rangewf, sumwf(v));
    wa(v) = 1 - w0(v) - sumwf(v);
    disp(['voxel ', num2str(v),'/',num2str(nVoxels)])   
end

fit.w0                  = w0;
fit.Qa                  = Qa;
fit.wa                  = wa;
fit.b                   = bval;
fit.D0                  = D0;

end

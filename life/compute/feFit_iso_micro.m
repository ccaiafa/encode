function [fit] = feFit_iso_micro(M, w, bval, bvecs, S0, dSig)
% 
%   
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F.  Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

%% The following part allows fitting isotropic weight (w0) ad microstructure (Qa)

% Compute signal contributed by fascicles (with mean)
Sf = M_times_w_with_mean(M,w);
nTheta = size(M.Dict,1);
Nvoxels = size(M.Phi,2);
Sf = reshape(Sf,[nTheta,Nvoxels]);

% Extract unique bval
bval = bval(1);

% Compute sum of weights per voxel
nFibers = size(M.Phi,3);
A = sparse(M.Phi.subs(:,2),M.Phi.subs(:,3),ones(nnz(M.Phi),1),Nvoxels,nFibers);
sumwf = A*w;
sumwf(sumwf>1) = 1;

D0 = 3; %[Alexander, 2001] water free isotropic diffusion

% Compute isotropic weight (w0) and microstructure tensor (Qa) per voxel
w0 = zeros(Nvoxels,1);
wa = zeros(Nvoxels,1);
Qa = cell(Nvoxels,1);
parfor v=1:Nvoxels
    rangewf = linspace(1e-9,1-sumwf(v)-1e-9,100);
    [w0(v), Qa{v}] = fit_iso_micro(S0(v), bval, D0, bvecs', dSig(:,v), Sf(:,v), rangewf, sumwf(v));
    wa(v) = 1 - w0(v) - sumwf(v);
    disp(['voxel ', num2str(v),'/',num2str(Nvoxels)])   
end

fit.w0                  = w0;
fit.Qa                  = Qa;
fit.wa                  = wa;
fit.b                   = bval;
fit.D0                  = D0;

end

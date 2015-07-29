function [fit, w, R2] = feFitModel(M,dSig,dSig_demeaned,S0,fitMethod,Niter,preconditioner)
% 
% feFitModel() function in LiFE but restricted to the
% 
% BBNNLS algorithm and using the Factorization model.
% M is the factorization model composed by:
%
%   M.D_demean    Dictionary
%   
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com


% Below are the old comments which are not longer correct in general
% Fit the LiFE model.
%
% Finds the weights for each fiber to best predict the directional
% diffusion signal (dSig_demeaned)
%
%  fit = feFitModel(M,dSig_demeaned,S0,fitMethod)
%
% dSig_demeaned:  The diffusion weighted signal measured at each
%        voxel in each direction. These are extracted from 
%        the dwi data at some white-matter coordinates.
% M:     The LiFE difusion model matrix, constructed
%        by feConnectomeBuildModel.m
%
% fitMethod: 
%  - 'bbnnls' - DEFAULT and best, faster large-scale solver.
%
% See also: feCreate.m, feConnectomeBuildModel.m, feGet.m, feSet.m
%
% Example:
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%
% Notes about the LiFE model:
%
% The rows of the M matrix are nVoxels*nBvecs. We are going to predict the
% diffusion signal in each voxel for each direction.
%
% The columns of the M matrix are nFibers + nVoxels.  The diffusion signal
% for each voxel is predicted as the weighted sum of predictions from each
% fibers that passes through a voxel plus an isotropic (CSF) term.
%
% In addition to M, we typically return dSig_demeaned, which is the signal measured
% at each voxel in each direction.  These are extracted from the dwi data
% and knowledge of the roiCoords.

% fit the model, by selecting the proper toolbox.

mycomputer = computer();
release = version('-release');

[nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
[nTheta]  = size(M.D_demean,1);
[nAtoms] = size(M.D_demean,2); %feGet(fe,'natoms');
[Nvoxels] = size(M.Phi,2); %feGet(fe,'nvoxels');

if strcmp(preconditioner,'preconditioner')
    switch strcat(mycomputer,'_',release)
        case {'GLNXA64_2015a','MACI64_2014b','MACI64_2015a'}
            h = compute_diag_mex(M.Phi.subs(:,1), M.Phi.subs(:,3), M.Phi.vals, M.D_demean,nFibers);
            vals = M.Phi.vals./h(M.Phi.subs(:,3));
            M.Phi = sptensor(M.Phi.subs,vals,size(M.Phi));
        otherwise
            h = dot(M.MmatrixM,M.MmatrixM,1);
            M.MmatrixM = M.MmatrixM*diag(h.^(-1));
    end
end

switch fitMethod
   case {'bbnnls'}
    
    tic
    fprintf('\nLiFE: Computing least-square minimization with BBNNLS...\n')
    opt = solopt;
    opt.maxit = Niter;
    opt.use_tolo = 1; 
    opt.tolg = 1e-4;
    
    switch strcat(mycomputer,'_',release)
        case {'GLNXA64_2015a'}
        out_data = bbnnls_GLNXA64(M,dSig_demeaned,zeros(nFibers,1),opt);
        case {'MACI64_2014b','MACI64_2015a'}
        out_data = bbnnls_MACI64(M,dSig_demeaned,zeros(nFibers,1),opt);
        otherwise
        sprintf('WARNING: currently LiFE is optimized for an efficient usage of memory \n using the Sparse Tucker Decomposition aproach (Caiafa&Pestilli, 2015) \n ONLY for Linux (MatlabR2015a) and MacOS (Matlab R2014b/R2015a). \n If you have a different system or version you can still \n use the old version of LiFE (memory intensive). \n\n')
        sprintf('\n Starting using old version of LiFE...\n')
        out_data = bbnnls_OLD(M.MmatrixM,dSig_demeaned,zeros(nFibers,1),opt);
    end
    
    if strcmp(preconditioner,'preconditioner')
        out_data.x = out_data.x./h;
    end
    
    fprintf('BBNNLS status: %s\nReason: %s\n',out_data.status,out_data.termReason);
    w = out_data.x;
    fprintf(' ...fit process completed in %2.3fminutes\n',toc/60)
    % Save the state of the random generator so that the stochasit cfit can be recomputed.
    defaultStream = RandStream.getGlobalStream; %RandStream.getDefaultStream;
    fit.randState = defaultStream.State;   
    
    % Save out some results 
    fit.results.R2        = [];
    fit.results.nParams   = size(M,2);
    fit.results.nMeasures = size(M,1);
    R2=[];

   otherwise
     error('Cannot fit LiFE model using method: %s.\n',fitMethod);
end

%% Compute isotropic weight per voxel (w0)
w0 = ttv(M.Phi,w,3);
w0 = double(sptenmat(w0,1));
w0 = ones(1,Nvoxels) - sum(w0,1)./S0';
w0 = w0';

%% Fit isotropic term per voxel (A0)
y = ttv(M.Phi,w,3); % This is memory efficient version of above
% The following makes ans \times_1 Dic
if nnz(y)
    y = ttm(y,M.Dict,1);
else
    y = sptensor([size(M.Dict,1),size(y,2)]); % In case empty tensor
end
y = tenmat(y,1); % sptensor -> tenmat
y = double(y); % tenmat -> sparse matrix

y = dSig - y;
y = y./repmat((w0.*S0)',nTheta,1);
y(y<0) = 0;
A0 = -log(mean(y))';

% Save output structure.
fit.weights             = w;
fit.params.fitMethod    = fitMethod;
fit.w0                  = w0;
fit.A0                  = A0;

end

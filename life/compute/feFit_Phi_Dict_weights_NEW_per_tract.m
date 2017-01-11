function [ fe, error_out ] = feFit_Phi_Dict_weights_NEW_per_tract( fe, error_trhesholdSig, error_threshold_weights, Niter, Niter_loop, varargin)

lambda = 0; % Parameter that control l2 regularization (maybe will be eliminated in the future)

nTheta  = feGet(fe,'nbvecs');
nVoxels = feGet(fe,'nvoxels');
weights_old = fe.life.fit.weights;

% Initialize with only one dictionary for all voxels
fe.life.M.Dictionaries{1} = fe.life.M.DictSig;
fe.life.M.ind_vox{1} = 1:nVoxels;
nDict = 1;
% Fit Dictionaries
fe =  feFitDictionaries(fe,lambda);

% Normalize Phi and w
fe = feNormalize(fe);

%fe = fe_DivDict(fe, 1, lambda ); % Dvide first dictionary
%nDict = 2;
%fe =  feFitDictionaries(fe,lambda); % Fit new Dictionaries
pred_full = feGet(fe,'predfull');
meas_full = feGet(fe,'dsigmeasured');
norm_full = norm(meas_full,'fro');
error_full =  norm(meas_full - pred_full,'fro')/norm_full;

fe.life.M.splits = [];

i = 1;
error_out(i) = error_full;
error_full_old = error_full;
delta_error = Inf;
%delta_error = 0; %% NOT TO ENTER IN THE GLOBAL DICTIONARY LEARNING AND JUMP INTO THE DICT PER TRACT OPTIMIZATION
delta_weights = Inf;
iFit =1;
while (delta_weights > error_threshold_weights) && (delta_error > error_trhesholdSig) && (iFit <= Niter_loop)
    
    % fit weights
    fit = feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'nopreconditioner',fe.life.fit.weights);
    fe = feSet(fe,'fit',fit);
    
    % Fit Dictionaries
    fe =  feFitDictionaries(fe,lambda);
    
    % Normalize Phi and w
    fe = feNormalize(fe);
    
    pred_full = feGet(fe,'predfull');
    error_full =  norm(meas_full - pred_full,'fro')/norm_full;
    delta_error = abs(error_full_old - error_full)/error_full_old;
    delta_weights = norm(weights_old - fe.life.fit.weights)/norm(weights_old);
    
    disp(['Number of Dictionaries=',num2str(nDict),' Fitting=',num2str(iFit)]);
    disp([' Delta weights=',num2str(100*delta_weights),'%', '  Delta Fit= ', num2str(100*delta_error), '%', ' error full=',num2str(100*error_full),'%']);
    
    weights_old = fe.life.fit.weights;
    error_full_old = error_full;
    iFit = iFit + 1;
    
    error_out(i) = error_full;
    i = i + 1;
end

if ~isempty('varargin')
    dataOutputPath = varargin{1};
    name = varargin{2};
    save(fullfile(dataOutputPath,sprintf('fe_struct_%s_per_tract.mat',name)), 'fe','error_out','-v7.3')
end



%% start fitting dictionaries per tract

fe.life.M.splits = i-1;
nTract = size(fe.life.M.tracts,2);
for n=1:nTract
    fe.life.M.tracts{n}.DictSig = fe.life.M.Dictionaries{1};
end

delta_error = Inf;
delta_weights = Inf;
iFit =1;
while (delta_weights > error_threshold_weights) && (delta_error > error_trhesholdSig) && (iFit <= Niter_loop)
    % fit weights
    fit = feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'nopreconditioner',fe.life.fit.weights);
    fe = feSet(fe,'fit',fit);
    
    % Fit Dictionaries
    fe =  feFitDictionaries_per_tract(fe,lambda,error_trhesholdSig,Niter_loop);
        
    % Normalize Phi and w
    fe = feNormalize(fe);
    
    pred_full = feGet(fe,'predfull');
    error_full =  norm(meas_full - pred_full,'fro')/norm_full;
    delta_error = abs(error_full_old - error_full)/error_full_old;
    delta_weights = norm(weights_old - fe.life.fit.weights)/norm(weights_old);
    
    disp(['Fitting Dicts per tract, Iter = ',num2str(iFit)]);
    disp([' Delta weights=',num2str(100*delta_weights),'%', '  Delta Fit= ', num2str(100*delta_error), '%', ' error full=',num2str(100*error_full),'%']);
    
    weights_old = fe.life.fit.weights;
    error_full_old = error_full;
    iFit = iFit + 1;
    
    error_out(i) = error_full;
    i = i + 1;
    
    if ~isempty('varargin')
        dataOutputPath = varargin{1};
        name = varargin{2};
        save(fullfile(dataOutputPath,sprintf('fe_struct_%s_per_tract.mat',name)), 'fe','error_out','-v7.3')
    end
    
end





end


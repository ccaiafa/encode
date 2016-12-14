function [ fe ] = feFit_Phi_Dict_weights( fe, error_trhesholdSig, error_threshold_weights, nDictMax, Niter, Niter_loop,dataOutputPath)

lambda = 0; % Parameter that control l2 regularization (maybe will be eliminated in the future)

nTheta  = feGet(fe,'nbvecs');
nVoxels = feGet(fe,'nvoxels');
weights_old = fe.life.fit.weights;

error_full_old = Inf;
% Initialize with only one dictionary for all voxels
fe.life.M.Dictionaries{1} = fe.life.M.DictSig;
fe.life.M.ind_vox{1} = 1:nVoxels;

fe = fe_DivDict(fe, 1, lambda ); % Dvide first dictionary
nDict = 2;
fe =  feFitDictionaries(fe,lambda); % Fit new Dictionaries
pred_full = feGet(fe,'predfull');
meas_full = feGet(fe,'dsigmeasured');
norm_full = norm(meas_full,'fro');
error_full =  norm(meas_full - pred_full,'fro')/norm_full;

while (error_full > error_trhesholdSig)&& (nDict < nDictMax)
    %%     
    iFit =1;
    delta_weights = Inf;
    while (delta_weights > error_threshold_weights) && (iFit <= Niter_loop)

        % fit weights
        fit = feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'nopreconditioner',fe.life.fit.weights);
        fe = feSet(fe,'fit',fit);
        
        % Fit Dictionaries
        fe =  feFitDictionaries(fe,lambda);
        
        pred_full = feGet(fe,'predfull');
        error_full =  norm(meas_full - pred_full,'fro')/norm_full;
        delta_error = abs(error_full_old - error_full)/error_full_old;
        delta_weights = norm(weights_old - fe.life.fit.weights)/norm(weights_old);
        
        disp(['Number of Dictionaries=',num2str(nDict),' Fitting=',num2str(iFit),' Delta weights=',num2str(100*delta_weights),'%',' error full=',num2str(100*error_full),'%']);
        
        weights_old = fe.life.fit.weights;
        error_full_old = error_full;
        iFit = iFit + 1;
    end
    
    % Try to divide last dictionary. To be modified to test all possible
    % Dictionary split and chose the one that fives better fit.
    
    for n=nDict:-1:1
        fe = fe_DivDict(fe, n, lambda ); % Dvide first dictionary
        nDict = nDict + 1;
    end
    
    fe =  feFitDictionaries(fe,lambda); % Fit new Dictionaries
    
    save(fullfile(dataOutputPath,sprintf('fe_structure_Dict_Learn_Dicts%s.mat',num2str(nDict))), 'fe','-v7.3')
    

end
end


function [ fe, error_out ] = feFit_Phi_Dict_weights_NEW( fe, error_trhesholdSig, error_threshold_weights, nDictMax, Niter, Niter_loop, Nitermax, varargin)

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
while ((delta_error > error_trhesholdSig) || (nDict <= nDictMax)) &&  (i < Nitermax) % while NOT CONVERGED OR Max number of Dict not reached
    %%     
    iFit =1;
    delta_weights = Inf;
    while ((delta_weights > error_threshold_weights) || (delta_error > error_trhesholdSig)) && (iFit <= Niter_loop) % While NOT Converged AND Max Iteration not reached

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
        
        if ~isempty('varargin')
            dataOutputPath = varargin{1};
            name = varargin{2};
            save(fullfile(dataOutputPath,sprintf('error_out_%s_nDicts%s.mat',name,num2str(nDict))),'error_out','-v7.3')
        end
    end
    
    % Save fe structure before doubling nDict
    if ~isempty('varargin')
        dataOutputPath = varargin{1};
        name = varargin{2};
        save(fullfile(dataOutputPath,sprintf('fe_struct_%s_nDicts%s.mat',name,num2str(nDict))), 'fe','error_out','-v7.3')
    end

    % Try to divide last dictionary. To be modified to test all possible
    % Dictionary split and chose the one that fives better fit.
    
    if nDict < nDictMax
        fe.life.M.splits = [fe.life.M.splits, i-1];
        for n=nDict:-1:1
            fe = fe_DivDict(fe, n, lambda ); 
            nDict = nDict + 1;
        end
        
        fe =  feFitDictionaries(fe,lambda); % Fit new Dictionaries
        
        % Normalize Phi and w
        fe = feNormalize(fe);       
        delta_error = Inf;      
    end
    

end
end


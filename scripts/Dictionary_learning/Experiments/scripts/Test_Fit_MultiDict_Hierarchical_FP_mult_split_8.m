dataRootPath = '/N/dc2/projects/lifebid/code/ccaiafa/demo_datasets/STN/';
subject = 'FP';
conn = 'NUM01'; % 
param = 'lmax10'; % {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8', ''}
alg = 'SD_PROB'; % {'SD_PROB', 'SD_STREAM','tensor'}

vista_soft_path = '/N/dc2/projects/lifebid/code/vistasoft/';
%vista_soft_path = '/Users/CesarMB13/SOFT/vistasoft/';
addpath(genpath(vista_soft_path));

%feStructurePath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/';

encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
%encode_path = '/Users/CesarMB13/SOFT/encode/';
addpath(genpath(encode_path));

dataOutputPath = '/N/dc2/projects/lifebid/code/ccaiafa/encode/scripts/Dictionary_learning/Experiments/Results/';

% Generate fe_strucure
%% Build the file names for the diffusion data, the anatomical MRI.
dwiFile       = deblank(ls(fullfile(dataRootPath,strcat('sub-',subject),'dwi','run01_fliprot_aligned_trilin.nii.gz')));
dwiFileRepeat = deblank(ls(fullfile(dataRootPath,strcat('sub-',subject),'dwi','run02_fliprot_aligned_trilin.nii.gz')));
t1File        = deblank(fullfile(dataRootPath,strcat('sub-',subject),'anatomy',  't1.nii.gz'));

fgFileName    = deblank(ls(fullfile(dataRootPath,strcat('sub-',subject),'tractography', 'run01_fliprot_aligned_trilin_csd_lmax10_wm_SD_PROB-NUM01-500000.tck')));
feFileName    = 'STN_FP_test'; 

L = 90; % Discretization parameter
Niter = 50;
%fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFile,t1File,L,[1.8,0.8]); % We set dwiFileRepeat =  run 01
% fit = feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'nopreconditioner');
% fe = feSet(fe,'fit',fit);
% 
% %% Save reference fe-structure
% save('fe_structure_105115_L90.mat', 'fe','-v7.3');


%% Load fe structure with single Dictionary learned from data
% load('/N/dc2/projects/lifebid/code/ccaiafa/Development/Dictionary_learning/Experiments/results/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01_LAMBDA05.mat')
% pred_full = feGet(fe,'predfull');
% meas_full = feGet(fe,'dsigmeasured');
% error_full_old =  norm(meas_full - pred_full,'fro')/norm(meas_full,'fro');
% disp(['Iter=',num2str(0),'  error full=',num2str(100*error_full_old),'%']);

fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFile,t1File,L,[1.8,0.8]); % We set dwiFileRepeat =  run 01
% Initialize weights
fe.life.fit.weights = ones(feGet(fe,'nfibers'),1);
normw = sum(fe.life.fit.weights);
fe.life.fit.weights = fe.life.fit.weights/normw;

   
%% Here we start with the main loop that divide Dictionaries corresponding to different voxels

lambda = 0; % Parameter that control l2 regularization (maybe will be eliminated in the future)

nTheta  = feGet(fe,'nbvecs');
nVoxels = feGet(fe,'nvoxels');
weights_old = fe.life.fit.weights;

error_full_old = Inf;
% Initialize with only one dictionary for all voxels
fe.life.M.Dictionaries{1} = fe.life.M.DictSig;
fe.life.M.ind_vox{1} = 1:nVoxels;

fe_new = fe_DivDict(fe, 1, lambda ); % Dvide first dictionary
nDict = 2;
fe_new =  feFitDictionaries(fe_new,lambda); % Fit new Dictionaries
pred_full = feGet(fe,'predfull');
meas_full = feGet(fe,'dsigmeasured');
error_full =  norm(meas_full - pred_full,'fro')/norm(meas_full,'fro');
%if error_full < 0.95* error_full_old
    fe = fe_new;
    exit_cond = false;
%else
%    exit_cond = True;
%    disp(['not improvement dividing dictionary ', num2str(1)])
%end


while (~exit_cond)&& (nDict <= 8)
    %%     
    iFit =1;
    delta_weights = Inf;
    while (delta_weights > 0.0001) && (iFit <= 20)

        % fit weights
        fit = feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'nopreconditioner',fe.life.fit.weights);
        fe = feSet(fe,'fit',fit);
        
        % Fit Dictionaries
        fe =  feFitDictionaries(fe,lambda);
        
        pred_full = feGet(fe,'predfull');
        meas_full = feGet(fe,'dsigmeasured');
        error_full =  norm(meas_full - pred_full,'fro')/norm(meas_full,'fro');
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
%     pred_full = feGet(fe,'predfull');
%     meas_full = feGet(fe,'dsigmeasured');
%     error_full =  norm(meas_full - pred_full,'fro')/norm(meas_full,'fro');
%     %if error_full < 0.95* error_full_old % To add a check on size of voxel set (avoid small voxel sets (overfitting))
%         fe = fe_new;
%         exit_cond = false;
%         error_full_old = error_full;
%     else
%         exit_cond = True;
%         disp(['not improvement dividing dictionary ', num2str(iter+1)])
%     end
    
    
    save(fullfile(dataOutputPath,sprintf('fe_structure_MultiHier_mult_split_%s_%s_%s_%s_conn%s_%sDicts_32.mat',subject,'STC_run01',char(alg),char(param),conn,num2str(nDict))), 'fe','-v7.3')
    

end


rmpath(genpath(vista_soft_path));
rmpath(genpath(encode_path));






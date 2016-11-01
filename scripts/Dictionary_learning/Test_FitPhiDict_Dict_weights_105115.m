dataRootPath = '/N/dc2/projects/lifebid/2t1/HCP/';
subject = '105115';
conn = 'NUM01'; % 
param = 'lmax10'; % {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8', ''}
alg = 'SD_PROB'; % {'SD_PROB', 'SD_STREAM','tensor'}

vista_soft_path = '/N/dc2/projects/lifebid/code/vistasoft/';
addpath(genpath(vista_soft_path));

feStructurePath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/';

encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
addpath(genpath(encode_path));

% Generate fe_strucure
%% Build the file names for the diffusion data, the anatomical MRI.
dwiFile       = deblank(ls(fullfile(dataRootPath,subject,'diffusion_data','*b2000_aligned*.nii.gz')));
dwiFileRepeat = deblank(ls(fullfile(dataRootPath,subject,'diffusion_data','*b2000_aligned*.nii.gz')));
t1File        = deblank(fullfile(dataRootPath,subject,'anatomy',  'T1w_acpc_dc_restore_1p25.nii.gz'));

fgFileName    = deblank(ls(fullfile(dataRootPath,subject,'fibers_new', strcat('*b2000*',char(param),'*',char(alg),'*',conn,'*','500000.tck'))));
feFileName    = 'HCP105115test'; 

L = 90; % Discretization parameter
Niter = 500;
%fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFile,t1File,L,[1.8,0.8]); % We set dwiFileRepeat =  run 01
% fit = feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'nopreconditioner');
% fe = feSet(fe,'fit',fit);
% 
% %% Save reference fe-structure
% save('fe_structure_105115_L90.mat', 'fe','-v7.3');


%% Create a new fe structure
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFile,t1File,L,[1.8,0.8]); % We set dwiFileRepeat =  run 01
% 
% % simulate diffusion signal
% pred_full =  feGet(fe,'predfull');
% meas_full = feGet(fe,'dsigmeasured');
% error_full =  norm(meas_full - pred_full,'fro')/norm(meas_full,'fro');
% disp(['Error predicting full signal=',num2str(error_full)]);
% fe.life.diffusion_signal_img = pred_full';

% Initialize weights
fe.life.fit.weights = ones(feGet(fe,'nfibers'),1);
normw = sum(fe.life.fit.weights);
fe.life.fit.weights = fe.life.fit.weights/normw;

% Initialize Tensor model
%fe = feSet(fe,'model tensor',[1.5 0.2 0.2]);

nTheta  = feGet(fe,'nbvecs');
nVoxels = feGet(fe,'nvoxels');
weights_old = fe.life.fit.weights;
for iter=1:16
%% Fix weights optimize Phi & Dict
lambda = 0;
[fe ]= feInitFitDictionaries(fe,lambda);
% Adjust weights according to the mean diffusion signal per fascicle
% A = ttv(Phi,ones(size(Phi,1),1),1); % sum over atoms
% A = ttv(A, ones(size(A,1),1),1); % sum over voxels
% a = double(A);
% a = full(a)/(nTheta*nVoxels);
% Phi = sptensor(Phi.subs, Phi.vals./a(Phi.subs(:,3)),size(Phi)); % Normalize tensor Phi

% fe = feSet(fe,'Indication Tensor',Phi);
% fe.life.M.DictSig = Dict;

%fe.life.M.AdpDict{iter} = Dict; 
% fe.life.M.Qtensors{iter} = Qtensors; 
% fe.life.M.Atom_mean{iter} = Atom_mean;
% fe.life.M.Tensor_fit_error{iter} = Tensor_fit_error;
%fe.life.M.Vox_per_atom = Vox_per_atom;

% fe.life.fit.weights = fe.life.fit.weights.*a;
%fe.life.M.weights{iter} = fe.life.fit.weights;


pred_full = feGet(fe,'predfull'); %%%%%%%%% Cambiar!!!!!!
meas_full = feGet(fe,'dsigmeasured');
error_full =  norm(meas_full - pred_full,'fro')/norm(meas_full,'fro');

error_weights = norm(weights_old - fe.life.fit.weights)/norm(weights_old); 
weights_old = fe.life.fit.weights;
rel_error = feGet(fe,'relative error');
disp(['Iter=',num2str(iter),' Delta weights=',num2str(100*error_weights),'%',' rel error=',num2str(100*rel_error),'%',' error full=',num2str(100*error_full),'%']);
disp(['nnz w=',num2str(nnz(fe.life.fit.weights)),' nnz Phi=',num2str(nnz(fe.life.M.Phi))])

results{iter}.AdpDict = Dict;
results{iter}.weights = fe.life.fit.weights;
results{iter}.error_weights = error_weights;
results{iter}.rel_errors = rel_error;
results{iter}.error_full = error_full;

fit = feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',500,'nopreconditioner',fe.life.fit.weights);
fe = feSet(fe,'fit',fit);

save('fe_structure_Dict_weights_105115_50iter_eps_lambda0.mat','fe','-v7.3');
save('results_Dict_weights_105115_50iter_eps_lambda0.mat','results','-v7.3');
end


rmpath(genpath(dataRootPath));
rmpath(genpath(vista_soft_path));
rmpath(genpath(encode_path));






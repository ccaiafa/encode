%% Example of script for the Concussion model
% The Concussion model decomposes the observed diffusion signal in three
% components (Isotropic + Fibers + Microstructure):
% 1) The Isotropic signal due free motion of water molecules
% 2) The Anisotropic signal due to fascicle's (directional features)
% 3) The Anisotropic signal due to other directional features not related
% with the fascicles (microstructure)
%
% After the model is fit, this script generates separated nifti files for
% each of the components.
% Also, for each one of the 20 major tracts, we compute its prediction and
% generate individual nifti files with the diffusion signal associated only
% to the fibers in those corresponding tracts.

%% Set the Path for your output
dataOutputPath = '/N/dc2/projects/lifebid/code/ccaiafa/Development/Concussion_model/results/';

%% Set the proper path for VISTASOFT 
vista_soft_path = '/N/dc2/projects/lifebid/code/vistasoft/';
addpath(genpath(vista_soft_path));

%% Set the proper path for the ENCODE with Concussion model
ENCODE_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode_concussion/';
addpath(genpath(ENCODE_path));

%% HCP dataset 
dataRootPath = '/N/dc2/projects/lifebid/2t1/HCP/';
subject = '105115';
conn = 'NUM01'; % 
param = 'lmax10'; % {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8', ''}
alg = 'SD_PROB'; % {'SD_PROB', 'SD_STREAM','tensor'}

% Generate fe_strucure
%% Build the file names for the diffusion data, the anatomical MRI.
dwiFile       = deblank(ls(fullfile(dataRootPath,subject,'diffusion_data','*b2000_aligned*.nii.gz')));
t1File        = deblank(fullfile(dataRootPath,subject,'anatomy',  'T1w_acpc_dc_restore_1p25.nii.gz'));

fgFileName    = deblank(ls(fullfile(dataRootPath,subject,'fibers_new', strcat('*b2000*',char(param),'*',char(alg),'*',conn,'*','500000.tck'))));
feFileName    = strcat(subject,'_',alg,'_',param,'_',conn);  

L = 360; % Discretization parameter
Niter = 500; % Number of iterations for the optimization in LiFE

%% Initialize the model
tic
fe_prob = feConnectomeInit(dwiFile,fgFileName,feFileName,[] ,[], t1File,L,[1,0],0); % We set dwiFileRepeat =  run 02
disp(' ')
disp(['Time for model construction ','(L=',num2str(L),')=',num2str(toc),'secs']);

%% Fit the LiFE model
fe_prob = feSet(fe_prob,'fit',feFitModel(feGet(fe_prob,'model'),feGet(fe_prob,'dsigdemeaned'),'bbnnls',Niter,'preconditioner'));

%% Fit the concussion model (Isotropic + Fibers + Microstructure)
fe_prob = feBuildDictionaries(fe_prob,L,L);
fit = feFit_iso_micro(feGet(fe_prob,'model'),feGet(fe_prob,'fiber weights'), feGet(fe_prob,'bvals'), feGet(fe_prob,'bvecs'), feGet(fe_prob,'s0_img'), feGet(fe_prob,'diffusionsignalinvoxel'));
fe_prob.life.fit.w0 = fit.w0;
fe_prob.life.fit.wa = fit.wa;
fe_prob.life.fit.Qa = fit.Qa;

%save('/N/dc2/projects/lifebid/code/ccaiafa/Development/Concussion_model/scripts/variables_105115_PROB_NEW.mat','-v7.3');

%% Generate  Iso signal prediction
prediso = feGet(fe_prob,'prediso');

%% Generate Fibers signal prediction
predfibersfull = feGet(fe_prob,'predfibersfull');

%% Generate Microstructure signal prediction
predmicro = feGet(fe_prob,'predmicro');

%% Total signal prediction
predfull = prediso + predfibersfull + predmicro;

%% Compute error of fit
meas_full = feGet(fe_prob,'dsigmeasured');
norm_full = norm(meas_full,'fro');
error_full =  norm(meas_full - predfull,'fro')/norm_full;

%% Generation of NIFTI files
% Copy original NIFTI file to a local folder
ni = niftiRead(dwiFile);
fName = fullfile(dataOutputPath,strcat(subject,'_PROB_original.nii.gz'));
niftiWrite(ni,fName);

%% Iso signal
name = fullfile(dataOutputPath,strcat(subject,'_PROB_iso.nii.gz'));
Generate_nifti(ni,name,fe_prob,prediso);

%% Fibers signal
name = fullfile(dataOutputPath,strcat(subject,'_PROB_fibers.nii.gz'));
Generate_nifti(ni,name,fe_prob,predfibersfull);

%% Micro
name = fullfile(dataOutputPath,strcat(subject,'_PROB_micro.nii.gz'));
Generate_nifti(ni,name,fe_prob,predmicro);

%% Error signal
name = fullfile(dataOutputPath,strcat(subject,'_PROB_error.nii.gz'));
Generate_nifti(ni,name,fe_prob,meas_full - predfull);

%% Predfull diffusion signal
name = fullfile(dataOutputPath,strcat(subject,'_PROB_predfull.nii.gz'));
Generate_nifti(ni,name,fe_prob,predfull);

%% Generate Tract preditions
feStructurePath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/';
TractsFile = deblank(ls(char(fullfile(feStructurePath,strcat('fe_structure_*',subject,'*_STC_','*run01_500000*',alg,'*',param,'*',conn,'_TRACTS-nocull.mat')))));
load(TractsFile);
%ind = find(fe_prob.life.fit.weights);
%if length(ind)~=length(classification.index)
%    error('WARNING!: nnz weights number does not match with the numer of fibers in the tracts')
%end

for tract=1:20
    ind_tracts = find(classification.index == tract);
    % We clean all major tracs

    [fg, keep] = mbaComputeFibersOutliers(fascicles(tract),3,3);
    ind_tracts = ind_tracts(keep);
    
    %nnz_ind     = find(feGet(fe_prob,'fiber weights'));
    %ind_tracts = nnz_ind(ind_tracts);  
    predtract = feGet(fe_prob,'predfibersfull',ind_tracts);
    %% Pred tract diffusion signal
    name = fullfile(dataOutputPath,strcat(subject,'_PROB_pred_',strrep(fascicles(tract).name,' ','_'),'_NEW.nii.gz'));
    disp(['saving nifti tract ', fascicles(tract).name])
    Generate_nifti(ni,name,fe_prob,predtract);
end

rmpath(genpath(vista_soft_path));
rmpath(genpath(ENCODE_path));


 



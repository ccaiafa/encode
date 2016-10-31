
%% Load data from disk on DC2
dataRootPath = '/N/dc2/projects/lifebid/2t1/HCP/';
subject      = '105115';
conn         = 'NUM01'; % 
param        = 'lmax10'; % {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8', ''}
alg          = 'SD_PROB'; % {'SD_PROB', 'SD_STREAM','tensor'}

feStructurePath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/';
dataOutputPath  = pwd;

% Copy original NIFTI file to a local folder
dwiFile       = deblank(ls(fullfile(dataRootPath,subject,'diffusion_data','*b2000_aligned*.nii.gz')));
ni            = niftiRead(dwiFile);
fName         = fullfile(dataOutputPath,strcat(subject,'_b2000_aligned_original.nii.gz'));
niftiWrite(ni,fName);

% Load fe_structure (LiFE processed dataset)
feFileName    = deblank(ls(fullfile(feStructurePath,subject,strcat('fe_structure_',subject,'_STC_run01_',alg,'*',param,'*',conn,'.mat'))));
load(feFileName);

%% Build the LiFE model.
L = 360;

% Recompute Dictionay adding isotropic and Full dictionaries
fe.life.modelTensor = [1.5,0.1,0.1];
fe = feBuildDictionaries(fe,L,L);

% Generate predicted and measured diffusion signal (anisotropci 
% plus isotropic).
pred_full  = feGet(fe,'predfull');
meas_full  = feGet(fe,'dsigmeasured');
error_full =  norm(meas_full - pred_full,'fro')/norm(meas_full,'fro');
disp(['Error predicting full signal=',num2str(error_full)]);

%% Build new nifti file with simulated data - all fascicles validated with LiFE (pred_full)
% Initialize a signal prediction NIFTI.
ni_sim       = ni;
ni_sim.fname = fullfile(dataOutputPath,strcat(subject,'_b2000_aligned_all_fascicles_sim.nii.gz'));
coords       = fe.roi.coords;

% load dwi structure
dwi       = feGet(fe, 'dwi');
bvecs     = dwi.bvecs;
bvals     = dwi.bvals;
indexes   = find(bvals ~= 0);
b0indexes = find(bvals == 0);

% Prepare data to substitute predicted signal into NIFTI 
% Copy original S0 values
b0_data = nan(size(b0indexes,1),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = ni_sim.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end

% Initialize the new data volume with NaNs
ni_sim.data = nan(size(ni.data));

% Replace NaNs in the voxels of tract using the measured B0's 
ni_sim.data = feReplaceImageValues(ni_sim.data,b0_data,coords,b0indexes);

% Replace Nans in voxels of fiber group with dw_vals
ni_sim.data = feReplaceImageValues(ni_sim.data,pred_full,coords,indexes);

% save nifti to disk
niftiWrite(ni_sim,ni_sim.fname);


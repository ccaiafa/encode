
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
coords       = fe.roi.coords;

%% Build new nifti file with simulated data - Arcuate/CST (pred_full)
% load classification file
tract1 = 3; % CST

TractsFile = deblank(ls(char(fullfile(feStructurePath,strcat('fe_structure_*',subject,'*_STC_','*run01_500000*',alg,'*',param,'*',conn,'_TRACTS.mat')))));
load(TractsFile);
ind = find(fe.life.fit.weights);
if length(ind)~=length(classification.index)
    error('WARNING!: nnz weights number does not match with the numer of fibers in the tracts')
end

ind_tracts1 = find(classification.index==tract1);
% We clean all major tracs
for iTract = [tract1]
    fprintf('\n Cleaning %s ...',classification.names{iTract})
    [~, keep{iTract}] = mbaComputeFibersOutliers(fascicles(iTract),3,3);
end
ind_tracts1 = ind_tracts1(keep{tract1});
nnz_ind     = find(feGet(fe,'fiber weights'));
ind_tracts1 = nnz_ind(ind_tracts1);

%% Construct diffusion signal for CST only
pred_CST  = feGet(fe,'psigfiberfull',ind_tracts1);
meas_full = feGet(fe,'dsigmeasured');
error_CST =  norm(meas_full - pred_CST,'fro')/norm(meas_full,'fro');
disp(['Error predicting signal CST=',num2str(error_CST)]); 

% Compute coords for CST voxels
[inds, ~]  = find(fe.life.M.Phi(:,:,ind_tracts1)); 
voxind_CST =  unique(inds(:,2));
coords_CST = coords(voxind_CST,:);
pred_CST   = pred_CST(:,voxind_CST);

% Build new nifti file with simulated data - CST
ni_sim = ni;
ni_sim.fname = fullfile(dataOutputPath,strcat(subject,'_b2000_aligned_CST_sim.nii.gz'));

% load dwi structure
dwi   = feGet(fe, 'dwi');
bvecs = dwi.bvecs;
bvals = dwi.bvals;
indexes   = find(bvals~=0);
b0indexes = find(bvals==0);

% Copy original S0 values
b0_data = nan(size(b0indexes,1),size(coords_CST,1));
for ivx = 1:size(coords_CST,1)
    b0_data(:,ivx) = ni_sim.data(coords_CST(ivx,1),coords_CST(ivx,2),coords_CST(ivx,3),b0indexes);
end
ni_sim.data = nan(size(ni.data));

% Replace Nans in voxels of fiber group with b0_data 
ni_sim.data = feReplaceImageValues(ni_sim.data,b0_data,coords_CST,b0indexes);
% Replace Nans in voxels of fiber group with dw_vals
ni_sim.data = feReplaceImageValues(ni_sim.data,pred_CST,coords_CST,indexes);

% save nifti to disk
niftiWrite(ni_sim,ni_sim.fname);

% Extract the fiber group from the FE structure:
fg_CST = feGet(fe,'fibers acpc');
fg_CST = fgExtract(fg_CST,,'')

% Make a figure of the diffusion properties in the tract:
SF = dtiComputeDiffusionPropertiesAlongFG(fg

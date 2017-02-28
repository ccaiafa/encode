function [] = Generate_nifti(ni_in,name,fe,dwisignal)

ni_out = ni_in;
ni_out.fname = name;
% Get the coordinates of the nodes in each voxel of the connectome
coords = fe.roi.coords;

% load dwi structure
dwi = feGet(fe, 'dwi');
bvecs = dwi.bvecs;
bvals = dwi.bvals;
indexes = find(bvals~=0);
b0indexes = find(bvals==0);

% Copy original S0 values
b0_data = nan(size(b0indexes,1),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = ni_out.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end
ni_out.data = nan(size(ni_in.data));

% Replace Nans with b0_data
ni_out.data = feReplaceImageValues(ni_out.data,b0_data,coords,b0indexes);

% Replace Nans with dw_vals
ni_out.data = feReplaceImageValues(ni_out.data,dwisignal,coords,indexes);

% save nifti to disk
niftiWrite(ni_out,ni_out.fname);
end
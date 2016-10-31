function [] = DspAtoms()

encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
addpath(genpath(encode_path));

load '/N/dc2/projects/lifebid/code/ccaiafa/Development/Dictionary_learning/fe_structure_Dict_weights_105115_50iter_eps.mat'
load '/N/dc2/projects/lifebid/code/ccaiafa/Development/Dictionary_learning/results_Dict_weights_105115_50iter_eps.mat'

%Niter = size(fe.life.M.DictSig,2);
Niter = 1;

D = results{Niter}.AdpDict{Niter};
%D0 = fe.life.M.DictSig;
% load('D0.mat')
% D=D0;
orient = fe.life.M.orient;
bvecs = fe.life.bvecs;

nAtoms = size(D,2);
nDir = size(bvecs,1);
v = zeros(nDir*nAtoms,4);
vn = zeros(nDir*nAtoms,3);

figure    
hold on
box(gca,'on');
axis(gca,'tight');
i=1;
for a=1:nAtoms
    %if fe.life.M.Vox_per_atom{Niter}{a} >1000
    v((i-1)*nDir +1: i*nDir,:) = rotate_atom(D(:,i), orient(:,i), bvecs);
    vn((i-1)*nDir +1: i*nDir,:) = v((i-1)*nDir +1: i*nDir,1:3);
    vn((i-1)*nDir +1: i*nDir,:) = vn((i-1)*nDir +1: i*nDir,:).*repmat(v((i-1)*nDir +1: i*nDir,4),1,3);
    plot_atom(v((i-1)*nDir +1: i*nDir,4), v((i-1)*nDir +1: i*nDir,1:3),[1 0 0]);
    i=i+1;
    %end
end
v(i*nDir+1:end,:) = [];
vn(i*nDir+1:end,:) = [];
end

function [v] = rotate_atom(d, orient, bvecs)
%new_bvec = zeros(size(bvecs));

[Rot,~, ~] = svd(orient);
new_bvecs = bvecs*Rot;

% [theta, phi, rho] = cart2sph(bvecs(:,1),bvecs(:,2),bvecs(:,3));
% 
% [dt, dp, dr] = cart2sph(orient(1), orient(2), orient(3));
% %dp = pi/2 - dp;
% 
% % Rotate directions
% theta = theta - dt;
% phi = phi -dp;
% 
% [new_bvec(:,1), new_bvec(:,2), new_bvec(:,3),]= sph2cart(theta, phi, rho);

v = [new_bvecs, d];

end
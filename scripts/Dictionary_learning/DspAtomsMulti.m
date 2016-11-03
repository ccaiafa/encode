function [] = DspAtomsMulti()

encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
addpath(genpath(encode_path));

load '/N/dc2/projects/lifebid/code/ccaiafa/encode/scripts/Dictionary_learning/fe_structure_Dict_weights_FP.mat'
load '/N/dc2/projects/lifebid/code/ccaiafa/encode/scripts/Dictionary_learning/results_Dict_weights_FP.mat'

[results, fe] = CleanDictionaries(results, fe);

LastIter = size(results,2);
nDict = size(fe.life.M.Dictionaries,2);

for n=1:nDict
    D{n} = results{LastIter}.AdpDict{n};
end

orient = fe.life.M.orient;
bvecs = fe.life.bvecs;
MaxAtoms = 30;


for n=1:nDict
    nAtoms = size(D{n},2);
    nDir = size(bvecs,1);
    new_orient = orient(:,fe.life.M.ind_atoms{n});
    v = zeros(nDir*nAtoms,4);

    figure
    hold on
    box(gca,'on');
    axis(gca,'tight');
    title(['Dictionary Number ', num2str(n)])

    for i=1:MaxAtoms
        v((i-1)*nDir +1: i*nDir,:) = rotate_atom(D{n}(:,i), new_orient(:,i), bvecs);
        plot_atom(v((i-1)*nDir +1: i*nDir,4), v((i-1)*nDir +1: i*nDir,1:3),[1 0 0]);
        i=i+1;
        %end
    end
    v(i*nDir+1:end,:) = [];
end

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

function [results, fe] = CleanDictionaries(results, fe)
Niter = size(results,2);
nDict = size(fe.life.M.Dictionaries,2);

for n=1:nDict
    ind_vox = fe.life.M.ind_vox{n};
    Phi = fe.life.M.Phi(:,ind_vox,:);
    [sub, ~] = find(Phi);
    
    ind_atoms = unique(sub(:,1));
    
    fe.life.M.Dictionaries{n} = fe.life.M.Dictionaries{n}(:,ind_atoms);
    fe.life.M.ind_atoms{n} = ind_atoms;
    
    for iter=1:Niter
        results{iter}.AdpDict{n} = results{iter}.AdpDict{n}(:,ind_atoms);
        results{iter}.ind_atoms{n} = ind_atoms;
    end
    
end



end
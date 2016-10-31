function [] = FitTensors2Atoms()

encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
addpath(genpath(encode_path));

load 'fe_structure_Dict_weights_105115_50iter.mat'
load 'results_Dict_weights_105115_50iter.mat'

%Niter = size(fe.life.M.DictSig,2);
Niter = 9;

D = results{Niter}.AdpDict{Niter};
orient = fe.life.M.orient;
bvecs = fe.life.bvecs;
bvals = fe.life.bvals;

nAtoms = size(D,2);
nDir = size(bvecs,1);

type = 'unconstrained';

correlation = zeros(1,nAtoms);
eigenvalues = zeros(3,nAtoms);
FA = zeros(1,nAtoms);
suba = [];
for a=1:nAtoms
    alpha = norm(D(:,a));
    f = D(:,a)/alpha;
    
    [Rot,~, ~] = svd(orient(:,a)); % Compute the eigen vectors of the kernel orientation
    %Q0 = Rot*diag([1, 0.5, 0.2])*Rot';
    
    Q0 = randQa(1);
    [ Qest, mu, error ] = Estimate_atom_tensor_demeaned( f, alpha, bvecs, bvals, type ,Q0);
    [U,l] = eig(Qest,'vector');
    
    [l,ind1] = sort(l,'descend');
    U = U(:,ind1);
    
%     v1 = l(1)*U(:,1);
%     v2 = l(2)*U(:,2);
%     v3 = l(3)*U(:,3);
%     
    
%     figure
%     hold on
%     plot3([0,v1(1)],[0,v1(2)],[0,v1(3)],'Color','r','LineWidth',4);
%     plot3([0,v2(1)],[0,v2(2)],[0,v2(3)],'Color','r');
%     plot3([0,v3(1)],[0,v3(2)],[0,v3(3)],'Color','r');
    
%    vorient = orient(:,a)*l(1);
%     plot3([0,vorient(1)],[0,vorient(2)],[0,vorient(3)],'Color','b','LineWidth',4);
%     set(gca, 'DataAspectRatio',[1,1,1])
    
    correlation(a) = dot(U(:,1),orient(:,a));
    eigenvalues(:,a) = l;
    
    FA(a) = sqrt(1/2)*sqrt(((l(1)-l(2))^2 + (l(1)-l(3))^2 + (l(2)-l(3))^2)/(l(1)^2 + l(2)^2 + l(3)^2));
    
    if (FA(a)<0.9) && (FA(a) > 0.5) && (l(1) > 1) && (l(1) < 1.5) && (correlation(a) > 0.85)
        suba = [suba, a];
    end
        
    
    
    disp(['mu=',num2str(mu), ' fit error= ',num2str(100*error), ' correlation=', num2str(correlation(a))]);
%     close(gcf)
end

figure
hold on
set(gca, 'DataAspectRatio',[1,1,1])
nDir = size(bvecs,1);
v = zeros(nDir*length(suba),4);
i=1;
for a=suba
    v((i-1)*nDir +1: i*nDir,:) = rotate_atom(D(:,i), orient(:,i), bvecs);
    plot_atom(v((i-1)*nDir +1: i*nDir,4), v((i-1)*nDir +1: i*nDir,1:3),[1 0 0]);
    i=i+1;
end


end

function [Qa] = randQa(sigma)
x = sigma*randn(6,1);
X = [x(1), x(4), x(6); 0, x(2), x(5); 0, 0 ,x(3)];
Qa =  X'*X;
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
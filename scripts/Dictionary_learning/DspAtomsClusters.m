
function [] = DspAtomsClusters()

load variables.mat


Nprof = 2;
Members = zeros(Nprof,1);
[idx,C,sumd,Dout] = kmeans(Mcoeffs',Nprof);
figure
for prof=1:Nprof
    Members(prof) = sum(idx == prof);
    subplot(1,Nprof,prof);
    %plotSphFunctionCoeffs(Fnm_tdes, 'complex', 5, 5, 'real', gca); view(3), title('Reconstructed function')
    plot(C(prof,:));
end


for prof=1:Nprof
%     dirs = grid2dirs(5, 5);
%     F = inverseSHT(C(prof,:), dirs, 'real');
%     Fgrid = Fdirs2grid(F, 5, 5, 1);
%     Fgrid = Fgrid + min(Fgrid);
%     d = leastSquaresSHTeven(Lmax, Fgrid, dirs, 'real', []);
%    figure
%    plotSphFunctionCoeffs(d, 'real', 5, 5, 'real', gca); view(3), title(['atom ', num2str(prof),', ', num2str(Members(prof)),' examples'])
    
    coef = complete(mup + sigmap.*C(prof,:)', Lmax);
    %coef = C(prof,:)';
    %%subplot(1,Nprof,prof);
    figure
    plotSphFunctionCoeffs(coef, 'real', 5, 5, 'real', gca); view(3), title(['atom ', num2str(prof),', ', num2str(Members(prof)),' examples'])
    xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-0.5 0.5]); 
    
    %F = inverseSHT(coef, dirs, basisType);
end

% fe_old = load('/N/dc2/projects/lifebid/code/ccaiafa/demo_datasets/HCP3T/sub-105115/fe_structures/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat');
% 
[vals, ind_fib] = sort(fe.life.fit.weights,'descend');
ind_fib = ind_fib(1:10000);

c = jet(Nprof);

%[subs, vals] = find(fe.life.M.Phi);
%Nfib = size(fe.life.M.Phi,3);
%class_noce = zeros(Nprof, Nfib);
nAtoms = size(fe.life.M.DictSig,2);    

for prof=1:Nprof
figure
hold on
    ind = find(idx == prof);
    %ind = unique(ind);
    
    %for ind_fib=1:Nfib
        position_at = position(ind);
        ind_vox_at = [];
        for n=1:nDict
            pos = position_at((position_at >= (n-1)*nAtoms + 1) & position_at <= n*nAtoms);
            pos = mod(pos,nAtoms);
            pos(pos==0) = nAtoms;
            subtensor = fe.life.M.Phi(pos, fe.life.M.ind_vox{n}, ind_fib);
            [subs, vals] = find(subtensor);
            if ~isempty(subs)
                ind_vox_at = [ind_vox_at , fe.life.M.ind_vox{n}(unique(subs(:,2)))];
            end
        end
        
        %class_node(prof,ind_fib) = length(ind_vox_at)/nnz(fe.life.M.Phi(:,:,ind_fib))
        scatter3(fe.roi.coords(ind_vox_at,1), fe.roi.coords(ind_vox_at,2), fe.roi.coords(ind_vox_at,3), 1*ones(length(ind_vox_at),1), repmat(c(prof,:),length(ind_vox_at),1), '.')
        view(3)
    %end
end


% figure 
% hold on
% c = jet(Nprof);
% for prof=1:Nprof
%     ind = find(idx == prof);
%     ind = unique(ind);
%     position(ind )
% 
%     scatter3(fe.roi.coords(ind,1), fe.roi.coords(ind,2), fe.roi.coords(ind,3), 25*ones(length(ind),1), repmat(c(prof,:),length(ind),1), '.')
%     view(3)
% end




end


function [rho, U, mu] = correlation(d, v, bvecs, bvals)
alpha = norm(d);
f = d/alpha;

[Rot,~, ~] = svd(v); % Compute the eigen vectors of the kernel orientation


Q0 = randQa(1);
[ Qest, mu, error ] = Estimate_atom_tensor_demeaned( f, alpha, bvecs, bvals, 'unconstrained' ,Q0);
[U,l] = eig(Qest,'vector');

[l,ind1] = sort(l,'descend');
U = U(:,ind1);

rho = dot(U(:,1),v);
%u1 = U(:,1);
end

function [Qa] = randQa(sigma)
x = sigma*randn(6,1);
X = [x(1), x(4), x(6); 0, x(2), x(5); 0, 0 ,x(3)];
Qa =  X'*X;
end

function [fe] = CleanDictionaries(fe)
nDict = size(fe.life.M.Dictionaries,2);

for n=1:nDict
    
    ind_vox = fe.life.M.ind_vox{n};
    Phi = fe.life.M.Phi(:,ind_vox,:);
    [sub, ~] = find(Phi);
    
    ind_atoms = unique(sub(:,1));
    
    fe.life.M.Dictionaries{n} = fe.life.M.Dictionaries{n}(:,ind_atoms);
    fe.life.M.ind_atoms{n} = ind_atoms;
    fe.life.M.orient_atoms{n} = fe.life.M.orient(:,ind_atoms);
    
    
end
end


function [c_out] = complete(c_in, Lmax)
c_out = zeros((Lmax+1)^2,1);
nin = 1;
nout = 1;

for l=0:Lmax
    if ~mod(l,2) % even
        c_out(nout:nout+2*l) = c_in(nin:nin+2*l);
        nin = nin + 2*l + 1;
    else
        c_out(nout:nout+2*l) = zeros(2*l+1,1);
    end
    
    nout = nout + 2*l + 1;
    
end

end

function [v] = rotate_atom(d, orient, bvecs)
%new_bvec = zeros(size(bvecs));

%[Rot,~, ~] = svd(orient);
Rot = orient;
new_bvecs = bvecs*Rot;

v = new_bvecs.*repmat(d,1,3);
%v = new_bvecs;

% 
% [dt, dp, dr] = cart2sph(orient(1), orient(2), orient(3));
% %dp = pi/2 - dp;

end
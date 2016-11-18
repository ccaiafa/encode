function [] = Test_SPHT()

encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
addpath(genpath(encode_path));

%load '/N/dc2/projects/lifebid/code/ccaiafa/encode/scripts/Dictionary_learning/Experiments/Results/fe_structure_MultiHier_mult_split_FP_STC_run01_SD_PROB_lmax10_connNUM01_8Dicts_32.mat'
load '/N/dc2/projects/lifebid/code/ccaiafa/Development/Dictionary_learning/Experiments/November14/results/fe_structure_Dict_Learn_FP_tensor__connNUM01_Dicts64.mat'

aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; %convert dire from Matlab azimuth-elevetion to azimuth-inclination

% Thresholding weights
% w = fe.life.fit.weights;
% med = median(w(w>0));
% fe.life.fit.weights(w<med) = 0;

fe = CleanDictionaries(fe);


nDict = size(fe.life.M.Dictionaries,2);
D = [];
orient = [];
ind_vox = [];
for n=1:nDict
    D = [D, fe.life.M.Dictionaries{n}];
    orient = [orient, fe.life.M.orient_atoms{n}];
    ind_vox = [ind_vox, fe.life.M.ind_vox{n}];
end

bvecs = fe.life.bvecs;
bvals = fe.life.bvals;

nAtoms = size(orient,2);
nAtomsTot = size(D,2);
nDir = size(bvecs,1);

% figure    
% hold on
% box(gca,'on');
% axis(gca,'tight');
% color = [0 0 1];


Lmax = 10;
Nharm = (Lmax+1)*(Lmax+2)/2;
%Nharm = (Lmax+1)^2;
Mcoeffs = zeros(Nharm, nAtomsTot);
lambdas = zeros(3, nAtomsTot);
FA = zeros(1, nAtomsTot);
position = zeros(1, nAtomsTot);
%h1 = figure;
%h2 = figure;
hold on;
%ind_at =1;
parfor a=1:nAtomsTot
    [rho, U, v, m] = correlation(D(:,a),orient(:,a),bvecs,bvals);
    lambdas(:,a) = v;
    if abs(rho) > 0.5
        FA(a) = sqrt(0.5)*sqrt((v(1) -v(2))^2 + (v(1) -v(3))^2 + (v(2) -v(3))^2)/norm(v);        %v = rotate_atom(D(:,a), orient(:,a), bvecs);
        v = rotate_atom(D(:,a), U, bvecs);
        %scatter3(v(:,1),v(:,2),v(:,3),'.','MarkerFaceColor',color);
        [azi, elev, rho] = cart2sph(v(:,1),v(:,2),v(:,3));
        [points] = aziElev2aziIncl([azi, elev]);
        
        %figure(h1);
        %plotSphFunctionTriangle(D(:,a), points, 'real', gca); daspect([1,1,1]); title(['measurements atom ',num2str(a)]); view(3);
        
        %cond_N_reg = checkCondNumberSHT(Lmax, points, 'real', [])
        
        %regular_weights = getVoronoiWeights(aziElev2aziIncl(points));
        %cond_N_reg_weighted = checkCondNumberSHT(Lmax, points, 'real', regular_weights)
        
        % perform SHT for the regular grid using weighted least-squares and complex SHs
        Mcoeffs(:,a) = leastSquaresSHTeven(Lmax, D(:,a), points, 'real', []);
        position(a) = a;
        %norm(Mcoeffs(:,ind_at))
        %     figure(h2)
        %     plot(coeffs_reg)
        
        %Y_N = getSHeven(6, points, 'real');
        %ind_at = ind_at + 1;
        
    end
   
end
ind = find(position==0);

Mcoeffs(:,ind) = [];
position(ind) = [];
lambdas(:,ind) = [];
FA(ind) = [];

[coeff,score,latent,tsquared,explained,mu] = pca(lambdas');

Ncolors = 10;
dot_color = hsv(Ncolors);
figure
h = histogram(score(:,1),Ncolors);
ind_color = discretize(score(:,1),get(h,'BinEdges'));
%ind_color = discretize(score(:,1),[-1, -0.2, 1.5]);

figure
scatter3(lambdas(1,:), lambdas(2,:), lambdas(3,:),  1*ones(1,size(lambdas,2)), dot_color(ind_color,:))
hold on

alpha = 50;

plot3([mu(1) mu(1) + alpha*latent(1)*coeff(1,1)],[mu(2) mu(2) + alpha*latent(1)*coeff(2,1)],[mu(3) mu(3) + alpha*latent(1)*coeff(3,1)],'r', 'LineWidth',2)
plot3([mu(1) mu(1) + alpha*latent(2)*coeff(1,2)],[mu(2) mu(2) + alpha*latent(2)*coeff(2,2)],[mu(3) mu(3) + alpha*latent(2)*coeff(3,2)],'g', 'LineWidth',2)
plot3([mu(1) mu(1) + alpha*latent(3)*coeff(1,3)],[mu(2) mu(2) + alpha*latent(3)*coeff(2,3)],[mu(3) mu(3) + alpha*latent(3)*coeff(3,3)],'b', 'LineWidth',2)

xlim([0 2]); ylim([0 2]); zlim([0 2]); 

%Mcoeffs = Mcoeffs./repmat(sqrt(sum(Mcoeffs.^2,1)),Nharm,1);

load('/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_500000_SD_PROB_lmax10_connNUM01_TRACTS-nocull.mat')
%load('/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/fe_structure_105115_STC_run01_500000_tensor__connNUM01_TRACTS-nocull.mat')


tract = 3;
ind_fib = find(classification.index==tract);


% [~, ind_fib] = sort(fe.life.fit.weights,'descend');
% ind_fib = ind_fib(1:10000);

nAtoms = size(fe.life.M.DictSig,2);    
Nprof = Ncolors;

figure
hold on
for prof=1:Nprof

    ind = find(ind_color == prof);
    position_at = position(ind);
    ind_vox_at = [];
    for n=1:nDict
        pos = position_at((position_at >= (n-1)*nAtoms + 1) & position_at <= n*nAtoms);
        pos = mod(pos,nAtoms);
        pos(pos==0) = nAtoms;
        subtensor = fe.life.M.Phi(pos, fe.life.M.ind_vox{n}, ind_fib);
        [subs, ~] = find(subtensor);
        if ~isempty(subs)
            ind_vox_at = [ind_vox_at , fe.life.M.ind_vox{n}(unique(subs(:,2)))];
        end
    end
    
    scatter3(fe.roi.coords(ind_vox_at,1), fe.roi.coords(ind_vox_at,2), fe.roi.coords(ind_vox_at,3), 2*ones(length(ind_vox_at),1), repmat(dot_color(prof,:),length(ind_vox_at),1), '.')
    
    view(3)
    
    
end


%% plot FAs
clear ind_vox_at;
figure
h = histogram(FA,Ncolors);
ind_color = discretize(FA,get(h,'BinEdges'));

for prof=1:Nprof
figure
hold on
    ind = find(ind_color == prof);
    position_at = position(ind);
    ind_vox_at = [];
    for n=1:nDict
        pos = position_at((position_at >= (n-1)*nAtoms + 1) & position_at <= n*nAtoms);
        pos = mod(pos,nAtoms);
        pos(pos==0) = nAtoms;
        subtensor = fe.life.M.Phi(pos, fe.life.M.ind_vox{n}, ind_fib(:));
        [subs, ~] = find(subtensor);
        if ~isempty(subs)
            ind_vox_at = [ind_vox_at , fe.life.M.ind_vox{n}(unique(subs(:,2)))];
        end
    end
    
    scatter3(fe.roi.coords(ind_vox_at,1), fe.roi.coords(ind_vox_at,2), fe.roi.coords(ind_vox_at,3), 2*ones(length(ind_vox_at),1), repmat(dot_color(prof,:),length(ind_vox_at),1), '.')
    
    view(3)
    
    
end

%% compute statistics for individual fibers
FAfib = zeros(length(ind_fib),3);
for i = 1:length(ind_fib)
    f = ind_fib(i);
    FAnode = [];
    for n=1:nDict
        subtensor = fe.life.M.Phi(:, fe.life.M.ind_vox{n}, f);
        [subs, ~] = find(subtensor);
        
        if ~isempty(subs)
            for a = subs(:,1)'
                [rho, U, v, m] = correlation(fe.life.M.Dictionaries{n}(:,a),fe.life.M.orient_atoms{n}(:,a),bvecs,bvals);
                FAnode = [FAnode, sqrt(0.5)*sqrt((v(1) -v(2))^2 + (v(1) -v(3))^2 + (v(2) -v(3))^2)/norm(v)];
            end
        end
        

    end
    FAfib(i,:) = [length(FAnode), mean(FAnode), std(FAnode)];
    disp(['Fib ',num2str(i), ' natoms=',num2str(FAfib(i,1)), ' mean=',num2str(FAfib(i,2)), ' std=',num2str(FAfib(i,3))])
end

figure
histogram(FAfib(:,2))
title('mean FA per fiber ')

figure
histogram(FAfib(:,3))
title('std FA per fiber')

figure
histogram(FAfib(:,3)./FAfib(:,2))
title('std/mean FA per fiber')



mup = mean(Mcoeffs,2);
sigmap = std(Mcoeffs,0,2);
Mcoeffs = Mcoeffs - repmat(mup,1,size(Mcoeffs,2));
Mcoeffs = Mcoeffs./repmat(sigmap,1,size(Mcoeffs,2));


Nprof = 15;
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


c = jet(Nprof);

for prof=1:Nprof
figure
hold on
    ind = find(idx == prof);
    ind = unique(ind);
    position_at = position(ind);
    ind_vox_at = [];
    for n=1:nDict
        pos = position_at((position_at >= (n-1)*nAtoms + 1) & position_at <= n*nAtoms);
        pos = mod(pos,n*nAtoms);
        [subs, vals] = find(fe.life.M.Phi(pos, fe.life.M.ind_vox{n}, :));
        if ~isempty(subs)
            ind_vox_at = [ind_vox_at , fe.life.M.ind_vox{n}(unique(subs(:,2)))];
        end
    end

    scatter3(fe.roi.coords(ind_vox_at,1), fe.roi.coords(ind_vox_at,2), fe.roi.coords(ind_vox_at,3), 25*ones(length(ind_vox_at),1), repmat(c(prof,:),length(ind_vox_at),1), '.')
    view(3)
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


function [rho, U, l, mu] = correlation(d, v, bvecs, bvals)
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
    
    %fe.life.M.Dictionaries{n} = fe.life.M.Dictionaries{n}(:,ind_atoms);
    fe.life.M.ind_atoms{n} = ind_atoms;
    %fe.life.M.orient_atoms{n} = fe.life.M.orient(:,ind_atoms);
    fe.life.M.orient_atoms{n} = fe.life.M.orient;
    
    
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
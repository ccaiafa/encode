encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
addpath(genpath(encode_path));

load 'fe_structure_Dict_weights_105115_50iter_eps.mat'
load 'results_Dict_weights_105115_50iter_eps.mat'

Niter = 6;

full_error = zeros(1,Niter);
weights_error = zeros(1,Niter);

weights_ant = ones(feGet(fe,'nFibers'),1)/feGet(fe,'nFibers');
for i=1:Niter-1

    full_error(i) =  results{i}.error_full;
    weights_error(i) = norm(results{i}.weights{i}-weights_ant)/norm(results{i}.weights{i});
    weights_ant = results{i}.weights{i};
    
    
end

% plot erro
figure
plot(full_error)

figure
plot(weights_error,'r')


%% normalize weights
% Adjust weights according to the mean diffusion signal per fascicle
nTheta  = feGet(fe,'nbvecs');
nVoxels = feGet(fe,'nvoxels');
Phi = fe.life.M.Phi;
A = ttv(Phi,ones(size(Phi,1),1),1); % sum over atoms
A = ttv(A, ones(size(A,1),1),1); % sum over voxels
a = double(A);
a = full(a)/(nTheta*nVoxels);
Phi = sptensor(Phi.subs, Phi.vals./a(Phi.subs(:,3)),size(Phi)); % Normalize tensor Phi
fe.life.M.Phi = Phi;
fe.life.fit.weights = fe.life.fit.weights.*a;

fe_new = fe;
w_new = fe_new.life.fit.weights;
%% Load LiFE fe (paper)
dataRootPath = '/N/dc2/projects/lifebid/2t1/HCP/';
dataOutputPath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/105115/';

subject = '105115';
conn = 'NUM01'; % 

param = 'lmax10'; % {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8', ''}
alg = 'SD_PROB'; % {'SD_PROB', 'SD_STREAM','tensor'}

% load results computed for the paper
%load(fullfile(dataOutputPath,sprintf('fe_structure_%s_%s_%s_%s_conn%s.mat',subject,'STC_run01',char(alg),char(param),conn)))

load('fe_structure_105115_L90.mat');

%% normalize weights
% Adjust weights according to the mean diffusion signal per fascicle
nTheta  = feGet(fe,'nbvecs');
nVoxels = feGet(fe,'nvoxels');
Phi = fe.life.M.Phi;
A = ttv(Phi,ones(size(Phi,1),1),1); % sum over atoms
A = ttv(A, ones(size(A,1),1),1); % sum over voxels
a = double(A);
a = full(a)/(nTheta*nVoxels);
Phi = sptensor(Phi.subs, Phi.vals./a(Phi.subs(:,3)),size(Phi)); % Normalize tensor Phi
fe.life.M.Phi = Phi;
fe.life.fit.weights = fe.life.fit.weights.*a;

w = fe.life.fit.weights;

%% Compare weights histograms
edges = [-Inf -18: 0.1: -4 -2];
figure
histogram(log(w),edges,'DisplayName','LiFE weights')
hold on
histogram(log(w_new),edges,'DisplayName','New weights')
legend(gca,'show');
title({'Fiber weights histograms'});

%% Compare nnz LiFE weights values against new weights
figure
scatter(log(w_new), log(w))
xlabel({'new weights'});
ylabel({'LiFE weights'});
title({'Relation between nonzero LiFE weights and new weights'});

%% Check new weights values for fibers with zero LiFE weights
ind = find(w==0); % indices to fibers with zero LiFE weights
figure
histogram(log(w_new(ind)),edges,'DisplayName','New weights')
legend(gca,'show');
title({'Fiber new weights histogram for fibers with zero LiFE weights '});

%% Check local weights in fibers
ind = find(w); % indices to fibers with nonzero LiFE weights
w_nnz = w(ind);
[y,ind_sort] = sort(w_nnz,1,'descend');


TotalFibers = length(ind);
zero_nodes = zeros(TotalFibers,1);

for i=1:TotalFibers

   %histogram(valA);
   A = fe.life.M.Phi(:,:,ind(ind_sort(i)));
   [indA, valA] = find(A);  
   
   A_new = fe_new.life.M.Phi(:,:,ind(ind_sort(i)));
   
   valA_new = zeros(size(valA));
   for  ia=1:length(valA)
       valA_new(ia) = A_new(indA(ia,1), indA(ia,2));
   end
   
   zero_nodes(i) = length(valA_new(valA_new<1E-8))/length(valA_new);
   disp(['% Zero nodes fiber ', num2str(i),'=',num2str(100*zero_nodes(i)),'%' ])
%    figure
%    plot(valA)
%    hold on
%    plot(valA_new,'r')
   
   
end




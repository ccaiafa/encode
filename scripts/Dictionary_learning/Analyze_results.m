encode_path = '/N/dc2/projects/lifebid/code/ccaiafa/encode/';
addpath(genpath(encode_path));

% load 'fe_structure_demo_with_bbnls_weights_Dict.mat'
% load 'results_demo_with_bbnls_weights_Dict.mat'
% 
% fe_wD = fe;
% 
% Niter = size(fe.life.M.AdpDict,2);
% 
% full_error_wD = zeros(1,Niter);
% weights_error_wD = zeros(1,Niter);
% 
% weights_ant = ones(feGet(fe,'nFibers'),1)/feGet(fe,'nFibers');
% for i=1:Niter-1
% 
%     full_error_wD(i) =  results{i}.error_full;
%     weights_error_wD(i) = norm(fe.life.M.weights{i}-weights_ant)/norm(fe.life.M.weights{i});
%     weights_ant = fe.life.M.weights{i};
%     
%     
% end

load 'fe_structure_demo_with_bbnls_Dict_weights.mat'
load 'results_demo_with_bbnls_Dict_weights.mat'

Niter = size(fe.life.M.AdpDict,2);

fe_Dw = fe;

full_error_Dw = zeros(1,Niter);
weights_error_Dw = zeros(1,Niter);

weights_ant = ones(feGet(fe,'nFibers'),1)/feGet(fe,'nFibers');
for i=1:Niter-1

    full_error_Dw(i) =  results{i}.error_full;
    weights_error_Dw(i) = norm(fe.life.M.weights{i}-weights_ant)/norm(fe.life.M.weights{i});
    weights_ant = fe.life.M.weights{i};
    
    
end

% plot erro
figure
% plot(full_error_wD)
% hold on
plot(full_error_Dw)

figure
% plot(weights_error_wD)
% hold on
plot(weights_error_Dw,'r')
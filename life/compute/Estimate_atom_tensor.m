function [ mu, Qa ] = Estimate_atom_tensor( f, bvecs, bvals )

nTheta = size(bvecs,1);

% Initialization
%alpha = 100;%2*max(abs(f));
%mu = 0.5;
Qa = eye(3);
signal = exp(- bvals .* diag(bvecs*Qa*bvecs'));
mu = mean(signal);
y = signal - ones(size(signal))*mu;


%y = signal - ones(size(signal))*mu;
%alpha = sqrt(1/sum(y.^2));
%alpha = 100;


Niter = 500;

for i=1:100
    %alpha = norm(y);
    l = log(f + mu);
    % estimate tensor Qa
    Qa = NewRaph_Qa(bvecs', l, bvals(1), Niter, Qa);
    
    signal = exp(- bvals .* diag(bvecs*Qa*bvecs'));
    mu = mean(signal);
    y = signal - ones(size(signal))*mu;
    
    disp(['iter=',num2str(i),' mu=',num2str(mu), ' Qa='])
    disp(num2str(Qa(1,:)));
    disp(num2str(Qa(2,:)));
    disp(num2str(Qa(3,:)));
    
end





end


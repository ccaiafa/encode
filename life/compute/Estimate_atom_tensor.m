function [ alpha, mu, Qa ] = Estimate_atom_tensor( f, bvecs, bvals )

nTheta = size(bvecs,1);

% Initialization
%alpha = 100;%2*max(abs(f));
%mu = 0.5;
Qa = eye(3);
signal = exp(- bvals .* diag(bvecs*Qa*bvecs'));
mu = mean(signal);
y = signal - ones(size(signal))*mu;
%alpha = sqrt(1/sum(y.^2));
alpha = 100;


Niter = 500;

for i=1:100
    l = log(f/alpha + mu);
    % estimate tensor Qa
    Qa = NewRaph_Qa(bvecs', l, bvals(1), Niter, Qa);
    
    signal = exp(- bvals .* diag(bvecs*Qa*bvecs'));
    
    % compute mean mu
    mu = mean(signal);
    
    % compute alpha
    y = signal - ones(size(signal))*mu;
    %alpha = (sqrt(1/sum(y.^2)) + f'*y/sum(f.^2))/2;
    %alpha = f'*y/sum(f.^2);


    %alpha = f'*y/sum(f.^2);
    
    disp(['iter=',num2str(i),' alpha=',num2str(alpha),' mu=',num2str(mu), ' Qa='])
    disp(num2str(Qa(1,:)));
    disp(num2str(Qa(2,:)));
    disp(num2str(Qa(3,:)));
    
end





end


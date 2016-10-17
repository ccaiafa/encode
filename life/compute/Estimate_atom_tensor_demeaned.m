function [ Qa, mu, error ] = Estimate_atom_tensor_demeaned( f, alpha, bvecs, bvals, type ,Q0)

nTheta = size(bvecs,1);

% Initialization
%sigma =0.5;
%Qa = randQa(sigma);
Qa = Q0;
signal = exp(- bvals .* diag(bvecs*Qa*bvecs'));
mu = mean(signal);
% while sum(alpha*f<=-mu)
%     sigma = sigma*0.9;
%     Qa = randQa(sigma);
%     signal = exp(- bvals .* diag(bvecs*Qa*bvecs'));
%     mu = mean(signal);  
% end

%y = (signal - ones(size(signal))*mu);


Niter = 100;

delta_error = Inf;
threshold = 1e-8;
error_ant = Inf;
i=1;
while (i<Niter)&&(delta_error > threshold)

    % estimate tensor Qa
    [Qa]= NewRaph_Qa_demeaned_new(bvecs', f, bvals(1), Niter, Qa, alpha, mu, type);
    
    signal = exp(- bvals .* diag(bvecs*Qa*bvecs'));
    mu = mean(signal);
    y = (signal - ones(size(signal))*mu);
    y = y/alpha;
    error=norm(f-y)/norm(f);
    delta_error = abs(error - error_ant);
    %disp(['Iter= ', num2str(i),'   Delta Error=', num2str(delta_error)])
    
    %disp(['iter=',num2str(i),' mu=',num2str(mu),' alpha=',num2str(alpha),  ' Qa='])
    %disp(num2str(Qa(1,:)));
    %disp(num2str(Qa(2,:)));
    %disp(num2str(Qa(3,:)));
    i = i+1;
    error_ant = error;
end


error = error/norm(f);


end


function [Qa] = randQa(sigma)
x = sigma*randn(6,1);
X = [x(1), x(4), x(6); 0, x(2), x(5); 0, 0 ,x(3)];
Qa =  X'*X;
end
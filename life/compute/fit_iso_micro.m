function [w0, Qa] = fit_iso_micro(S0, b, D0, theta, S, S_f, rangewf, sumwf)
%% Initialization
Qa_ap = eye(3);

Niter = 5000; 

eror_w0 = zeros(1,Niter);
error_Qa = zeros(1,Niter);


obj = size(1,length(rangewf));

Qa_ap = cell(length(rangewf));

iter = 1;
for w0_ap = rangewf
    Qa_ap{iter} = eye(3);
    l = Compute_l(S,w0_ap,S0,b,D0,S_f,sumwf);
    [Qa_ap{iter}] = Compute_Qa(theta, l, b, 100, 1, Qa_ap{iter});
    obj(iter) = Compute_obj(S,S0,b,D0,S_f,sumwf,theta,Qa_ap{iter},w0_ap);
    iter = iter + 1;
end

[val,ind] = min(obj);

w0 = rangewf(ind);
Qa = Qa_ap{ind};

% f1 = figure;
% semilogy(rangewf,obj)
% hold on
% plot(w0,obj(ind),'.r')
% pause(0.1)
% clear f1

end

%%
function [ l ] = Compute_l(S,w0,S0,b,D0,S_f,sumwf)

Ntheta = size(S,1);

S_iso = w0*S0*exp(-b*D0);

l = S - S_iso - S_f;

l = l/(S0*( 1 - w0 - sumwf));

l(l<=0) = 1e-3;

l = log(l);

end


%% 
function [ Qa_ap ] = Compute_Qa( theta, l, b, Niter, nu , Qa_ap)
X = chol(Qa_ap);
x = [X(1,1),X(2,2),X(3,3),X(1,2),X(2,3),X(1,3)]';
grad = zeros(6,1);
dH = zeros(6,6);

lambda = 0;
flag =1;


normgrad = [];
error = [];
obj = [];
cond = 0;
Qa_ant = Qa_ap;
iter = 1;
while (cond==0) && (iter < Niter)
    
    D = X*theta;   
    m = l + b*(dot(D,D))';
    
    if flag
        grad(1) = sum(4*b*(m.*D(1,:)'.*theta(1,:)'));
        grad(2) = sum(4*b*(m.*D(2,:)'.*theta(2,:)'));
        grad(3) = sum(4*b*(m.*D(3,:)'.*theta(3,:)'));
        grad(4) = sum(4*b*(m.*D(1,:)'.*theta(2,:)'));
        grad(5) = sum(4*b*(m.*D(2,:)'.*theta(3,:)'));
        grad(6) = sum(4*b*(m.*D(1,:)'.*theta(3,:)'));

        dH(1,1) = 4*b*sum(2*b*(D(1,:).^2)'.*(theta(1,:).^2)' + m.*(theta(1,:).^2)');
        dH(1,2) = 4*b*sum(2*b*D(2,:)'.*D(1,:)'.*theta(2,:)'.*theta(1,:)');
        dH(1,3) = 4*b*sum(2*b*D(3,:)'.*D(1,:)'.*theta(3,:)'.*theta(1,:)');
        dH(1,4) = 4*b*sum(2*b*(D(1,:).^2)'.*theta(2,:)'.*theta(1,:)' + m.*theta(1,:)'.*theta(2,:)');
        dH(1,5) = 4*b*sum(2*b*D(2,:)'.*D(1,:)'.*theta(3,:)'.*theta(1,:)');
        dH(1,6) = 4*b*sum(2*b*(D(1,:).^2)'.*theta(3,:)'.*theta(1,:)' + m.*theta(3,:)'.*theta(1,:)');

        dH(2,1) = dH(1,2);
        dH(2,2) = 4*b*sum(2*b*(D(2,:).^2)'.*(theta(2,:).^2)' + m.*(theta(2,:).^2)');
        dH(2,3) = 4*b*sum(2*b*D(3,:)'.*D(2,:)'.*theta(3,:)'.*theta(2,:)');
        dH(2,4) = 4*b*sum(2*b*D(1,:)'.*D(2,:)'.*theta(2,:)'.*theta(2,:)');
        dH(2,5) = 4*b*sum(2*b*(D(2,:).^2)'.*theta(3,:)'.*theta(2,:)' + m.*theta(3,:)'.*theta(2,:)');
        dH(2,6) = 4*b*sum(2*b*D(1,:)'.*D(2,:)'.*theta(3,:)'.*theta(2,:)');

        dH(3,1) = dH(1,3);
        dH(3,2) = dH(3,2);
        dH(3,3) = 4*b*sum(2*b*(D(3,:).^2)'.*(theta(3,:).^2)' + m.*(theta(3,:).^2)');
        dH(3,4) = 4*b*sum(2*b*D(1,:)'.*D(3,:)'.*theta(2,:)'.*theta(3,:)');
        dH(3,5) = 4*b*sum(2*b*D(2,:)'.*D(3,:)'.*theta(3,:)'.*theta(3,:)');
        dH(3,6) = 4*b*sum(2*b*D(1,:)'.*D(3,:)'.*theta(3,:)'.*theta(3,:)');

        dH(4,1) = dH(1,4);
        dH(4,2) = dH(2,4);
        dH(4,3) = dH(3,4);
        dH(4,4) = 4*b*sum(2*b*(D(1,:).^2)'.*(theta(2,:).^2)' + m.*(theta(2,:).^2)');
        dH(4,5) = 4*b*sum(2*b*D(2,:)'.*D(1,:)'.*theta(3,:)'.*theta(2,:)');
        dH(4,6) = 4*b*sum(2*b*(D(1,:).^2)'.*theta(3,:)'.*theta(2,:)' + m.*theta(3,:)'.*theta(2,:)');

        dH(5,1) = dH(1,5);
        dH(5,2) = dH(2,5);
        dH(5,3) = dH(3,5);
        dH(5,4) = dH(4,5);
        dH(5,5) = 4*b*sum(2*b*(D(2,:).^2)'.*(theta(3,:).^2)' + m.*(theta(3,:).^2)');
        dH(5,6) = 4*b*sum(2*b*D(1,:)'.*D(2,:)'.*theta(3,:)'.*theta(3,:)');

        dH(6,1) = dH(1,6);
        dH(6,2) = dH(2,6);
        dH(6,3) = dH(3,6);
        dH(6,4) = dH(4,6);
        dH(6,5) = dH(5,6);
        dH(6,6) = 4*b*sum(2*b*(D(1,:).^2)'.*(theta(3,:).^2)' + m.*(theta(3,:).^2)');
    end
    
    % modified newton iteration
    
    delta = (dH+lambda*eye(6))\(-grad);
    xn = x + delta;
    Xn = [xn(1), xn(4), xn(6); 0, xn(2), xn(5); 0, 0 ,xn(3)];
    Dn = Xn*theta;
    
    if sum((l + b*(dot(Dn,Dn))').^2) < sum((l + b*(dot(D,D))').^2)
        lambda = 0.1*lambda;
        x = xn;
        D = Dn;
        
        flag = 1;
        
        X = [x(1), x(4), x(6); 0, x(2), x(5); 0, 0 ,x(3)];
    
        Qa_ap = X'*X;
        relchange = norm(Qa_ant - Qa_ap)/norm(Qa_ap);
        cond = relchange < 1e-4;
        %error = [error, norm(Qa-Qa_ap)];
        %normgrad = [normgrad, norm(grad)];
        obj = [obj, sum((l + b*(dot(D,D))').^2)];

        %plot(obj)

        %pause

        Qa_ant = Qa_ap;
    
    else
        if (lambda==0)
            lambda = 0.00001;
        else
            lambda = 10*lambda;
        end
        flag = 0;
    end
    
    
    %x = x - nu*pinv(dH)*grad;
    
    
    iter = iter + 1;
end

 if nnz(eigs(Qa_ap)<=0)
    disp('error: Qa_ap not positive definite')
    %Qa_ap = eye(3);
    Qa_ap = NaN;
end


end

%%
function [ obj ] = Compute_obj(S,S0,b,D0,S_f,sumwf,theta,Qa,w0)

Ntheta = size(theta,2);
obj = zeros(Ntheta,1);

obj = S;
for i=1:Ntheta
    % Isotropic term
    S_iso = w0*S0*exp(-b*D0); 
    
    
    % Microstructure (anisotropic)
    wa = 1 - w0 - sumwf;
    S_a = wa*S0*exp(-b*theta(:,i)'*Qa*theta(:,i));
    
    obj(i) = S(i) - S_iso - S_f(i) - S_a;   
end
obj = sum(obj.^2);
end



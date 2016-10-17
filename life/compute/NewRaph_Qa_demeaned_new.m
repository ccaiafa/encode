%% This function fit a tensor matrix Qa such that l = -b*theta'*Qa*theta
function [ Qa_ap ] = NewRaph_Qa_demeaned_new( theta, f, b, Niter, Qa_ap, alpha, mu, type)
% check if Qa_ap is positive definite

lambda = eig(Qa_ap);
while ~isempty(lambda(lambda < 1e-4))
    Qa_ap=Qa_ap+0.1*eye(3);
    lambda = eig(Qa_ap);
end
X = chol(Qa_ap);
x = [X(1,1),X(2,2),X(3,3),X(1,2),X(2,3),X(1,3)]';
grad = zeros(6,1);
dH = zeros(6,6);

lambda = 0;
flag =1;
cond = 0;
Qa_ant = Qa_ap;
iter = 1;
while (cond==0) && (iter < Niter)
    
    l = log(mu + alpha*f);
    while ~isreal(l)
        mu = -1.1*alpha*min(f);
        l = log(mu + alpha*f);
    end
    
    D = X*theta;   
    m = l + b*(dot(D,D))';
    
    if flag
        %factor = f./(mu*ones(size(f)) + alpha*f);
        
        grad(1) = sum(4*b*(m.*D(1,:)'.*theta(1,:)'));
        grad(2) = sum(4*b*(m.*D(2,:)'.*theta(2,:)'));
        grad(3) = sum(4*b*(m.*D(3,:)'.*theta(3,:)'));
        grad(4) = sum(4*b*(m.*D(1,:)'.*theta(2,:)'));
        grad(5) = sum(4*b*(m.*D(2,:)'.*theta(3,:)'));
        grad(6) = sum(4*b*(m.*D(1,:)'.*theta(3,:)'));
        %grad(7) = sum(2*m.*factor);

        dH(1,1) = 4*b*sum(2*b*(D(1,:).^2)'.*(theta(1,:).^2)' + m.*(theta(1,:).^2)');
        dH(1,2) = 4*b*sum(2*b*D(2,:)'.*D(1,:)'.*theta(2,:)'.*theta(1,:)');
        dH(1,3) = 4*b*sum(2*b*D(3,:)'.*D(1,:)'.*theta(3,:)'.*theta(1,:)');
        dH(1,4) = 4*b*sum(2*b*(D(1,:).^2)'.*theta(2,:)'.*theta(1,:)' + m.*theta(1,:)'.*theta(2,:)');
        dH(1,5) = 4*b*sum(2*b*D(2,:)'.*D(1,:)'.*theta(3,:)'.*theta(1,:)');
        dH(1,6) = 4*b*sum(2*b*(D(1,:).^2)'.*theta(3,:)'.*theta(1,:)' + m.*theta(3,:)'.*theta(1,:)');
        %dH(1,7) = sum(4*b*(m.*D(1,:)'.*theta(1,:)'.*factor));

        dH(2,1) = dH(1,2);
        dH(2,2) = 4*b*sum(2*b*(D(2,:).^2)'.*(theta(2,:).^2)' + m.*(theta(2,:).^2)');
        dH(2,3) = 4*b*sum(2*b*D(3,:)'.*D(2,:)'.*theta(3,:)'.*theta(2,:)');
        dH(2,4) = 4*b*sum(2*b*D(1,:)'.*D(2,:)'.*theta(2,:)'.*theta(2,:)');
        dH(2,5) = 4*b*sum(2*b*(D(2,:).^2)'.*theta(3,:)'.*theta(2,:)' + m.*theta(3,:)'.*theta(2,:)');
        dH(2,6) = 4*b*sum(2*b*D(1,:)'.*D(2,:)'.*theta(3,:)'.*theta(2,:)');
        %dH(2,7) = sum(4*b*(m.*D(2,:)'.*theta(2,:)'.*factor));

        dH(3,1) = dH(1,3);
        dH(3,2) = dH(3,2);
        dH(3,3) = 4*b*sum(2*b*(D(3,:).^2)'.*(theta(3,:).^2)' + m.*(theta(3,:).^2)');
        dH(3,4) = 4*b*sum(2*b*D(1,:)'.*D(3,:)'.*theta(2,:)'.*theta(3,:)');
        dH(3,5) = 4*b*sum(2*b*D(2,:)'.*D(3,:)'.*theta(3,:)'.*theta(3,:)');
        dH(3,6) = 4*b*sum(2*b*D(1,:)'.*D(3,:)'.*theta(3,:)'.*theta(3,:)');
        %dH(3,7) = sum(4*b*(m.*D(3,:)'.*theta(3,:)'.*factor));

        dH(4,1) = dH(1,4);
        dH(4,2) = dH(2,4);
        dH(4,3) = dH(3,4);
        dH(4,4) = 4*b*sum(2*b*(D(1,:).^2)'.*(theta(2,:).^2)' + m.*(theta(2,:).^2)');
        dH(4,5) = 4*b*sum(2*b*D(2,:)'.*D(1,:)'.*theta(3,:)'.*theta(2,:)');
        dH(4,6) = 4*b*sum(2*b*(D(1,:).^2)'.*theta(3,:)'.*theta(2,:)' + m.*theta(3,:)'.*theta(2,:)');
        %dH(4,7) = sum(4*b*(m.*D(1,:)'.*theta(2,:)'.*factor));

        dH(5,1) = dH(1,5);
        dH(5,2) = dH(2,5);
        dH(5,3) = dH(3,5);
        dH(5,4) = dH(4,5);
        dH(5,5) = 4*b*sum(2*b*(D(2,:).^2)'.*(theta(3,:).^2)' + m.*(theta(3,:).^2)');
        dH(5,6) = 4*b*sum(2*b*D(1,:)'.*D(2,:)'.*theta(3,:)'.*theta(3,:)');
        %dH(5,7) = sum(4*b*(m.*D(2,:)'.*theta(3,:)'.*factor));

        dH(6,1) = dH(1,6);
        dH(6,2) = dH(2,6);
        dH(6,3) = dH(3,6);
        dH(6,4) = dH(4,6);
        dH(6,5) = dH(5,6);
        dH(6,6) = 4*b*sum(2*b*(D(1,:).^2)'.*(theta(3,:).^2)' + m.*(theta(3,:).^2)');
        %dH(6,7) = sum(4*b*(m.*D(1,:)'.*theta(3,:)'.*factor));
        
%         dH(7,1) = dH(1,7);
%         dH(7,2) = dH(2,7);
%         dH(7,3) = dH(3,7);
%         dH(7,4) = dH(4,7);
%         dH(7,5) = dH(5,7);
%         dH(7,6) = dH(6,7);
        %dH(7,7) = sum(2*(factor.^2).*(ones(size(factor)) - m));
    end
    
    % modified newton iteration
    delta = (dH+lambda*eye(6))\(-grad);
    xn = x + delta;
    Xn = [xn(1), xn(4), xn(6); 0, xn(2), xn(5); 0, 0 ,xn(3)];
    %alpha_n = max(xn(7),-xn(7));
    Dn = Xn*theta;
    mu_n = mean(exp(-  b*(dot(Dn,Dn))));
    ln = log(mu_n + alpha*f);
    while ~isreal(ln)
        mu_n = -1.1*alpha*min(f);
        ln = log(mu_n + alpha*f);
    end
    %if ~isreal(ln)
    %    keyboard
    %end
    if (sum((ln + b*(dot(Dn,Dn))').^2) < sum((l + b*(dot(D,D))').^2))
        lambda = 0.1*lambda;
        x = xn;
        D = Dn;
        
        flag = 1;
        
        X = Xn;
        %alpha_ap = alpha_n;
        mu = mu_n;
        Qa_ap = X'*X;
        relchange = norm(Qa_ant - Qa_ap)/norm(Qa_ap);
        cond = relchange < 1e-8;
        %error = [error, norm(Qa-Qa_ap)];
        %normgrad = [normgrad, norm(grad)]; 
        %obj = [obj, sum((ln + b*(dot(Dn,Dn))').^2)];
        %disp(['Iter=', num2str(iter), 'Obj=',num2str(obj(end))])
        %Qa_ap

        %plot(obj)

        %pause(0.1)

        Qa_ant = Qa_ap;
        %l = log(mu + alpha_ap*f);
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

if strcmp(type,'constrained')
    Qa_ap = Qa_ap/norm(Qa_ap);
end

if nnz(eigs(Qa_ap)<=0)
    %disp('error: Qa_ap not positive definite')
    %Qa_ap = eye(3);
    Qa_ap = Qa_ap + 0.1*eye(3);
end


end
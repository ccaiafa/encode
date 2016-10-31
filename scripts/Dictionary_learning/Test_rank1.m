% Test rank-1 left non-negative factor retrieval
function[] = Test_rank1
load('D0.mat')

da = D0(:,1);
ba = rand(90,1);

sigma = 0.1;
Noise = randn(size(da,1),size(ba,1));

E =  da*ba' + sigma*Noise;

[U,S,V] = svd(E);

da_ap = U(:,1);
%ba_ap = S(1,1)*V(:,1);

ba_ap = mean(E,1)/mean(da_ap);

if dot(da_ap,da) < 0
    da_ap = -da_ap;
    ba_ap = - ba_ap;
end

ba_ap_norm = ba_ap/norm(ba_ap);
ba_norm = ba/norm(ba);
da_ap_norm = da_ap/norm(da_ap);
da_norm = da/norm(da);

figure
plot(ba_norm)
hold on
plot(ba_ap_norm)

figure
plot(da_norm)
hold on
plot(da_ap_norm)

error_da_svd = norm(da - da_ap)
error_ba_svd = norm(ba - ba_ap)
error_svd = norm(E - da_ap*ba_ap','fro')/norm(E,'fro')

lambda = 0;

[da_ap, ba_ap] = right_nn_svds(E, lambda);
ba_ap = ba_ap';
%ba_ap = ba;
ba_ap_norm = ba_ap/norm(ba_ap);
da_ap_norm = da_ap/norm(da_ap);

figure
plot(ba_norm)
hold on
plot(ba_ap_norm)

figure
plot(da_norm)
hold on
plot(da_ap_norm)

error_da_nnsvd = norm(da - da_ap)
error_ba_nnsvd = norm(ba - ba_ap)
error_svd = norm(E - da_ap*ba_ap','fro')/norm(E,'fro')

end



function [da,ba] = right_nn_svds(E,lambda)
tol = 1e-8;

%[u1,s1,v1] = svds(E,1);
v1 = ones(size(E,2),1);
v1 = v1/norm(v1);
u1 = E*v1/(v1'*v1);
s1 = norm(E);
% try positive sing vectors
da = u1;
ba = s1*v1;
ba(ba<0) = eps;
error = norm(E-da*ba','fro');
delta = Inf;
normE = norm(E,'fro');
while delta > tol && sum(ba)
    ba = da'*E/(da'*da+2*lambda);
    ba(ba<0) = eps;
    if sum(ba) 
        da = E*ba'/(ba*ba');
    end
    errorn = norm(E-da*ba,'fro')/normE;
    delta = abs(errorn - error);
    error = errorn;
end
error_pos = error;
%error_pos = norm(atom - da);
da_pos = da;
ba_pos = ba;

% try neg sing vectors
da = -u1;
ba = -s1*v1;
ba(ba<0) = eps;
error = norm(E-da*ba','fro');
delta = Inf;
while delta > tol && sum(ba) 
    ba = da'*E/(da'*da+2*lambda);
  
    ba(ba<0) = eps;
    
    if sum(ba) 
        da = E*ba'/(ba*ba');
    end
    errorn = norm(E-da*ba,'fro')/normE;
    delta = abs(errorn - error);
    error = errorn;
end
error_neg = error;
%error_neg = norm(atom - da);
da_neg = da;
ba_neg = ba;

if error_pos < error_neg
    da = da_pos;
    ba = ba_pos;    
else
    da = da_neg;
    ba = ba_neg;
end

lambda = norm(da);
if lambda
    da = da/lambda;
    ba = ba*lambda;
end

end



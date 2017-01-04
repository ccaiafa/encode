function [D, B] = update_DandB(dSig, D, B,lambda)
% Compute DB
DB = D*B;
nAtoms = size(D,2);
[ind1, ind2] = find(B);
%% Update B(a,v) and Dict(theta,a)
%h = waitbar(0,'Udapting atoms ...');
disp('Updating atoms...')
for a=1:nAtoms
    %waitbar(a/nAtoms, h,['Adapting atoms (',num2str(a),'/',num2str(nAtoms),') ...']);
    pos = find(ind1==a);
    if ~isempty(pos)
        
        DB = DB - D(:,a)*B(a,:);
        E = dSig - DB;
        cols = ind2(pos);
        E = E(:,cols);
        [da,ba] = right_nn_svds(E, lambda);

        % update Dict and B
        D(:,a) = da;
        B(a,cols) = ba';
        
        % uptdate DB
        DB = DB + D(:,a)*B(a,:);
        
        %error_post = norm(dSig-D*B,'fro')/norm(dSig,'fro')
    end
    
    %disp(['atom ',num2str(a),' of ',num2str(nAtoms)]);
end
%close(h)

%disp(['Fit Dic and B, iter',num2str(iter),' error=',num2str(100*norm(dSig-M.DictSig*B,'fro')/norm(dSig,'fro'))])

end


function [da,ba] = right_nn_svds(E,lambda)
tol = 1e-4;

[u1,s1,v1] = svds(E,1);
% [U,S,V] = svd(E,'econ');
% u1 = U(:,1);
% v1 = V(:,1);
% s1 = S(1,1);
% try positive sing vectors
da = u1;
ba = s1*v1;

ba(ba<0) = eps;

normE = norm(E,'fro');
error = norm(E-da*ba','fro')/normE;
delta = Inf;
while delta > tol && sum(ba)
    %ba = da'*E/(da'*da+2*lambda);
    ba = E'*da/(da'*da+2*lambda);
    ba(ba<0) = eps;
    if sum(ba) 
        %da = E*ba'/(ba*ba');
        da = E*ba/(ba'*ba);
    end
    %errorn = norm(E-da*ba,'fro')/normE;
    errorn = norm(E-da*ba','fro')/normE;
    delta = abs(errorn - error);
    error = errorn;
end
error_pos = error;
da_pos = da;
ba_pos = ba;

% try neg sing vectors
da = -u1;
ba = -s1*v1;
ba(ba<0) = eps;
error = norm(E-da*ba','fro')/normE;
delta = Inf;
while delta > tol && sum(ba) 
    %ba = da'*E/(da'*da+2*lambda);
    ba = E'*da/(da'*da+2*lambda);
    ba(ba<0) = eps;
    if sum(ba) 
        %da = E*ba'/(ba*ba');
        da = E*ba/(ba'*ba);
    end
    %errorn = norm(E-da*ba,'fro')/normE;
    errorn = norm(E-da*ba','fro')/normE;
    delta = abs(errorn - error);
    error = errorn;
end
error_neg = error;
da_neg = da;
ba_neg = ba;

%disp([num2str(error_pos),'  ',num2str(error_neg)])
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

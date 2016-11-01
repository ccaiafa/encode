function [fe] = feInitFitDictionaries(varargin)
fe = varargin{1};
lambda = varargin{2};


% feFitPhi() function that given a fix vector of weights w optimize the tensor Phi
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

%load variables.mat

%[nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
[nAtoms] = feGet(fe,'natoms');
[nTheta] = size(fe.life.M.DictSig,1);
[nVoxels] = feGet(fe,'nvoxels');

w = fe.life.fit.weights;

dSig = reshape(feGet(fe,'dsigdemeaned'),[nTheta,nVoxels]);

B = ttv(fe.life.M.Phi,w,3);
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);

ind_vox = 1:nVoxels;
level = 1;
epsilon = 0.3;
while (length(ind_vox) > 1000) && level < 8
    [Dict, B(:,ind_vox)] = update(dSig, fe.life.M.DictSig, B(:,ind_vox), lambda);
    e = sum((dSig - Dict*B(:,ind_vox)).^2,1)./sum(dSig.^2,1);
    ind_voxA = find(e < epsilon); % indices to voxels with low error
    fe.life.M.Dictionaries{level} = Dict;
    fe.life.M.ind_vox{level} = ind_vox(ind_voxA);
    ind_vox = find(e >= epsilon); %indices to voxels with high error need to be fitted with a new dictionary (next level)
    dSig = dSig(:,ind_vox);
    level = level + 1;
    
end

%% Compute Phi compatible with B(a,v)
[sub, ~] = find(fe.life.M.Phi);
ind = sub2ind(size(B), sub(:,1), sub(:,2));
b = B(:);

A = sptensor(sub, ones(size(sub,1),1), size(fe.life.M.Phi));
A = ttv(A,w.^2,3);
[subA, val] = find(A);
A = sparse(subA(:,1),subA(:,2),val,nAtoms,nVoxels);
a = A(:);
div = a(ind);
[ind1,~,val] = find(div);

newval = full(w(sub(:,3)).*b(ind));
newval(ind1) = newval(ind1)./val;

fe.life.M.Phi = sptensor(sub, newval, size(fe.life.M.Phi));

end

function [D, B] = update(dSig, D, B,lambda)
% Compute DB
DB = D*B;
nAtoms = size(D,2);
[ind1, ind2] = find(B);
%% Update B(a,v) and Dict(theta,a)
h = waitbar(0,'Adapting atoms ...');
for a=1:nAtoms
    waitbar(a/nAtoms, h,['Adapting atoms (',num2str(a),'/',num2str(nAtoms),') ...']);
    pos = find(ind1==a);
    if ~isempty(pos)
        
        DB = DB - D(:,a)*B(a,:);
        E = dSig - DB;
        cols = ind2(pos);
        E = E(:,cols);
        [da,ba] = right_nn_svds(E, lambda);
        
        % update Dict and B
        D(:,a) = da;
        B(a,cols) = ba;
        
        % uptdate DB
        DB = DB + D(:,a)*B(a,:);
        
        %error_post = norm(dSig-D*B,'fro')/norm(dSig,'fro')
    end
    
    %disp(['atom ',num2str(a),' of ',num2str(nAtoms)]);
end
close(h)

%disp(['Fit Dic and B, iter',num2str(iter),' error=',num2str(100*norm(dSig-M.DictSig*B,'fro')/norm(dSig,'fro'))])

end


function [da,ba] = right_nn_svds(E,lambda)
tol = 1e-4;

[u1,s1,v1] = svds(E,1);
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

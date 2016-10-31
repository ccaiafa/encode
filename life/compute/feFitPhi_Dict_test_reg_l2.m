function [Phi, Dict, Vox_per_atom] = feFitPhi_Dict_test_reg_l2(varargin)
M = varargin{1};
w = varargin{2};
dSig = varargin{3};
Niter = varargin{4};
bvecs = varargin{5};
bvals = varargin{6};
%preconditioner = varargin{5};

% feFitPhi() function that given a fix vector of weights w optimize the tensor Phi
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

%load variables.mat

[nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
[nAtoms] = size(M.DictSig,2); %feGet(fe,'natoms');
[nTheta] = size(M.DictSig,1);
[nVoxels] = size(M.Phi,2); %feGet(fe,'nvoxels');

dSig = reshape(dSig,[nTheta,nVoxels]);

B = ttv(M.Phi,w,3);
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);

for iter=1:Niter
    % Compute DB
    DB = M.DictSig*B;
    
    %% Update B(a,v) and Dict(theta,a)
    h = waitbar(0,'Adapting atoms ...');
    for a=1:nAtoms
        waitbar(a/nAtoms, h,['Adapting atoms (',num2str(a),'/',num2str(nAtoms),') ...']);
        
        pos = find(ind(:,1)==a);
        Vox_per_atom{a} = pos;
        if ~isempty(pos)      
            
            DB = DB - M.DictSig(:,a)*B(a,:);
            E = dSig - DB;
            cols = ind(pos,2);     
            E = E(:,cols);
            [da,ba] = right_nn_svds(E, M.DictSig(:,a));  
            
            % update Dict and B
            M.DictSig(:,a) = da;
            B(a,cols) = ba;           
            
            % uptdate DB
            DB = DB + M.DictSig(:,a)*B(a,:);              

            %error_post = norm(dSig-M.DictSig*B,'fro')/norm(dSig,'fro')
        end
        
        %disp(['atom ',num2str(a),' of ',num2str(nAtoms)]);
    end
    close(h)
    
    disp(['Fit Dic and B, iter',num2str(iter),' error=',num2str(100*norm(dSig-M.DictSig*B,'fro')/norm(dSig,'fro'))])
end

%% Compute Phi compatible with B(a,v)
[sub, ~] = find(M.Phi);
ind = sub2ind(size(B), sub(:,1), sub(:,2));
b = B(:);

A = sptensor(sub, ones(size(sub,1),1), size(M.Phi));
A = ttv(A,w.^2,3);
[subA, val] = find(A);
A = sparse(subA(:,1),subA(:,2),val,nAtoms,nVoxels);
a = A(:);
div = a(ind);
[ind1,~,val] = find(div);

newval = full(w(sub(:,3)).*b(ind));
newval(ind1) = newval(ind1)./val;

Phi = sptensor(sub, newval, size(M.Phi));
Dict = M.DictSig;


end



function [da,ba] = right_nn_svds(E,atom)
tol = 1e-2;
atom = atom/norm(atom);
lambda = 0.5;

% E(isnan(E))=0;
% E(isinf(E))=0;
[u1,s1,v1] = svds(E,1);
% try positive sing vectors
da = u1;
ba = s1*v1;
m = mean(ba(ba>0));
if ~isnan(m)
    ba(ba<0) = m;
else
    ba(ba<0) = 0;
end
error = norm(E-da*ba','fro');
delta = Inf;
normE = norm(E,'fro');
while delta > tol && sum(ba)
    ba = da'*E/(1+2*lambda);
    m = mean(ba(ba>0));
    if ~isnan(m)
        ba(ba<0) = m;
    else
        ba(ba<0) = 0;
    end
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
ba(ba<0) = mean(ba(ba>0));
error = norm(E-da*ba','fro');
delta = Inf;
while delta > tol && sum(ba) 
    ba = da'*E/(1+2*lambda);
  
    m = mean(ba(ba>0));
    if ~isnan(m)
        ba(ba<0) = m;
    else
        ba(ba<0) = 0;
    end
    
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

% if error_pos < 0.1*error_neg
%     da = da_pos;
%     ba = ba_pos;
% else
%     if error_neg < 0.1*error_pos
%         da = da_neg;
%         ba = ba_neg;
%     else
%         if atom'*da_pos/norm(da_pos) > atom'*da_neg/norm(da_neg)
%             da = da_pos;
%             ba = ba_pos;
%         else
%             da = da_neg;
%             ba = ba_neg;
%         end
%     end
% 
% end
end

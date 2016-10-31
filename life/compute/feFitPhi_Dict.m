function [Phi, Dict, Vox_per_atom] = feFitPhi_Dict(varargin)
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

Qtensors = {};
Atom_mean = zeros(nAtoms,1);
Tensor_fit_error = zeros(nAtoms,1);

dSig = reshape(dSig,[nTheta,nVoxels]);


for iter=1:Niter
    B = ttv(M.Phi,w,3);
    [ind, val] = find(B);
    B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
    
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  plot_atom(M.DictSig(:,a),bvecs,[1 0 0]);
%                                  hold on
%                                  %plot_atom(M.DictFull(:,a),bvecs,[0 1 0]);
%                                  plot_atom(da,bvecs,[0 0 1]);
%                                  title(['Nr of voxels=',num2str(length(ba))]);
%             
%                                  plot3([0,M.orient(1,a)],[0,M.orient(2,a)],[0,M.orient(3,a)],'Color','r','LineWidth',4)
%                                  set(gca, 'DataAspectRatio',[1,1,1])
% plot(M.DictSig(:,a))
% hold on
% plot(da,'r')
            
            
            % Fit a demeaned tensor to the atom
            %                     alpha = norm(M.DictSig(:,a));
            %                     f = M.DictSig(:,a)/alpha;
            %                     [ Qest, mu, error ] = Estimate_atom_tensor_demeaned( f, alpha, bvecs, bvals);
            %                     %[U,S,V] = svd(Qest);
            %                     [U,l] = eig(Qest,'vector');
            %
            %                     [l,ind1] = sort(l,'descend');
            %                     U = U(:,ind1);
            %
            %
            %                     v1 = l(1)*U(:,1);
            %                     v2 = l(2)*U(:,2);
            %                     v3 = l(3)*U(:,3);
            %
            %                     plot3([0,v1(1)],[0,v1(2)],[0,v1(3)],'Color','r','LineWidth',4);
            %                     plot3([0,v2(1)],[0,v2(2)],[0,v2(3)],'Color','r');
            %                     plot3([0,v3(1)],[0,v3(2)],[0,v3(3)],'Color','r');
            %
            %                     FA = sqrt(1/2)*sqrt(((l(1)-l(2))^2 + (l(1)-l(3))^2 + (l(2)-l(3))^2)/(l(1)^2 + l(2)^2 + l(3)^2));
            %                     disp(['Original Dict mu=',num2str(mu), ' fit error= ',num2str(100*error),'%','  FA=',num2str(FA)]);
            %                     disp(['Lambdas=',num2str(l')])
            
            
            % Fit a demeaned tensor to the adapted atom
%             alpha = std(da);
%             f = da/alpha;
%             
%             [Rot,~, ~] = svd(M.orient(:,a)); % Compute the eigen vectors of the kernel orientation
%             Q0 = Rot*diag([1, 0.5, 0.2])*Rot';
%             %[ Qest, mu, error ] = Estimate_atom_tensor_demeaned( f, alpha, bvecs, bvals,'unconstrained',Q0);
%             [ Qest, ~ ] = NewRaph_Qa_demeaned_alpha_mu( bvecs', f, bvals(1), 100, Q0, alpha, 'unconstrained');
%             mu = mean(exp(- bvals .* diag(bvecs*Qest*bvecs')));
%             error = norm(alpha*f -  exp(- bvals .* diag(bvecs*Qest*bvecs')) + mean(exp(- bvals .* diag(bvecs*Qest*bvecs'))))/norm(alpha*f);
%             
%             
%             Qtensors{a} = Qest;
%             Atom_mean(a) = mu;
%             Tensor_fit_error(a) = error;
            
            %disp(['Atom ',num2str(a),'/',num2str(nAtoms)])
            
            %[U,S,V] = svd(Qest);
            %                     [U,l] = eig(Qest,'vector');
            %
            %                     [l,ind1] = sort(l,'descend');
            %                     U = U(:,ind1);
            
            %                     ang_error_unconst(a) = acos(abs(U(:,1)'*M.orient(:,a)));
            %                     medias_unconst(a) = mu;
            %                     errores_unconst(a) = error;
            
            %                     v1 = l(1)*U(:,1);
            %                     v2 = l(2)*U(:,2);
            %                     v3 = l(3)*U(:,3);
            
            %                     plot3([0,v1(1)],[0,v1(2)],[0,v1(3)],'Color','b','LineWidth',4);
            %                     plot3([0,v2(1)],[0,v2(2)],[0,v2(3)],'Color','b');
            %                     plot3([0,v3(1)],[0,v3(2)],[0,v3(3)],'Color','b');
            
            %                     FA = sqrt(1/2)*sqrt(((l(1)-l(2))^2 + (l(1)-l(3))^2 + (l(2)-l(3))^2)/(l(1)^2 + l(2)^2 + l(3)^2));
            %                     disp(['Atom= ',num2str(a), ' Adapted Dict mu=',num2str(mu), ' fit error= ',num2str(100*error),'%','  FA=',num2str(FA)]);
            %                     disp(['Lambdas=',num2str(l')])
            
            %                      [ Qest, mu, error ] = Estimate_atom_tensor_demeaned( f, alpha, bvecs, bvals,'constrained');
            %                     %[U,S,V] = svd(Qest);
            %                     [U,l] = eig(Qest,'vector');
            %
            %                     [l,ind1] = sort(l,'descend');
            %                     U = U(:,ind1);
            %
            %                     ang_error_const(a) = acos(abs(U(:,1)'*M.orient(:,a)));
            %                     medias_const(a) = mu;
            %                     errores_const(a) = error;
            %
            %                     FA = sqrt(1/2)*sqrt(((l(1)-l(2))^2 + (l(1)-l(3))^2 + (l(2)-l(3))^2)/(l(1)^2 + l(2)^2 + l(3)^2));
            %                     disp(['(constrained) Atom= ',num2str(a), ' Adapted Dict mu=',num2str(mu), ' fit error= ',num2str(100*error),'%','  FA=',num2str(FA)]);
            %                     disp(['Lambdas=',num2str(l')])
            
            %%%%%%%%%%%%%%%%%%%%%%
            
            M.DictSig(:,a) = da;
            B(a,cols) = ba;
            
            % uptdate DB
            DB = DB + M.DictSig(:,a)*B(a,:);
            
            %error_post = norm(dSig-M.DictSig*B,'fro')/norm(dSig,'fro')
        end
        
        %disp(['atom ',num2str(a),' of ',num2str(nAtoms)]);
    end
    close(h)
    
    %disp(['Fit Dic and B, iter',num2str(iter),' error=',num2str(100*norm(dSig-M.DictSig*B,'fro')/norm(dSig,'fro'))])
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
tol = 1e-8;
atom = atom/norm(atom);

% E(isnan(E))=0;
% E(isinf(E))=0;
[u1,s1,v1] = svds(E,1);
% try positive sing vectors
da = u1;
ba = s1*v1;
ba(ba<0) = 0;
error = norm(E-da*ba','fro');
delta = Inf;
while delta > tol && sum(ba)
    ba = da'*E;
    ba(ba<0) = 0;
    if sum(ba) 
        da = E*ba'/(ba*ba');
    end
    errorn = norm(E-da*ba,'fro');
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
ba(ba<0) = 0;
error = norm(E-da*ba','fro');
delta = Inf;
while delta > tol && sum(ba) 
    ba = da'*E;
    ba(ba<0) = 0;
    if sum(ba) 
        da = E*ba'/(ba*ba');
    end
    errorn = norm(E-da*ba,'fro');
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

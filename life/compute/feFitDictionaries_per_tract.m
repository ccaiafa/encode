function [fe] = feFitDictionaries_per_tract(varargin)
fe = varargin{1};
lambda = varargin{2};
error_trhesholdSig = varargin{3};
nIterMax = varargin{4};

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
normSig = norm(dSig,'fro');

pSig = reshape(feGet(fe,'psigfiber'),[nTheta,nVoxels]);
Error = dSig - pSig;

%% Make loop updating each dictionary at a time
nTract = size(fe.life.M.tracts,2);

delta_error = Inf;
error_old = norm(Error,'fro')/normSig;
nIter = 1;
while (delta_error > error_trhesholdSig) && (nIter <= nIterMax)
    disp(['Iter=',num2str(nIter),' Error=',num2str(error_old), '  delta_error=', num2str(delta_error)])
    for n=1:nTract
        Phi_n = fe.life.M.Phi(:,:,fe.life.M.tracts{n}.ind);
        w_n = w(fe.life.M.tracts{n}.ind);
        D_n = fe.life.M.tracts{n}.DictSig;
        
        pSig_n = M_times_w(Phi_n.subs(:,1),Phi_n.subs(:,2),Phi_n.subs(:,3),Phi_n.vals, D_n, w_n, nTheta, nVoxels); % tract n prediction
        pSig_n = reshape(pSig_n,[nTheta,nVoxels]);
        dSig_n = Error + pSig_n;
        
        Bmat_n = ttv(Phi_n, w_n,3);
        [ind, val] = find(Bmat_n);
        Bmat_n = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
        
        disp(['Tract=',num2str(n)])
        [fe.life.M.tracts{n}.DictSig, Bmat_n] = update_DandB(dSig_n, D_n, Bmat_n, lambda);
        
        %% Compute Phi_n compatible with Bmat_n(a,v)
        [sub, ~] = find(Phi_n);
        ind = sub2ind(size(Bmat_n), sub(:,1), sub(:,2));
        b = Bmat_n(:);
        
        A = sptensor(sub, ones(size(sub,1),1), size(Phi_n));
        A = ttv(A,w_n.^2,3);
        [subA, val] = find(A);
        A = sparse(subA(:,1),subA(:,2),val,nAtoms,nVoxels);
        a = A(:);
        div = a(ind);
        [ind1,~,val] = find(div);
        
        newval = full(w_n(sub(:,3)).*b(ind));
        newval(ind1) = newval(ind1)./val;
        
        Phi_n = sptensor(sub, newval, size(Phi_n));
        % Update tensor Phi with the block of tract n, Phi_n
        fe.life.M.Phi(:,:,fe.life.M.tracts{n}.ind') = Phi_n;
        
        % Update Error by substracting the new tract_n prediction (new Phi_n and new D_n)
        Error = dSig_n - reshape(M_times_w(Phi_n.subs(:,1),Phi_n.subs(:,2),Phi_n.subs(:,3),Phi_n.vals, fe.life.M.tracts{n}.DictSig, w_n, nTheta, nVoxels),[nTheta,nVoxels]);
    end
    nIter = nIter + 1;
    delta_error = abs(norm(Error,'fro')/normSig - error_old);
    error_old = norm(Error,'fro')/normSig;
end

end

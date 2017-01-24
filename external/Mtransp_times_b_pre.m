function [ w ] = Mtransp_times_b_pre(A, dSig)
[nFibers] = size(A.Phi,3);
if ~isfield(A,'Dictionaries') % If there are not adaptive dictionaries available then use original dictionary
    w = Mtransp_times_b(A.Phi.subs(:,1),A.Phi.subs(:,2),A.Phi.subs(:,3),A.Phi.vals,A.DictSig,dSig,nFibers);
elseif isfield(A,'tracts') && isfield(A.tracts{1},'DictSig')
    % If there are dictionary associated to tracts, use them
    nTract = size(A.tracts,2);
    w = zeros(nFibers, 1);
    for n=1:nTract
        Phi_sub = A.Phi(:,:,A.tracts{n}.ind);
        w(A.tracts{n}.ind) = Mtransp_times_b(Phi_sub.subs(:,1),Phi_sub.subs(:,2),Phi_sub.subs(:,3),Phi_sub.vals,A.tracts{n}.DictSig,dSig,length(A.tracts{n}.ind));
    end
    %for n=1:nTract
    %    w(A.tracts{n}.ind) = w_aux{n};
    %end
    
else
    nDict = size(A.Dictionaries,2);
    w = zeros(nFibers, 1);
    Phi = A.Phi;
    for n=1:nDict
        Phi_sub = Phi(:,A.ind_vox{n},:);
        w = w + Mtransp_times_b(Phi_sub.subs(:,1),Phi_sub.subs(:,2),Phi_sub.subs(:,3),Phi_sub.vals,A.Dictionaries{n},dSig(:,A.ind_vox{n}),nFibers);
    end
end
end



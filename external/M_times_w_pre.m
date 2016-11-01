function b = M_times_w_pre(A,w)
[nTheta]  = size(A.DictSig,1);
[nVoxels] = size(A.Phi,2);
if ~isfield(A,'Dictionaries') % If there are not adaptive dictionaries available then use original dictionary
    b = M_times_w(A.Phi.subs(:,1),A.Phi.subs(:,2),A.Phi.subs(:,3),A.Phi.vals,A.DictSig,w,nTheta,nVoxels);
    b = reshape(b,[nTheta, nVoxels]);
else % It there are adaptive dictionaries available, use them
    nDict = size(A.Dictionaries,2);
    b = zeros(nTheta, nVoxels);
    Phi = A.Phi;
    for n=1:nDict
        Phi_sub = Phi(:,A.ind_vox{n},:);
        sub_val = M_times_w(Phi_sub.subs(:,1),Phi_sub.subs(:,2),Phi_sub.subs(:,3),Phi_sub.vals,A.Dictionaries{n},w,nTheta,length(A.ind_vox{n}));
        sub_val =  reshape(sub_val,[nTheta, length(A.ind_vox{n})]);
        b(:,A.ind_vox{n}) = sub_val;
    end
    b = b(:);
end


end
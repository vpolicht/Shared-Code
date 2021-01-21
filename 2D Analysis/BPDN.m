function [z,end_state] = BPDN(basis_A, target)
assert(size(basis_A,1) == size(target,1)); %ensure the basis is matched to the target
assert(size(target,2) == 1); %ensure the target is a column vector

opts = spgSetParms('verbosity',0,'iscomplex',0);
sigma = std(target,[],1)/3;
[z,~,~,info] = spg_bpdn(basis_A, target, sigma, opts);
end_state = info.stat;
z = z';
end
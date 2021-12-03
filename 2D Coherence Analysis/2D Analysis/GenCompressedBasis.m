function [basis_A] = GenCompressedBasis(F,T)
gen_cos_matrix = @(f,t) cos(repmat(f',[size(t,1) 1]).*repmat(t,[1 size(f,1)]));
gen_sin_matrix = @(f,t) sin(repmat(f',[size(t,1) 1]).*repmat(t,[1 size(f,1)]));
% gen_poly_matrix = @(t) cat(2,ones(size(t,1),1),t,t.^2);
% basis_A = cat(2,gen_poly_matrix(T-T(length(T)/2)),gen_cos_matrix(F,T),gen_sin_matrix(F,T));
basis_A = cat(2,gen_cos_matrix(F,T),gen_sin_matrix(F,T));

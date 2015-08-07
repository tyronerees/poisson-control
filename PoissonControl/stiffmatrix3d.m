% stiffmatrixbdy.m
%
% The stiffness matrix for the standard (Laplace) problem including the
% boundary

function [K] = stiffmatrix3d(h,cub)

N = 1/h;
NN = (N+1)^3;

K = sparse(NN,NN);

% Local Stiffness Matrix
Ka = (h/12)*toeplitz([4 0 -1 -0]);         
Kb = (h/12)*toeplitz([0 -1 -1 -1]);
Kk = [Ka Kb; Kb Ka];

% Global Stiffness Matrix
lc = length(cub(:,1));
ntriplets = lc*64;
I = zeros(ntriplets,1);
J = zeros(ntriplets,1);
X = zeros(ntriplets,1);
ntriplets = 1;

% Global Stiffness Matrix
k=1:lc;
for krow = 1:8                                                         
    for kcol = 1:8
        I(ntriplets:ntriplets+lc-1) = (cub(k,krow));
        J(ntriplets:ntriplets+lc-1) = (cub(k,kcol));
        X(ntriplets:ntriplets+lc-1) = Kk(krow,kcol);
        ntriplets = ntriplets+lc;
    end
end

K = sparse(I,J,X,NN,NN);
% for k=1:length(bound.cub(:,1))
%   K(bound.cub(k,1:8),bound.cub(k,1:8))=K(bound.cub(k,1:8),bound.cub(k,1:8))+Kk(1:8,1:8); 
% end
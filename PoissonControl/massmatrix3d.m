% massmatrix.m
%
% This is the code that generates the mass matrix for the standard problem

function [M] = massmatrix3d(h,cub)        

N = 1/h;
NN = (N+1)^3;


% Local Mass Matrix
Ma = (h^3/108)*toeplitz([4 2 1 2]);     
Mb = (h^3/216)*toeplitz([4 2 1 2]);
Mk = [Ma Mb; Mb Ma];

lc = length(cub(:,1));
ntriplets = lc*64;
I = zeros(ntriplets,1);
J = zeros(ntriplets,1);
X = zeros(ntriplets,1);
ntriplets = 1;

% Global Mass Matrix
k=1:lc;
for krow = 1:8                                                         
    for kcol = 1:8
        I(ntriplets:ntriplets+lc-1) = (cub(k,krow));
        J(ntriplets:ntriplets+lc-1) = (cub(k,kcol));
        X(ntriplets:ntriplets+lc-1) = Mk(krow,kcol);
        ntriplets = ntriplets+lc;
    end
end

M = sparse(I,J,X,NN,NN);

% %Global Mass Matrix
% for k=1:length(cub(:,1))
%       M(cub(k,1:8),cub(k,1:8))=(M(cub(k,1:8),cub(k,1:8))+Mk(1:8,1:8));  
% end

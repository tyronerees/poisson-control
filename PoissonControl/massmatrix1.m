% massmatrix.m
%
% This is the code that generates the mass matrix of size (N+1)^2 for the standard problem

function [M] = massmatrix1(h)        % removed the output M2
N=1/h;
NN = (N+1)^2;

% see Tim Davis - http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/

M = sparse(NN,NN);

[x y quad] = mesh_bdy(N);

% Local Mass Matrix
Mk = ((h^2)/36)*toeplitz([4 2 1 2]);

lq = length(quad(:,1));
ntriplets = lq*16;
I = zeros(ntriplets,1);
J = zeros(ntriplets,1);
X = zeros(ntriplets,1);
ntriplets = 1;

% Global Stiffness Matrix
k=1:lq;
for krow = 1:4                                                          
    for kcol = 1:4
        I(ntriplets:ntriplets+lq-1) = (quad(k,krow));
        J(ntriplets:ntriplets+lq-1) = (quad(k,kcol));
        X(ntriplets:ntriplets+lq-1) = Mk(krow,kcol);
        ntriplets = ntriplets+lq;
    end
end

M = sparse(I,J,X,NN,NN);


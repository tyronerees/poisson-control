% bdymassmatrix.m
%
% This is the code that generates the mass matrix of size (N+1)^2 for the boundary problem

function [M] = bdymassmatrix(h)        % removed the output M2
N=1/h;
NN = (N+1)^2;

% see Tim Davis - http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/

% M = sparse(NN,NN);

[x y quad bdy] = mesh_bdy(N);

% Local Mass Matrix
Mk = ((h)/6)*[2 1; 1 2];

lq = length(bdy(:,1));
ntriplets = lq*4;
I = zeros(ntriplets,1);
J = zeros(ntriplets,1);
X = zeros(ntriplets,1);
ntriplets = 1;

% Global Stiffness Matrix
k=1:lq;
for krow = 1:2                                                          
    for kcol = 1:2
        I(ntriplets:ntriplets+lq-1) = (bdy(k,krow));
        J(ntriplets:ntriplets+lq-1) = (bdy(k,kcol));
        X(ntriplets:ntriplets+lq-1) = Mk(krow,kcol);
        ntriplets = ntriplets+lq;
    end
end

M = sparse(I,J,X,NN,NN);


function [K] = stiffmatrixbdy3(h)
N=1/h;
NN = (N+1)^2;

% see Tim Davis - http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/

[x y quad] = mesh_bdy(N);

% Local Stiffness Matrix
Kk = (1/6)*toeplitz([4 -1 -2 -1]);

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
        X(ntriplets:ntriplets+lq-1) = Kk(krow,kcol);
        ntriplets = ntriplets+lq;
    end
end

K = sparse(I,J,X,NN,NN);
        

function [P]=prol(nelx,x,bc)

% nel = 8;                % number of elements in the rows
% 
% newel = nel/2;          % number of elements in the rows after multigrid
% 
% P = sparse((nel+1)^2,(newel+1)^2,(3*newel+1)^2);
% 
% nelx = 4; nely = 4;
% x = (0:0.25:1)'; y=x;

ecx=nelx/2; 


tx = x(2:nelx+1) - x(1:nelx);
bx = zeros(nelx,1);
bx(1:2:nelx-1) = x(3:2:nelx+1) - x(1:2:nelx-1);
bx(2:2:nelx)   = bx(1:2:nelx-1);
dx = tx./bx;

p = spalloc(2*ecx+1,ecx+1,3*ecx+1);
for j=2:ecx
   p(2*j-2,j) = dx(2*j-3);
   p(2*j-1,j) = 1.0;
   p(2*j  ,j) = dx(2*j);
end
j=1;
   p(2*j-1,1) = 1.0;           
   p(2*j  ,1) = dx(2*j);
j=ecx+1;
   p(2*j-2,j) = dx(2*j-3);
   p(2*j-1,j) = 1.0;
   
P = kron(p,p);

cd PoissonControl
switch bc
    case 'mixed'
        dirn = unique(robbdy(nelx));
        dirn2 = unique(robbdy(nelx/2));
    case 'dirichlet'
        [dir,neu] = robbdy(nelx);
        dirn = unique(dir); neun = unique(neu);
        dirn = unique([dirn;neun]);
        [dir2,neu2] = robbdy(nelx/2);
        dirn2 = unique(dir2); neun2 = unique(neu2);
        dirn2 = unique([dirn2;neun2]);
    case 'neumann'
        dirn = unique(neubdy(nelx));
        dirn2 = unique(neubdy(nelx/2));
    case 'neumannbc'
        dirn = [];
        dirn2 = [];
end
cd ..
P(dirn,:) = [];
P(:,dirn2) = [];

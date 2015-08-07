function [P,bdn2]=prol(nel,v,bdn)

% nel = 8;                % number of elements in the rows
% 
% newel = nel/2;          % number of elements in the rows after multigrid
% 
% P = sparse((nel+1)^2,(newel+1)^2,(3*newel+1)^2);
% 
% nelx = 4; nely = 4;
% x = (0:0.25:1)'; y=x;

ec=nel/2; 
%

t = v(2:nel+1) - v(1:nel);
b = zeros(nel,1);
b(1:2:nel-1) = v(3:2:nel+1) - v(1:2:nel-1);
b(2:2:nel)   = b(1:2:nel-1);
% keyboard
d = t./b;

p = spalloc(2*ec+1,ec+1,3*ec+1);
for j=2:ec
   p(2*j-2,j) = d(2*j-3);
   p(2*j-1,j) = 1.0;
   p(2*j  ,j) = d(2*j);
end
j=1;
   p(2*j-1,1) = 1.0;           
   p(2*j  ,1) = d(2*j);
j=ec+1;
   p(2*j-2,j) = d(2*j-3);
   p(2*j-1,j) = 1.0;
   
p1 = kron(p,p);    
P = sparse(kron(p,p1));

% keyboard

P(bdn,:) = [];
N=nel/2;
% k = 1:N;
% % bdy = [k' k'+1; k'*(N+1) (k'+1)*(N+1); (N+1)^2-k'+1 (N+1)^2-k'];
% k = N:-1:1;
% bdy = [bdy; k'*(N+1)+1 (k'-1)*(N+1)+1];
% bdn2 = unique(bdy);
% P(:,bdn2) = [];

% create bdy - a matrix where each row defines the nodes lying on the
% boundary for some side of an element 

M = [(1:N)', (2:N+1)', (N+3:2*N+2)', (N+2:2*N+1)'];
bdy3 = kron((N+1)*(0:N-1)',ones(N,4))+ kron(ones(N,1),M);     % bottom of cube
bdy3 = [bdy3;bdy3+N*(N+1)^2];                                 % top of cube
M = [(1:N)', (2:N+1)', (((N+1)^2+2):((N+1)^2+N+1))', (((N+1)^2+1):((N+1)^2+N))'];
bdy = kron(ones(N,1),M)+kron((N+1)^2*(0:N-1)',ones(N,4));     % front of cube
bdy3 = [bdy3; bdy];   
bdy3 = [bdy3;bdy+N*(N+1)]   ;                                 % add back of cube (maybe need to change ordering?!)
M = [(1:N+1:N*(N+1))', (N+2:N+1:(N+1)^2-N)', ((N+1)*(N+2)+1:N+1:(2*(N+1)^2-N))', (((N+1)^2+1):N+1:((N+1)*(2*N+1)))'];
bdy = kron(ones(N,1),M)+kron((N+1)^2*(0:N-1)',ones(N,4));
bdy3 = [bdy3; bdy];                                           % add the left side of the cube (again, may need to reorder?)
bdy3 = [bdy3; bdy+N];                                         % and the right...(ditto)

bdn2 = unique(bdy3);

P(:,bdn2) = [];



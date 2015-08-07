% mesh_bdy.m
%
% Generates the mesh for Q1 approximation
%
% Input  - size.N is the number of elements per edge
%          size.h is the length of an element
%
% Output - bound.x, bound.y and bound.z are vectors so that node i is 
%             [x(i) y(i) z(i)]
%          bound.cub is the connectivity matrix  - the eight numbers in 
%             each row define the nodes that make up an element
%          bound.bdy is the matrix where each row defines the nodes 
%             lying on the boundary for some side of an element 
%          bound.bdn is the vector of nodes lying on the boundary
%          bound.h is the length of an element

function [x y z cub bdy] = mesh3d_bdy(N)

h = 1/N;

% Form x,y and z so that node i is [x(i) y(i) z(i)], etc. so we go along x
% first, then y, and finally z

vec = 0:h:1;
x1 = kron(ones(1,N+1),vec);
y1 = kron(h*(0:N),ones(1,N+1));
x = kron(ones(1,N+1),x1);
y = kron(ones(1,N+1),y1);
z = kron(h*(0:N),ones(1,(N+1)^2));

% Form the connectivity matrix cub - the eight numbers in each row define
% the nodes that make up an element

% quad=[];
% for l = 0:N-1
%     for i=0:N-1
%         for k = 1:N
%             quad = [quad;...
%             l*(N+1)^2+i*(N+1)+k, l*(N+1)^2+i*(N+1)+k+1,...
%             l*(N+1)^2+(i+1)*(N+1)+k+1, l*(N+1)^2+(i+1)*(N+1)+k,...
%             (l+1)*(N+1)^2+i*(N+1)+k, (l+1)*(N+1)^2+i*(N+1)+k+1,...
%             (l+1)*(N+1)^2+(i+1)*(N+1)+k+1, (l+1)*(N+1)^2+(i+1)*(N+1)+k];
%         end
%     end 
% end
% M = [(1:N)', (2:N+1)', (N+3:2*N+2)', (N+2:2*N+1)'];
% E = ones(N,4);
% quad = kron( (N+1)*(0:N-1)' ,E) + kron(ones(N,1),M );

M = [(1:N)', (2:N+1)', (N+3:2*N+2)', (N+2:2*N+1)',(1:N)', (2:N+1)', (N+3:2*N+2)', (N+2:2*N+1)'];
cub2 = kron((N+1)*(0:N-1)' ,ones(N,8)) + kron(ones(N,1),M )+[zeros(N^2,4),(N+1)^2*ones(N^2,4)];
cub = kron(ones(N,1),cub2) + kron((0:N-1)',(N+1)^2*ones(N^2,8));

% create bdy3 - a matrix where each row defines the nodes lying on the
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

bdy = bdy3;

bdn = unique(bdy3);

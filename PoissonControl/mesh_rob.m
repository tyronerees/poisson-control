% mesh_bdy.m
%
% Generates the mesh for Q1 approximation

function [x y quad dir neu] = mesh_rob(N)

h = 1/N;                             % define the mesh size

z = 0:h:1;
x = kron(ones(1,N+1),z);             % define the x vector
y = kron(h*(0:N),ones(1,N+1));       % define the y vector


M = [(1:N)', (2:N+1)', (N+3:2*N+2)', (N+2:2*N+1)'];
E = ones(N,4);
quad = kron( (N+1)*(0:N-1)' ,E) + kron(ones(N,1),M );   % define quad, which gives you the node
                                                        % which gives you the node numbers on each element 
%norm(quad-quado)                                       % on each element
[dir,neu] = robbdy(N);

% keyboard
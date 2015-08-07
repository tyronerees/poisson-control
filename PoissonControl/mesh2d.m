% 2dmesh.m
%
% Generates the mesh for Q1 approximation in 2D

function [x y quad] = mesh2d(N)

h = 1/N;                             % define the mesh size

z = 0:h:1;
x = kron(ones(1,N+1),z);             % define the x vector
y = kron(h*(0:N),ones(1,N+1));       % define the y vector


M = [(1:N)', (2:N+1)', (N+3:2*N+2)', (N+2:2*N+1)'];
E = ones(N,4);
quad = kron( (N+1)*(0:N-1)' ,E) + kron(ones(N,1),M );
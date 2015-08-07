% mesh_bdy.m
%
% Generates the mesh for Q1 approximation

function [x y quad bdy] = mesh_bdy(N)

h = 1/N;                             % define the mesh size

z = 0:h:1;
x = kron(ones(1,N+1),z);             % define the x vector
y = kron(h*(0:N),ones(1,N+1));       % define the y vector


%for i = 1:N
%    x = [x 0:h:1];
%end

%y = zeros(1,N+1);
%for j = 1:N
%    y1=y(end)+h;
%    y = [y y1*ones(1,N+1)];
%end

%quado=[];
%for i=0:N-1
%    for k = 1:N
%        quado = [quado; i*(N+1)+k i*(N+1)+k+1 (i+1)*(N+1)+k+1 (i+1)*(N+1)+k];
%    end 
%end

M = [(1:N)', (2:N+1)', (N+3:2*N+2)', (N+2:2*N+1)'];
E = ones(N,4);
quad = kron( (N+1)*(0:N-1)' ,E) + kron(ones(N,1),M );   % define quad, which gives you the node
                                                        % which gives you the node numbers on each element 
%norm(quad-quado)                                       % on each element

k = 1:N;
bdy = [k' k'+1; k'*(N+1) (k'+1)*(N+1); (N+1)^2-k'+1 (N+1)^2-k'];
k = N:-1:1;
bdy = [bdy; k'*(N+1)+1 (k'-1)*(N+1)+1];                 % define bdy, which gives you the node numbers
                                                        % touching the boundary

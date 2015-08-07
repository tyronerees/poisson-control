
% pdecoDriver.m
%
% Generates the saddle point matrix, the stiffness matrix and mass matrix
% for the pde test problem 1. 
% 
% [A,B,K,M,dirn,ubdy] = PDEDriverrob(N,beta,bc,ob)
% 
% INPUT - blank returns default value
% N -    number of elements on one axis (1/h, where h is the mesh size)
% beta - Regularization parameter
% bc -   What boundary conditions to enforce 'dirichlet', 'neumann' or
%        'mixed'
% ob - the objective function - 1,2,3...
%
% OUTPUT
% A - saddle point matrix
% b - right hand side
% K - stiffness matrix
% M - mass matrix/home/rees/Thesis/Thesis code/PoissonControl

function [A,b,bdy_set,ubdy,uhat,def_setup,prob_setup] = pdecoDriver_3d(def_setup)%(N,beta,bc,ob)

%% set N and h
N=2^def_setup.pow;
h=1/N;

%%



switch def_setup.bc
    case 'dirichlet'
        [x y z cub bdy] = mesh3d_bdy(N);
        bdy_set.dirn = unique(bdy');
        bdy_set.neun = [];
        bdy_set.bdn = bdy_set.dirn;
        bdy_set.dbynodes = 1:6*N^2+2; % set nodes on boundary with dirichlet bcs
    case 'neumann'
        fprintf('Neumann bcs for 3d problems not supported\n')
%         return
%         [dir,neu] = neubdy3d(N);
%         bdy_set.dirn = unique(dir);
%         bdy_set.neun = unique(neu);
%         bdy_set.bdn = unique([bdy_set.dirn;bdy_set.neun]);
%         bdy_set.dbynodes = [4*N]; % set nodes on boundary with dirichlet bcs
    case 'mixed'
        fprintf('Mixed bcs for 3d problems not supported\n')
        return
%         [dir,neu] = robbdy3d(N);
%         bdy_set.dirn = unique(dir);
%         bdy_set.neun = unique(neu);
%         bdy_set.bdn = unique([bdy_set.dirn;bdy_set.neun]);
%         bdy_set.dbynodes = [1:N+1,N+2:2:3*N];
end
[b K M ubdy uhat]=setupmat_3d(h,bdy_set,def_setup.ob);
prob_setup.yelt = 'q1';
prob_setup.uelt = 'q1';
prob_setup.dim = 3;
prob_setup.eqn = 'poiss';
prob_setup.bc = 'poiss';

l=length(M);

prob_setup.nu = length(M); % size of control
prob_setup.ny = length(K); % size of state

A = [def_setup.beta*M sparse(l,l) -M;...                           
    sparse(l,l) M K';...           
    -M K sparse(l,l)];                  



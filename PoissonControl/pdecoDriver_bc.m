
% pdecoDriver_bc.m
%
% Generates the saddle point matrix, the stiffness matrix and mass matrix
% for the pde test problem with boundary control. 
% 
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

function [A,b,bdy_set,ubdy,uhat,def_setup,prob_setup] = pdecoDriver_bc(def_setup)%(N,beta,bc,ob)

%% set N and h
N=2^def_setup.pow;
h=1/N;

%%


% [dir,neu] = neubdy(N);
% bdy_set.dirn = unique(dir);
% bdy_set.neun = unique(neu);
% bdy_set.bdn = unique([bdy_set.dirn;bdy_set.neun]);
% bdy_set.dbynodes = [4*N]; % set nodes on boundary with dirichlet bcs

[dir,neu] = neubdy(N);
dirn = unique(dir);
neun = unique(neu);
bdy_set.bdn = unique([dirn;neun]);
bdy_set.dirn = [];
bdy_set.neun = bdy_set.bdn;
bdy_set.dbynodes = []; % set nodes on boundary with dirichlet bcs

[b K Mu My Mhat ubdy uhat]=setupmat_bound(h,bdy_set,def_setup.ob,def_setup.plots);
prob_setup.yelt = 'q1';
prob_setup.uelt = 'q1';
prob_setup.dim = 2;
prob_setup.eqn = 'poiss';
prob_setup.bc = 'poiss';

def_setup.bc = 'neumannbc';

prob_setup.nu = length(Mu); % size of control
prob_setup.ny = length(K); % size of state

A = [def_setup.beta*Mu sparse(prob_setup.nu,prob_setup.ny) -Mhat';...                           
    sparse(prob_setup.ny,prob_setup.nu) My K';...           
    -Mhat K sparse(prob_setup.ny,prob_setup.ny)];                  



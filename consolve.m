% poissconsolve.m
%
% Solves the poisson optimal control problem
% 
% W = poissconsolve(A,b,def_setup)
% 
% INPUT - blank returns default value
% A -    saddle point matrix
% b -    Right hand size
% def_setup - Default parameters from matrix setup
% 
%
% OUTPUT
% W - solution

function [W,iter,t2] = consolve(A,b,prob_setup,def_soln,def_setup,bdn_set)%(N,beta,bc,ob)


if strcmp(def_soln.method,'backslash') == 1
    fprintf('Solving with backslash...\n')
    tic,
    W = A\b;
    t2 = toc;
    iter = NaN;
elseif strcmp(def_soln.method,'minres') == 1
    fprintf('Solving with MINRES...\n')
%     mmethod = 'chebit';
%     kmethod = 'amg';
    [W,iter,t2] = pdeminres(A,b,prob_setup,def_soln,def_setup);
    fprintf('Solving mass matrices with %s...\n ',def_soln.mmethod);
    fprintf('Solving stiffness matrices with %s...\n',def_soln.kmethod);
    fprintf('Converged in %d iterations and %d s \n',iter, t2);
elseif strcmp(def_soln.method,'gmres') == 1
    fprintf('Solving with GMRES...\n')
%     mmethod = 'chebit';
%     kmethod = 'amg';
    [W,iter,t2] = pdegmres(A,b,prob_setup,def_soln,def_setup);
    fprintf('Solving mass matrices with %s...\n ',def_soln.mmethod);
    fprintf('Solving stiffness matrices with %s...\n',def_soln.kmethod);
    fprintf('Converged in %d iterations and %d s \n',iter, t2);
elseif strcmp(def_soln.method,'bpcg') == 1
    fprintf('Solving with BPCG...\n')
    [W,iter,t2] = pdebpcg(A,b,prob_setup,def_soln,def_setup);
    fprintf('Solving stiffness matrices with %s...\n',def_soln.kmethod);
    fprintf('Converged in %d iterations and %d s \n',iter, t2);
elseif strcmp(def_soln.method,'ppcg') == 1
    fprintf('Solving with PPCG...\n')
%     mmethod = 'chebit';
%     kmethod = 'amg';
    [W,iter,t2] = pdeppcg(A,b,prob_setup,def_soln,def_setup);
    fprintf('Solving mass matrices with %s...\n ',def_soln.mmethod);
    fprintf('Solving stiffness matrices with %s...\n',def_soln.kmethod);
    fprintf('Converged in %d iterations and %d s \n',iter, t2);
else
    error('Invalid solution method')
end


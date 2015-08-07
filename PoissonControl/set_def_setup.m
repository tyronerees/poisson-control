function def_setup = set_def_setup

def_setup.pow = 5;          % No of elts on edge is 2^pow
def_setup.bc = 'dirichlet'; % dirichlet, neumann or mixed
def_setup.beta = 1e-2;      % regularization parameter
def_setup.ob = 1;           % objective function -- 1,2 or 3
def_setup.plots = 1;        % 0 no plots,1 for plots with subfigures or 2 for no subfigures
def_setup.type = 'dist2d';  % dist2d, dist3d, bound2d


function def_soln = set_def_soln(prob_setup)


def_soln.ucheb = 20;  % number of iterations for the control mass matrix
def_soln.ycheb = 20;  % number of iterations for the state mass matrix
def_soln.vcyc = 2;    % number of v-cycles
def_soln.method = 'minres';    % method - minres, bpcg or ppcg  
def_soln.mmethod = 'chebit';   % method for approx mass matrix
def_soln.kmethod = 'amg';      % method for approx stiff matrix
def_soln.tol = 1e-6;           % tolerance for iterative methods
def_soln.conmethod = 'pblktri'; % method for constraint preconditioner
                                % 'pblktri' for psycologically block
                                % triangular, 'diag' for diagonal
def_soln.scale = 0.9;
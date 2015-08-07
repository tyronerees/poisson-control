function def_soln = set_def_soln(prob_setup)

switch lower(prob_setup.eqn)
    case 'poiss'
        def_soln.ucheb = 20;  % number of iterations for the control mass matrix
        def_soln.ycheb = 20;  % number of iterations for the state mass matrix
        def_soln.vcyc = 2;    % number of v-cycles
        def_soln.method = 'minres';    % method - minres, bpcg or ppcg  
        def_soln.mmethod = 'chebit';   % method for approx mass matrix
        def_soln.kmethod = 'amg';      % method for approx stiff matrix
        def_soln.tol = 1e-6;  % tolerance for method
        def_soln.conmethod = 'pblktri'; % method for constraint preconditioner
        def_soln.resvec = 0;
                                % 'pblktri' for psycologically block
                                % triangular, 'diag' for diagonal
        def_soln.chebit_ppcg = 30;   % number of iterations for ppcg mass approx  
        def_soln.twonorm = 0;   % method for testing convergence
        def_soln.gmgpre = 2;   % number of pre-smoothing steps for gmg
        def_soln.gmgpost = 2;   % number of post-smoothing steps for gmg
        def_soln.loadprol = 1;  % 0 to generate multdata for gmg, 1 to load
        def_soln.scale = 0.9;
    case 'stokes'
        
    case 'advdiff'
        
end
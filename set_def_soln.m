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
                                % 'pblktri' for psycologically block
                                % triangular, 'diag' for diagonal
        def_soln.chebit_ppcg = 30;   % number of iterations for ppcg mass approx  
        def_soln.twonorm = 0;   % method for testing convergence
        def_soln.gmgpre = 2;   % number of pre-smoothing steps for gmg
        def_soln.gmgpost = 2;   % number of post-smoothing steps for gmg
        def_soln.loadprol = 1;  % 0 to generate multdata for gmg, 1 to load
    case 'stokes'
        def_soln.velcheb = 20;  % number of iterations for the control mass matrix
        def_soln.pcheb = 20;
        def_soln.sprecheb = 20;
        def_soln.ycheb = 20;  % number of iterations for the state mass matrix
        def_soln.vcyc = 2;    % number of v-cycles
        def_soln.method = 'minres';    % method - minres, bpcg or ppcg  
        def_soln.mmethod = 'chebit';   % method for approx mass matrix
        def_soln.kmethod = 'ltri';      % method for approx stiff matrix
        def_soln.tol = 1e-6;  % tolerance for method
        def_soln.conmethod = 'pblktri'; % method for constraint preconditioner
                                % 'pblktri' for psycologically block
                                % triangular, 'diag' for diagonal
        def_soln.chebit_ppcg = 30;   % number of iterations for ppcg mass approx  
        def_soln.twonorm = 0;   % method for testing convergence
        def_soln.gmgpre = 2;   % number of pre-smoothing steps for gmg
        def_soln.gmgpost = 2;   % number of post-smoothing steps for gmg
        def_soln.loadprol = 1;  % 0 to generate multdata for gmg, 1 to load
        def_soln.alp = 3/5; % Scaling parameter for ine. uzawa
        def_soln.Uits = 2; % No of steps of ine. Uzawa
        def_soln.kmethod_stokes = 'bslash';
        def_soln.mmethod_stokes = 'chebit';
        def_soln.scale = 0.9;
        def_soln.exactschur = 0;  % 0 for Mp approx to Schur complement

        
        
    case 'advdiff'
        def_soln.ucheb = 20;  % number of iterations for the control mass matrix
        def_soln.ycheb = 20;  % number of iterations for the state mass matrix
        def_soln.vcyc = 2;    % number of v-cycles
        def_soln.method = 'backslash';    % method - minres, bpcg or ppcg  
        def_soln.mmethod = 'chebit';   % method for approx mass matrix
        def_soln.kmethod = 'amg';      % method for approx stiff matrix
        def_soln.tol = 1e-6;  % tolerance for method
        def_soln.conmethod = 'pblktri'; % method for constraint preconditioner
                                % 'pblktri' for psycologically block
                                % triangular, 'diag' for diagonal
        def_soln.chebit_ppcg = 30;   % number of iterations for ppcg mass approx  
        def_soln.twonorm = 0;   % method for testing convergence
        def_soln.gmgpre = 2;   % number of pre-smoothing steps for gmg
        def_soln.gmgpost = 2;   % number of post-smoothing steps for gmg
        def_soln.loadprol = 1;  % 0 to generate multdata for gmg, 1 to load
end
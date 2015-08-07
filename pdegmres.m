function [x, iit,t2] = pdegmres(A,b,prob_setup,def_soln,def_setup)

%  -- Iterative template routine --
%     Based on code by 
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%     Modified by Tyrone Rees 
%
% [x, error, iter, flag] = gmres( A, x, b, M, restrt, max_it, tol )
%
% gmres.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector (TR - hard codeded)
%         b        REAL right hand side vector
%         P        REAL preconditioner matrix (TR - hard codeded)
%         restrt   INTEGER number of iterations between restarts (TR - hard codeded)
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         iit      INTEGER number of the inner iterations performed on the
%                  last restart
%         t2       REAL time taken to solve the system
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TR 26/05/2009
tic,
	
tol = def_soln.tol;
nu = prob_setup.nu;
ny = prob_setup.ny;
lb = 2*ny+nu;

x = sparse(lb,1); % set intial guess = 0
restrt = 500;             % set the restart parameter;
max_it =200;

% get the preconditioner

% Write A = [bM1 0 M3;
%            0 M2 -K1;
%           M4 -K2 0];

% Mu = A(1:nu,1:nu);
% My = A(nu+1:nu+ny,nu+1:nu+ny);
% Muy = A(nu+ny+1:2*ny+nu,1:nu);
% K = A(nu+ny+1:2*ny+nu,nu+1:nu+ny);



% P = blkdiag(M1, M2, (M4*(M1\M3) + K2*(M2\K1) ) );
% P = blkdiag(M1, M2, K2*(M2\K1) );
%
switch lower(def_soln.kmethod)
    case 'gmg'
        if def_soln.loadprol == 0
            multdata = genmultdata(def_setup.pow,def_setup.bc);
        elseif def_soln.loadprol == 1
            switch def_setup.bc
                case 'dirichlet'
                    load multdataQ1_dir
                case 'neumann'
                    load multdataQ1_neu
                case 'mixed'
                    load multdataQ1_mix
            end
        end
        multdata(1).pow = def_setup.pow;
        multdata(1).dim = prob_setup.dim;
        multdata(1).spre = def_soln.gmgpre;
        multdata(1).spost = def_soln.gmgpost;
        if strcmp(prob_setup.eqn,'advdiff')==1
            % if advection diffusion, use galerkin rest. and gauss seidel
            multdata(1).smoother = 'gs';
            multdata = genA(def_setup,multdata,A(nu+ny+1:2*ny+nu,nu+1:nu+ny));
            switch def_setup.meth
                case 'dto'
                    multdata = genAt(def_setup,multdata,A(nu+1:nu+ny,nu+ny+1:2*ny+nu));
                case 'otd'
                    def_setup2 = def_setup;
                    def_setup2.c = -def_setup.c;
                    def_setup2.s = -def_setup.s;
                    multdata = genAadj(def_setup2,multdata,A(nu+1:nu+ny,nu+ny+1:2*ny+nu));
            end
        end
    case 'amg'
        control = mi20_control;
        control.v_iterations = def_soln.vcyc;
        inform = mi20_setup(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),control);
end
if strcmp(def_soln.kmethod,'gmg') == 0
    multdata = [];
end
% switch lower(method)
%     case 'multigrid'
        %         load multdataQ1_nbdy
        %         smoother = 'gs';
        %         Amat2 = genA(pow,supg,multdata,smoother,K2);
%         %         supgadj = supg;
%         %         keyboard
%         supga = supg;
%         supga.c = -supg.c;
%         supga.s = -supg.s;
%         Amat1 = genA(pow,supga,multdata,smoother,K1);
%         %         Amat1 = Amat2;
%         %         for i = 1:pow+1
%         %             Amat1(i).A = Amat1(i).A';
%         %         end
%     case 'backslash'
%         multdata = 1;  % give a phoney value to multdata
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter = 0;                                         % initialization
flag = 0;


bnrm2 = norm( b );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

% r = P \ ( b-A*x );
% 
% Change the above line to make it problem-specific
rres = ( b-A*x );
kapp = 1;
r = presolve(A,rres,prob_setup,def_soln,multdata);

error = norm( r ) / bnrm2;
if ( error < tol ), return, end
[n,n] = size(A);                                  % initialize workspace
m = restrt;
V(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(n,1);
e1(1) = 1.0;

for iter = 1:max_it,                              % begin iteration

    %     r = P \ ( b-A*x );
    %
    % Change the above line to make it problem-specific
    rres = ( b-A*x );
    
    r = presolve(A,rres,prob_setup,def_soln,multdata);

    V(:,1) = r / norm( r );
    s = norm( r )*e1;
    iit = 0;
    for i = 1:m,                                   % construct orthonormal
        %         w = P \ (A*V(:,i));                         % basis using Gram-Schmidt
        %
        % Change the above line to make it problem-specific
        AV = A*V(:,i);
        
        w = presolve(A,AV,prob_setup,def_soln,multdata);
        %
        for k = 1:i,
            H(k,i)= w'*V(:,k);
            w = w - H(k,i)*V(:,k);
        end
        H(i+1,i) = norm( w );
        V(:,i+1) = w / H(i+1,i);
        for k = 1:i-1,                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = cs(i)*s(i);                        % approximate residual norm
        s(i+1) = -sn(i)*s(i);
        s(i)   = temp;
        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
        error  = abs(s(i+1)) / bnrm2;
        iit = iit+1;
        if ( error <= tol ),                        % update approximation
            y = H(1:i,1:i) \ s(1:i);                 % and exit
            x = x + V(:,1:i)*y;
            break;
        end
        
    end

    if ( error <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    x = x + V(:,1:m)*y;                            % update approximation
    %     r = P \ ( b-A*x );                              % compute residual %EEK!! CHECK THIS
    %
    % Change the above line to make it problem-specific
    rres = ( b-A*x );
    
    r = presolve(A,rres,prob_setup,def_soln,multdata);
    
    s(i+1) = norm(r);
    error = s(i+1) / bnrm2;                        % check convergence
    if ( error <= tol ), break, end;
    iter
end

if ( error > tol ), flag = 1; end;                 % converged

t2 = toc;

switch lower(def_soln.kmethod)
    case 'amg'
        inform = mi20_finalize;
    case 'amgsch'
        inform = mi20_finalize;
end 

%%%%%%%%%%%%%%%%%%%%%%%%
function z = presolve(A,v,prob_setup,def_soln,multdata)
nu = prob_setup.nu;
ny = prob_setup.ny;
lb = 2*ny+nu;

switch lower(def_soln.mmethod)
    case 'chebit'
        zed1 = chebsemiit(A(1:nu,1:nu),v(1:nu),def_soln.ucheb,prob_setup.dim,prob_setup.yelt);
        zed2 = chebsemiit(A(nu+1:nu+ny,nu+1:nu+ny),v(nu+1:nu+ny),def_soln.ycheb,prob_setup.dim,prob_setup.uelt);
    case 'bslash'
        zed1 = A(1:nu,1:nu)\v(1:nu);
        zed2 =A(nu+1:nu+ny,nu+1:nu+ny)\v(1+nu:nu+ny); 
end

switch lower(def_soln.kmethod)
    case 'gmg'
        if strcmp(prob_setup.eqn,'advdiff')==1 
            switch def_soln.schapp
                case 'mass';
                    
                case 'stiff';
                    zed3a = sparse(ny,1);
                    for i = 1:def_soln.vcyc
                        zed3a = mgvcyc_sd1(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),v(nu+ny+1:end),multdata,zed3a,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
                    end
                    zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
                    zed3 = sparse(ny,1);
                    for i = 1:def_soln.vcyc
                        zed3 = mgvcyc_sd2(A(nu+1:nu+ny,nu+ny+1:2*ny+nu),zed3b,multdata,zed3,multdata(1).pow,multdata(1).spost,multdata(1).spre,multdata(1).dim);
                    end
            end
        else
            zed3a = sparse(ny,1);
            for i = 1:def_soln.vcyc
                zed3a = mgvcyc(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),v(nu+ny+1:end),multdata,zed3a,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
            end
            zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
            zed3 = sparse(ny,1);
            for i = 1:def_soln.vcyc
                zed3 = mgvcyc(A(nu+1:nu+ny,nu+ny+1:2*ny+nu),zed3b,multdata,zed3,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
            end
        end
    case 'amg'
        zed3a = mi20_precondition(v(nu+ny+1:end));
        zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
        zed3 = mi20_precondition(zed3b);
    case 'bslash'
        switch def_soln.schapp
            case 'mass';
                zed3 = A(1:nu,nu+ny+1:2*ny+nu)\(A(1:nu,1:nu)*(A(nu+ny+1:2*ny+nu,1:nu)\v(nu+ny+1:end)));
            case 'stiff';
                zed3 = A(nu+1:nu+ny,nu+ny+1:2*ny+nu)\(A(nu+1:nu+ny,nu+1:nu+ny)*(A(nu+ny+1:2*ny+nu,nu+1:nu+ny)\v(nu+ny+1:end)));
        end
end

z = [zed1; zed2; zed3]; 

%%
% END of gmres.m
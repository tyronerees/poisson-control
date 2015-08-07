function [x, error, iter,iit,t2, flag] = gmressd( A, b, max_it, tol,pow,method,supg )

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



x = sparse(length(A),1); % set intial guess = 0
restrt = 500;             % set the restart parameter;

% get the preconditioner

% Write A = [bM1 0 M3;
%            0 M2 -K1;
%           M4 -K2 0];

l = length(A)/3;

M1 = A(1:l,1:l);
M2 = A(l+1:2*l,l+1:2*l);

K1 = -A(l+1:2*l,2*l+1:3*l);
K2 = -A(2*l+1:3*l,l+1:2*l);

M3 = A(1:l,2*l+1:3*l);
M4 = A(2*l+1:3*l,1:l);

% P = blkdiag(M1, M2, (M4*(M1\M3) + K2*(M2\K1) ) );
% P = blkdiag(M1, M2, K2*(M2\K1) );
% 
switch lower(method)
    case 'multigrid'
        load multdataQ1_nbdy
        smoother = 'gs';
        Amat2 = genA(pow,supg,multdata,smoother,K2);
        %         supgadj = supg;
        %         keyboard
        supga = supg;
        supga.c = -supg.c;
        supga.s = -supg.s;
        Amat1 = genA(pow,supga,multdata,smoother,K1);
        %         Amat1 = Amat2;
        %         for i = 1:pow+1
        %             Amat1(i).A = Amat1(i).A';
        %         end
    case 'backslash'
        multdata = 1;  % give a phoney value to multdata
end

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
r = Psolve(M1,M2,M3,M4,K1,K2,rres,method,kapp,l,multdata,pow,supg,Amat1,Amat2,smoother);
% switch method
%     case 'backslash'
%         r1 = (M1\rres(1:l));
%         r2 = (M2\rres(l+1:l+l));
%         kapp = 1;
%         switch kapp
%             case 1
%                 r3 = K1\(M2*(K2\rres(l+l+1:end)));
%             case 0
%                 %         r3 = M3\(M1*(M4\rres(l+l+1:end)));
%                 r3 = M4\rres(l+l+1:end);
%         end
%     case 'multigrid'
%         r1 = chebit2(M1,rres(1:l),20);%./(2*beta);
%         r2 = chebit2(M2,rres(l+1:l+l),20);
% %         zed1 = M1\v1(1:l);%./(2*beta);
% %         zed2 = M2\v1(1+l:l+l);
%         r3a = sparse(l,1); 
% %         [zed3a,i2,resvec2] = sdmgsolve(K2,v1(l+l+1:end),multdata,zed3a,pow,supg,Amat2,smoother);
%         r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
%         r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
%         r3b = M2*r3a;
%         r3 = sparse(l,1);
% %         [zed3,i1,resvec1] = sdmgsolve(K1,zed3b,multdata,zed3,pow,supg,Amat1,smoother);
%         r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
%         r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
% % keyboard
%         
% end
% r = [r1;r2;r3];
%

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
    r = Psolve(M1,M2,M3,M4,K1,K2,rres,method,kapp,l,multdata,pow,supg,Amat1,Amat2,smoother);
%     switch method
%         case 'backslash'
% 
%             r1 = (M1\rres(1:l));
%             r2 = (M2\rres(l+1:l+l));
%             switch kapp
%                 case 1
%                     r3 = K1\(M2*(K2\rres(l+l+1:end)));
%                 case 0
%                     %             r3 = M3\(M1*(M4\rres(l+l+1:end)));
%                     r3 = M4\rres(l+l+1:end);
%             end
%         case 'multigrid'
%             %         zed1 = pcg4(M2,v2(1:m))./(2*beta);
% %         zed2 = pcg4(M,v2(m+1:m+l));
% %         r1 = chebit2(M1,rres(1:l),20);%./(2*beta);
% %         r2 = chebit2(M2,rres(l+1:l+l),20);
%         r1 = M1\v2(1:l)./(2*beta);
%         r2 = M2\v2(l+1:2*l);
%         
%         r3a = sparse(l,1);
%         r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
%         r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
%         r3b = M2*r3a;
%         r3 = sparse(l,1);
%         r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
%         r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
%     end
%     r = [r1;r2;r3];
%     % 

    V(:,1) = r / norm( r );
    s = norm( r )*e1;
    iit = 0;
    for i = 1:m,                                   % construct orthonormal
        %         w = P \ (A*V(:,i));                         % basis using Gram-Schmidt
        %
        % Change the above line to make it problem-specific
        AV = A*V(:,i);
        
        w = Psolve(M1,M2,M3,M4,K1,K2,AV,method,kapp,l,multdata,pow,supg,Amat1,Amat2,smoother);
%         w1 = (M1\AV(1:l));
%         w2 = (M2\AV(l+1:l+l));
%         switch kapp
%             case 1
%                 w3 = K1\(M2*(K2\AV(l+l+1:end)));
%             case 0
% %                 w3 = M3\(M1*(M4\AV(l+l+1:end)));
%                 w3 = M4\AV(l+l+1:end);
%         end
%         w = [w1;w2;w3];
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
    r = Psolve(M1,M2,M3,M4,K1,K2,rres,method,kapp,l,multdata,pow,supg,Amat1,Amat2,smoother);
    %     r1 = (M1\rres(1:l));
    %     r2 = (M2\rres(l+1:l+l));
    %     switch kapp
    %         case 1
    %             r3 = K1\(M2*(K2\rres(l+l+1:end)));
    %         case 0
    % %             r3 = M3\(M1*(M4\rres(l+l+1:end)));
    %             r3 = M4\rres(l+l+1:end);
    %     end
%     switch method
%         case 'backslash'
% 
%             r1 = (M1\rres(1:l));
%             r2 = (M2\rres(l+1:l+l));
%             switch kapp
%                 case 1
%                     r3 = K1\(M2*(K2\rres(l+l+1:end)));
%                 case 0
%                     %             r3 = M3\(M1*(M4\rres(l+l+1:end)));
%                     r3 = M4\rres(l+l+1:end);
%             end
%         case 'multigrid'
%             %         zed1 = pcg4(M2,v2(1:m))./(2*beta);
%             %         zed2 = pcg4(M,v2(m+1:m+l));
%             r1 = chebit2(M1,rres(1:l),20);%./(2*beta);
%             r2 = chebit2(M2,rres(l+1:l+l),20);
%             %         zed1 = M1\v2(1:l);%./(2*beta);
%             %         zed2 = M2\v2(l+1:2*l);
% 
%             r3a = sparse(l,1);
%             r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
%             r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
%             r3b = M2*r3a;
%             r3 = sparse(l,1);
%             r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
%             r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
%     end
%     r = [r1;r2;r3];
    %
    s(i+1) = norm(r);
    error = s(i+1) / bnrm2;                        % check convergence
    if ( error <= tol ), break, end;
end

if ( error > tol ), flag = 1; end;                 % converged

t2 = toc;

%% solve with preconditioner
function r = Psolve(M1,M2,M3,M4,K1,K2,rres,method,kapp,l,multdata,pow,supg,Amat1,Amat2,smoother)

switch method
        case 'backslash'

            r1 = (M1\rres(1:l));
            r2 = (M2\rres(l+1:l+l));
            switch kapp
                case 1
                    r3 = K1\(M2*(K2\rres(l+l+1:end)));
                case 0
                    %             r3 = M3\(M1*(M4\rres(l+l+1:end)));
                    r3 = M4\rres(l+l+1:end);
            end
        case 'multigrid'
            %         zed1 = pcg4(M2,v2(1:m))./(2*beta);
            %         zed2 = pcg4(M,v2(m+1:m+l));
            r1 = chebit2(M1,rres(1:l),20);%./(2*beta);
%             r2 = chebit2(M2,rres(l+1:l+l),20);
%                     zed1 = M1\v2(1:l);%./(2*beta);
            r2 = M2\rres(l+1:2*l);

            r3a = sparse(l,1);
            r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
            r3a = multigrid_sd(K2,rres(l+l+1:end),multdata,r3a,pow,supg,Amat2,smoother);
            r3b = M2*r3a;
            r3 = sparse(l,1);
            r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
            r3 = multigrid_sd(K1,r3b,multdata,r3,pow,supg,Amat1,smoother);
end
r = [r1;r2;r3];
%%
% END of gmres.m
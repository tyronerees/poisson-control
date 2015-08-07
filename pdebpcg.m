function  [x,iter,t2]=pdebpcg(A,b,prob_setup,def_soln,def_setup)

% CG algorithm with inner product defined by H and preconditioners for the
% PDE constrained optimisation problem: A1 for the 1,1 block
%                                       A2 for the 2,2 block
%                                       S0 Schur complement preconditioner
%
% A1,A2,S0 are assumed to be function handles so that multgrid or Chebychev
% semi-iteration can be used as a preconditioner.
%
% M.Stoll 2008
% Adapted by T. Rees 2010

nu = prob_setup.nu;
ny = prob_setup.ny;
lb = 2*ny+nu;

tic

initksolve

x0 = zeros(lb,1); % intial guess
r0 = (b-A*x0);  % initial residual
n2b = norm(r0); % norm of inital residual


% pj2 = M1\r0; % direction
%steps = 10; % number of csi steps
%scale = .9; % scaling for inner product

pj = presolve(A,r0,prob_setup,def_soln,multdata);

% pj(1:nu,1) = chebsemiit(A(1:nu,1:nu),v(1:nu),def_soln.ucheb,prob_setup.dim,prob_setup.yelt);
% % pj(1:n1,1) = csi(scale*A(1:n1,1:n1),r0(1:n1,1),steps); % direction
% pj(nu+1:nu+ny,1) = 
% pj(n1+1:n1+n2,1) = csi(scale*A(n1+1:n1+n2,n1+1:n1+n2),r0(n1+1:n1+n2,1),steps); % direction
% 
% pj(n1+n2+1:n1+n2+m,1) = S0(A(n1+n2+1:n1+n2+m,1:n1+n2)*pj(1:n1+n2,1)-r0(n1+n2+1:n1+n2+m,1),A,n1,n2,m); % direction

bnd =find(b(1:nu)==0);

n2bb = norm(pj);
urj = r0; % unpreconditioned residual
rj = pj; % residual
xj  = x0; % solution

% tol = 1e-6; % tolerance
j   = 1;

resj = r0;
resvec(j) = norm(resj);

relres  = resvec(j)/n2b; % relative residual
itmax   = min(10000,lb); % maximal iteration number
tol = 1e-6; % desired tolerance for rel res
it = 1; % startcounter

tmp1 = A*rj;

rHrj = rj(1:nu+ny)'*tmp1(1:nu+ny)-rj(1:nu+ny)'*urj(1:nu+ny)-rj(nu+ny+1:lb)'*urj(nu+ny+1:lb);


% rHrj = [A(1:n1,1:n1)*rj(1:n1,1)-urj(1:n1,1);... % first block
%     A(n1+1:n1+n2,n1+1:n1+n2)*rj(n1+1:n1+n2,1)-urj(n1+1:n1+n2,1);... % second block
%     A(n1+n2+1:n1+n2+m,1:n1+n2)*rj(1:n1+n2,1)-urj(n1+n2+1:n1+n2+m,1)]'*rj % last block and multiplication


apj1     = tmp1;
while (relres > tol) & (it < itmax)
%     apj1     = A*pj; % matrix vector multiplication
    % P^{-1}Ap_j
    apj2 = presolve(A,apj1,prob_setup,def_soln,multdata);
%     apj2(1:n1,1) = csi(scale*A(1:n1,1:n1),apj1(1:n1,1),steps); % first block
%     apj2(n1+1:n1+n2,1) = csi(scale*A(n1+1:n1+n2,n1+1:n1+n2),apj1(n1+1:n1+n2,1),steps); % second block
%     apj2(n1+n2+1:n1+n2+m,1) = S0(A(n1+n2+1:n1+n2+m,1:n1+n2)*apj2(1:n1+n2,1)-apj1(n1+n2+1:n1+n2+m,1),A,n1,n2,m); % third block
    
    alphaj  = rHrj/(apj1(1:nu+ny)'*apj2(1:nu+ny)-pj(1:nu+ny)'*apj1(1:nu+ny)-pj(nu+ny+1:lb)'*apj1(nu+ny+1:lb));
    
    xjplus1 = xj+alphaj*pj; % update current solution
    urjplus1 = urj-alphaj*apj1; % update unpreconditioned residual
    rjplus1 = rj-alphaj*apj2; % update preconditioned residual

    tmp1 = A*rjplus1;
%     tmp1(1:n1+n2) = A(1:n1+n2,1:n1+n2)*rjplus1(1:n1+n2);
    rHrjplus1 = rjplus1(1:nu+ny)'*tmp1(1:nu+ny)-rjplus1(1:nu+ny)'*urjplus1(1:nu+ny)-rjplus1(nu+ny+1:lb)'*urjplus1(nu+ny+1:lb);
      
%     tmp1(1:n1+n2) = tmp1(1:n1+n2)+A(1:n1+n2,n1+n2+1:n1+n2+m)*rjplus1(n1+n2+1:n1+n2+m);
%     tmp1(n1+n2+1:n1+n2+m) = A(n1+n2+1:n1+n2+m,1:n1+n2)*rjplus1(1:n1+n2);
    
    betaj   = rHrjplus1/rHrj; % Compute beta

    pjplus1 = rjplus1+betaj*pj; % update direction

    apj1 = tmp1+betaj*apj1;
    % Update for next iteration
    pj = pjplus1;
    xj = xjplus1;
    rj = rjplus1;
    urj = urjplus1;
    rHrj = rHrjplus1;
    j  = j+1;
    resvec(j) = norm(urj); % Residual norm

    relres  = resvec(j)/n2bb; % Relative residual
    %     relres  = resvec(j)/n2b;
    it = it +1;

end
x = xj;
flag = [];
iter = it;
t2 = toc;
% fprintf('bpcg converged in %d iterations and took %4.3f seconds\n',iter,t2);
switch lower(def_soln.kmethod)
    case 'amg'
        inform = mi20_finalize;
    case 'amgsch'
        inform = mi20_finalize;
end 
if strcmp(prob_setup.eqn,'stokes')==1
    switch def_soln.kmethod_stokes
        case 'amg'
            inform = mi20_finalize;
    end
end

function z = presolve(A,v,prob_setup,def_soln,multdata)
nu = prob_setup.nu;
ny = prob_setup.ny;
lb = 2*ny+nu;

if strcmp(prob_setup.eqn,'stokes')==1
    nvel = prob_setup.nvel;
    nvel2 = nvel/2;
    np = prob_setup.np;
    zed1a = chebsemiit2(def_soln.scale*A(nu+1:nu+nvel2,nu+1:nu+nvel2),v(nu+1:nu+nvel2),def_soln.velcheb,prob_setup.dim,prob_setup.velelt);
    zed1b = chebsemiit2(def_soln.scale*A(nu+nvel2+1:nu+nvel,nu+nvel2+1:nu+nvel),v(nu+nvel2+1:nu+nvel),def_soln.velcheb,prob_setup.dim,prob_setup.velelt);
    zed1c = chebsemiit2(def_soln.scale*A(nu+nvel+1:nu+nvel+np,nu+nvel+1:nu+nvel+np),v(nu+nvel+1:nu+nvel+np),def_soln.pcheb,prob_setup.dim,prob_setup.pelt);
    zed2 = [zed1a;zed1b;zed1c];
    zed1 = chebsemiit2(def_soln.scale*A(1:nu,1:nu),v(1:nu),def_soln.ycheb,prob_setup.dim,prob_setup.uelt);
else
    zed1 = chebsemiit(def_soln.scale*A(1:nu,1:nu),v(1:nu),def_soln.ucheb,prob_setup.dim,prob_setup.yelt);
    zed2 = chebsemiit(def_soln.scale*A(nu+1:nu+ny,nu+1:nu+ny),v(nu+1:nu+ny),def_soln.ycheb,prob_setup.dim,prob_setup.uelt);
end

rhs = A(nu+ny+1:end,1:nu+ny)*[zed1;zed2]-v(nu+ny+1:end);
switch lower(def_soln.kmethod)
    case 'gmg'
        zed3a = sparse(ny,1);
        for i = 1:def_soln.vcyc
            zed3a = mgvcyc(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),rhs,multdata,zed3a,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
        end
        zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
        zed3 = sparse(ny,1);
        for i = 1:def_soln.vcyc
            zed3 = mgvcyc(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),zed3b,multdata,zed3,multdata(1).pow,multdata(1).spost,multdata(1).spre,multdata(1).dim);
        end
    case 'amg'
        zed3a = mi20_precondition(rhs);
        zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
        zed3 = mi20_precondition(zed3b);
    case 'bslash'
        zed3 = A(nu+ny+1:2*ny+nu,nu+1:nu+ny)\(A(nu+1:nu+ny,nu+1:nu+ny)*(A(nu+1:nu+ny,nu+ny+1:2*ny+nu)\rhs));
    case 'stokes'
        nvel = prob_setup.nvel;
        np = prob_setup.np;
        zed3a = stokesapprox(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),rhs,A(nu+nvel+1:nu+nvel+np,nu+nvel+1:nu+nvel+np)*def_soln.alp/prob_setup.delta,def_soln,prob_setup);
        zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
        zed3 = stokesapproxt(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),zed3b,A(nu+nvel+1:nu+nvel+np,nu+nvel+1:nu+nvel+np)*def_soln.alp/prob_setup.delta,def_soln,prob_setup);       
end

z = [zed1; zed2; zed3]; 

% function [y2] = csi(A,b,i)
% 
% % Generates the ith chebyshev semi iteration vector using Jacobi iteration.
% %
% % Input : A - Matrix to be solved
% %       : b - Right hand side
% %       : i - number of iterations of the Cheyshev semi iteration
% %
% % Output : y2 - approximation to y satisfying A*y = b;
% 
% dim = 2;                     % Switch for the dimesion of the problem 
% 
% % tic,
% l = length(A);
% 
% I = speye(l);                     % define identity
% % invM = diag(1./diag(A));                % Compute inverse once
% invM = spdiags(1./diag(A), 0, l, l);
% 
% % keyboard
% switch dim
%     case 2
%         w = 4/5;                                % define the relaxation parameter
%     case 3
%         w = 4/7;
% end
%         
% g = w*invM*b;                           % define relaxed Jacobi g
% S = (I-w*invM*A);                       % define the iteration Matrix
% 
% y1 = zeros(length(b),1);
% y2 = S*y1 + g;                          % x_(m+1) = Sx_m + g
% 
% switch dim
%     case 2
%         mu = 5/4;                               % largest eig for this relaxation parameter 
%     case 3                                   % using the mass matrix as A (see
%         mu = 14/13;                            % notes 30-11-07)
% end
% 
% ck0 = 2;
% ck1 = mu;
% 
% for j = 2:i
%    ck2 = mu*ck1-0.25*ck0;
% 
%      w = mu*ck1/ck2;
%     y3 = w*(S*y2+g-y1)+y1;  
%     
%     y1 = y2;    y2 = y3;
%    ck0 = ck1;
%    ck1 = ck2;
% end
% % t = toc;
% function [y] = S0(x,A,n1,n2,m)
% % 
% % beta = A(n1+1,n1+1)/A(1,1);
% 
% % y=x;
% % 
% % y=A(n1+n2+1:n1+n2+m,1:n1)\x;
% % y=A(1:n1,1:n1)*y;
% % y=A(n1+n2+1:n1+n2+m,1:n1)\y;
% x =  mi20_precondition(x);
% x =  A(1:n1,1:n1)*x;
% y =  mi20_precondition(x);
% 
% % % % matrix = [ A(1:n1,1:n1) -A(n1+n2+1:n1+n2+m,1:n1)';-A(n1+n2+1:n1+n2+m,1:n1) -beta*A(1:n1,1:n1)];
% % % % rhs = [zeros(n1,1);x];
% % % % [xx,flag,relres,iter,resvec]=bpcgcsiamg2(matrix,rhs,n1,m,1e-5,.9);
% % % % y=-xx(n1+1:n1+m,1);


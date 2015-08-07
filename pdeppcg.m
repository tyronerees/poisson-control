function [xx,iter,t2] = pdeppcg(A,b,prob_setup,def_soln,def_setup)
%PPCG Projected Preconditioned Conjugate Gradients
%   My implementation of the PPCG algorithm as presented by Gould, Hribar
%   and Nocedal (2000 - Algorithm II)

nu = prob_setup.nu;
ny = prob_setup.ny;
lb = 2*ny+nu;


tic,

initksolve

H = A(1:nu+ny,1:nu+ny);
%   Choose a point x satisfying Ax=d, i.e. -Mx_1 = d with x_2 = 0
x = sparse(nu+ny,1);
% x(1:n) = -M\b(end-n+1:end);
% x(1:nu) = -pcg4(-A(nu+ny+1:2*ny+nu,1:nu),b(nu+ny+1:end));

switch def_soln.mmethod
    case 'bslash'
        x(1:nu) = A(nu+ny+1:2*ny+nu,1:nu)\b(nu+ny+1:end);
    case 'chebit'
        x(1:nu) = chebsemiit(A(nu+ny+1:2*ny+nu,1:nu),b(nu+ny+1:end),def_soln.chebit_ppcg,prob_setup.dim,prob_setup.yelt);
end


%   Compute r = Hx-b
r = -b(1:nu+ny);
r(1:nu) = r(1:nu)+A(1:nu,1:nu)*x(1:nu);
iter = 0;
% keyboard
%   Schilders factorization
rt = r;


% P = [sparse(n,n) sparse(n,n) -speye(n,n);...

% % % % yy = P\[r;sparse(n,1)];            % change to a Schilders factorization
% yy = schfact([r;sparse(n,1)],M,K,beta);

switch lower(def_soln.conmethod)
    case 'diag'
        fprintf('Solving with diagonal (1,1) block  chosen \n')
        yy = prediag(A,r,prob_setup);
    case 'pblktri'
        yy = preblktri(A,r,prob_setup,def_soln,def_setup.beta,multdata);
end


g = yy(1:nu+ny); v = yy(nu+ny+1:end);
p=-g; y=v; 

r = r-A(1:nu+ny,nu+ny+1:end)*y;

sig = r'*g;
nrmtol = def_soln.tol*sig;
sig1 = sig*10;
sig;
% keyboard
if def_soln.twonorm == 1 % test on the residual of 2 norm
    test = 2;
else
    test = sig; % test on the the ppg norm
end
if def_soln.resvec == 1
    resvec = [];
end
% while (res>=def_soln.tol)
while (test>=nrmtol)

    Hp = H*p;
    alp = r'*g/(p'*Hp);
    x = x+alp*p;
    r1 = r+alp*Hp;
    switch lower(def_soln.conmethod)
        case 'diag'
            yy = prediag(A,r1,prob_setup);
        case 'pblktri'
            yy = preblktri(A,r1,prob_setup,def_soln,def_setup.beta,multdata);
    end

    g1 = yy(1:nu+ny);
    v1 = yy(nu+ny+1:end);
    sig1 = r1'*g1;
    bet = sig1/sig;
    p = -g1+bet*p;
    g = g1;
    r = r1-A(1:nu+ny,nu+ny+1:end)*v1;
    iter = iter+1;
    sig = sig1;
%     pause
% res = 
    if def_soln.twonorm == 1
       test = norm(b-A*[x;v1]);  % test on 2 norm
    else
       test = sig;       % test on ppcg norm
    end
    if def_soln.resvec == 1
        resvec = [resvec; iter sig norm(b-A*[x;v1])];
    end
end
xx = [x;v1];

t2 = toc;

if def_soln.resvec == 1
    figure(17)
    semilogy(resvec(:,1),resvec(:,2),'-','markersize',5,'linewidth',2);
    hold on
    semilogy(resvec(:,1),resvec(:,3),'--','markersize',5,'linewidth',2);
    hold off
    legend('PPCG norm','2 norm','location','bestOutside')
    cd ..
    cd Poisson
    print('-f17','-depsc','ppcg_resvec');
    cd ..
    cd 'Thesis code'
end

switch lower(def_soln.kmethod)
    case 'amg'
        inform = mi20_finalize;
end 

%% solution functions
function yy = preblktri(A,r,prob_setup,def_soln,beta,multdata)
nu = prob_setup.nu;
ny = prob_setup.ny;

r1 = r(1:nu);
r2 = r(nu+1:nu+ny);
r3 = sparse(ny,1);

switch def_soln.mmethod
    case 'bslash'
        y3 = A(nu+ny+1:2*ny+nu,1:nu)\r1;
    case 'chebit'
        y3 = chebsemiit(A(nu+ny+1:2*ny+nu,1:nu),r1,def_soln.chebit_ppcg,prob_setup.dim,prob_setup.yelt);
end

switch def_soln.kmethod
    case 'bslash'
        y2 = 1/beta*(A(nu+1:nu+ny,nu+ny+1:2*ny+nu)\(A(nu+1:nu+ny,nu+1:nu+ny)*(A(nu+ny+1:2*ny+nu,nu+1:nu+ny)\(r2-A(nu+1:nu+ny,nu+ny+1:2*ny+nu)*y3))));
    case 'amg'
        r2a = r2-A(nu+1:nu+ny,nu+ny+1:2*ny+nu)*y3;
        y2a = mi20_precondition(r2a);
        y2b = beta*A(nu+1:nu+ny,nu+1:nu+ny)*y2a;
        y2 = mi20_precondition(y2b);
    case 'gmg'
        r2a = r2-A(nu+1:nu+ny,nu+ny+1:2*ny+nu)*y3;
        y2a = sparse(ny,1);
        for i = 1:def_soln.vcyc
            y2a = mgvcyc(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),r2a,multdata,y2a,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
        end
        y2b = beta*A(nu+1:nu+ny,nu+1:nu+ny)*y2a;
        y2 = sparse(ny,1);
        for i = 1:def_soln.vcyc
            y2 = mgvcyc(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),y2b,multdata,y2,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
        end 
end

switch def_soln.mmethod
    case 'bslash'
        y1 = A(nu+ny+1:2*ny+nu,1:nu)\(r3-A(nu+ny+1:2*ny+nu,nu+1:nu+ny)*y2);
    case 'chebit'
        r1a = (r3-A(nu+ny+1:2*ny+nu,nu+1:nu+ny)*y2);
        y1 = chebsemiit(A(nu+ny+1:2*ny+nu,1:nu),r1a,def_soln.chebit_ppcg,prob_setup.dim,prob_setup.yelt);
end
yy = [y1;y2;y3];
%%
function yy = prediag(A,r,prob_setup)

nu = prob_setup.nu;
ny = prob_setup.ny;

P = A;
P(1:nu+ny,1:nu+ny) = diag(diag(A(1:nu+ny,1:nu+ny)));
 
% vv = [r;sparse(n,1)];

yy = P\[r;sparse(ny,1)];


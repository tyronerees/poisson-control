function [u,iter,t2] = pdeminres(A,b,prob_setup,def_soln,def_setup,bdy_set)

nu = prob_setup.nu;
ny = prob_setup.ny;
lb = 2*ny+nu;

% Mu = A(1:nu,1:nu);
% My = A(nu+1:nu+ny,nu+1:nu+ny);
% Muy = A(nu+ny+1:2*ny+nu,1:nu);
% K = A(nu+ny+1:2*ny+nu,nu+1:nu+ny);

tic,

initksolve;

gam0 = 1;
u = sparse(lb,1);

% eK = sort(real(eig(full(A(nu+ny+1:2*ny+nu,nu+1:nu+ny)))));
% eM = sort(real(eig(full(A(nu+1:nu+ny,nu+1:nu+ny)))));
% Kl = eK(1);
% KL = eK(end);
% Ml = eM(1);
% ML = eM(end);
% def_soln.DEL = (ML/Kl)^2;
% def_soln.del = (Ml/KL)^2;

v0 = sparse(lb,1);
w0 = sparse(lb,1); w1 = sparse(lb,1);
v1 = b-A*u;

z1 = presolve(A,v1,prob_setup,def_soln,multdata,def_setup);

gam1 = sqrt(z1'*v1);
eta = gam1; s0 = 0; s1 = 0; c0 = 1; c1 = 1;
eta0 = eta;
iter = 0;
if def_soln.twonorm == 1 % test on the residual of 2 norm
    test = 2;
else
    test = eta; % test on the the ppg norm
end
while abs(test)>def_soln.tol*eta0      
    z1 = z1/gam1;
    Az1 = A*z1;
    del = Az1'*z1;  
    v2 = Az1 - (del/gam1)*v1-(gam1/gam0)*v0;

    z2 = presolve(A,v2,prob_setup,def_soln,multdata,def_setup);

    gam2 = sqrt(z2'*v2);
    alp0 = c1*del-c0*s1*gam1;
    alp1 = sqrt(alp0^2+gam2^2);
    invse = 1/alp1;
    alp2 = s1*del + c0*c1*gam1;
    alp3 = s0*gam1;
    c0 = c1; s0 = s1; c1 = alp0*invse; s1 = gam2*invse;
    w2 = (z1 - alp3*w0-alp2*w1)*invse;
    u = u + c1*eta.*w2;
    eta = -s1*eta;
    v0 = v1; v1=v2; z0=z1; z1=z2; w0=w1; w1=w2;  gam0 = gam1; gam1 = gam2;
    iter = iter +1;
    if def_soln.twonorm == 1
       test = norm(b-A*u);  % test on 2 norm
    else
       test = eta;       % test on ppcg norm
    end
end
t2 = toc;
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
%%%%%%%%%%%%%%%%%%%%%%%%
function z = presolve(A,v,prob_setup,def_soln,multdata,def_setup)
nu = prob_setup.nu;
ny = prob_setup.ny;
lb = 2*ny+nu;

switch lower(def_soln.mmethod)
    case 'chebit'
        if strcmp(prob_setup.eqn,'stokes')==1
            nvel = prob_setup.nvel;
            nvel2 = nvel/2;
            np = prob_setup.np;
            zed2a = chebsemiit2(A(nu+1:nu+nvel2,nu+1:nu+nvel2),v(nu+1:nu+nvel2),def_soln.velcheb,prob_setup.dim,prob_setup.velelt);
            zed2b = chebsemiit2(A(nu+nvel2+1:nu+nvel,nu+nvel2+1:nu+nvel),v(nu+nvel2+1:nu+nvel),def_soln.velcheb,prob_setup.dim,prob_setup.velelt);
            zed2c = chebsemiit2(A(nu+nvel+1:nu+nvel+np,nu+nvel+1:nu+nvel+np),v(nu+nvel+1:nu+nvel+np),def_soln.pcheb,prob_setup.dim,prob_setup.pelt);
            zed2 = [zed2a;zed2b;zed2c];
            zed1a = chebsemiit2(A(1:nu/2,1:nu/2),v(1:nu/2),def_soln.ycheb,prob_setup.dim,prob_setup.uelt);
            zed1b = chebsemiit2(A(nu/2+1:nu,nu/2+1:nu),v(nu/2+1:nu),def_soln.ycheb,prob_setup.dim,prob_setup.uelt);
            zed1 =[zed1a;zed1b];
        else
            zed1 = chebsemiit2(A(1:nu,1:nu),v(1:nu),def_soln.ucheb,prob_setup.dim,prob_setup.uelt);
            zed2 = chebsemiit2(A(nu+1:nu+ny,nu+1:nu+ny),v(nu+1:nu+ny),def_soln.ycheb,prob_setup.dim,prob_setup.yelt);
        end
    case 'bslash'
        zed1 = A(1:nu,1:nu)\v(1:nu);
        zed2 =A(nu+1:nu+ny,nu+1:nu+ny)\v(1+nu:nu+ny);
end

switch lower(def_soln.kmethod)
    case 'gmg'
        if strcmp(prob_setup.eqn,'advdiff')==1
            zed3a = sparse(ny,1);
            for i = 1:def_soln.vcyc
                zed3a = mgvcyc_sd1(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),v(nu+ny+1:end),multdata,zed3a,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
            end
            zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
            zed3 = sparse(ny,1);
            for i = 1:def_soln.vcyc
                zed3 = mgvcyc_sd2(A(nu+1:nu+ny,nu+ny+1:2*ny+nu),zed3b,multdata,zed3,multdata(1).pow,multdata(1).spost,multdata(1).spre,multdata(1).dim);
            end
        else
            
            zed3a = sparse(ny,1);
            for i = 1:def_soln.vcyc
                zed3a = mgvcyc(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),v(nu+ny+1:end),multdata,zed3a,multdata(1).pow,multdata(1).spre,multdata(1).spost,multdata(1).dim);
            end
            zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
            zed3 = sparse(ny,1);
            for i = 1:def_soln.vcyc
                zed3 = mgvcyc(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),zed3b,multdata,zed3,multdata(1).pow,multdata(1).spost,multdata(1).spre,multdata(1).dim);
            end
        end
    case 'amg'
        zed3a = mi20_precondition(v(nu+ny+1:end));
        zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
        zed3 = mi20_precondition(zed3b);
    case 'bslash'
        zed3 = A(nu+1:nu+ny,nu+ny+1:2*ny+nu)\(A(nu+1:nu+ny,nu+1:nu+ny)*(A(nu+ny+1:2*ny+nu,nu+1:nu+ny)\v(nu+ny+1:end)));
    case 'stokes'
        nvel = prob_setup.nvel;
        np = prob_setup.np;
        if def_soln.exactschur == 1;
            S0 = A(nu+ny+nvel+1:ny+nu+nvel+np,nu+1:nu+nvel)*(A(nu+ny+1:nu+ny+nvel, nu+1:nu+nvel)\A(nu+1:nu+nvel,nu+ny+nvel+1:ny+nu+nvel+np));
            %keyboard
        else
            S0 = A(nu+nvel+1:nu+nvel+np,nu+nvel+1:nu+nvel+np)*def_soln.alp/prob_setup.delta;
        end
        zed3a = stokesapprox(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),v(nu+ny+1:end),S0,def_soln,prob_setup);
        zed3b = A(nu+1:nu+ny,nu+1:nu+ny)*zed3a;
        zed3 = stokesapproxt(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),zed3b,S0,def_soln,prob_setup);       
end


oo = 1;
switch oo
    case 1
        alp = 1;
    case 2
        alp = 1+(1/(2*def_setup.beta))*(def_soln.DEL+def_soln.del);
    case 3
        alp = 1/def_setup.beta;
    case 4
        alp = (1+2/(2*def_setup.beta));  
    case 5
        ll = 0.01;
        kk=(def_setup.beta)^0.5;
        alp = (1/((kk+ll)*def_setup.beta))*(ll*def_soln.DEL+kk*def_soln.del);
end
zed3 = zed3/alp;

z = [zed1; zed2; zed3]; 


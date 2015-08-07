function [u,iter,t2] = stokesminres(A,b,M,st_setup)

nvel = st_setup.vel;
np = st_setup.p;
lb = 2*nvel+np;

% K1 = A(1:nvel,1:nvel);
% K = A(nu+1:nu+ny,nu+1:nu+ny);
% Muy = A(nu+ny+1:2*ny+nu,1:nu);
% K = A(nu+ny+1:2*ny+nu,nu+1:nu+ny);

tic,

switch lower(st_setup.ksolve)
    case 'gmg'
        if st_setup.loadprol == 0
                    multdata = genmultdata(st_setup.pow+1,'dirichlet');
        else
             load multdataQ1_dir
        end
        multdata(1).pow = st_setup.pow;
    case 'amg'
        control = mi20_control;
        control.v_iterations = st_setup.vcyc;
        inform = mi20_setup(A(nvel+1:2*nvel,nvel+1:2*nvel),control);
end
if strcmp(st_setup.ksolve,'gmg') == 0
    multdata = [];
end



gam0 = 1;
u = sparse(lb,1);


v0 = sparse(lb,1);
w0 = sparse(lb,1); w1 = sparse(lb,1);
v1 = b-A*u;

z1 = presolve(A,v1,M,st_setup,multdata);

gam1 = sqrt(z1'*v1);
eta = gam1; s0 = 0; s1 = 0; c0 = 1; c1 = 1;
eta0 = eta;
iter = 0;
test = eta; % test on the the ppg norm

while abs(test)>st_setup.tol*eta0      
    z1 = z1/gam1;
    Az1 = A*z1;
    del = Az1'*z1;  
    v2 = Az1 - (del/gam1)*v1-(gam1/gam0)*v0;

    z2 = presolve(A,v2,M,st_setup,multdata);

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
    test = eta;       % test on ppcg norm
end
t2 = toc;
switch lower(st_setup.ksolve)
    case 'amg'
        inform = mi20_finalize;
    case 'amgsch'
        inform = mi20_finalize;
end 

%%%%%%%%%%%%%%%%%%%%%%%%
function z = presolve(A,v,M,prob_setup,multdata)
nvel = prob_setup.vel;
np = prob_setup.p;
lb = 2*nvel+np;

switch lower(prob_setup.ksolve)
    case 'bslash'
        zed1 = A(1:nvel,1:nvel)\v(1:nvel);
        zed2 =A(nvel+1:2*nvel,nvel+1:2*nvel)\v(1+nvel:2*nvel); 
    case 'gmg'
        zed1 = sparse(nvel,1);
        for i = 1:prob_setup.vcyc
            zed1 = mgvcyc(A(1:nvel,1:nvel),v(1:nvel),multdata,zed1,multdata(1).pow+1,prob_setup.pre,prob_setup.pre,2);
        end
        zed2 = sparse(nvel,1);
        for i = 1:prob_setup.vcyc
            zed2 = mgvcyc(A(nvel+1:2*nvel,nvel+1:2*nvel),v(nvel+1:2*nvel),multdata,zed2,multdata(1).pow+1,prob_setup.pre,prob_setup.pre,2);
        end
    case 'amg'
        zed1 = mi20_precondition(v(1:nvel));
        zed2 = mi20_precondition(v(nvel+1:2*nvel));
    
end

switch lower(prob_setup.msolve)
    case 'chebit'
        zed3 = chebsemiit2(M,v(2*nvel+1:end),prob_setup.cheb,2,prob_setup.pelt);
    case 'bslash'
        zed3 = M\v(2*nvel+1:end);
end


z = [zed1; zed2; zed3]; 


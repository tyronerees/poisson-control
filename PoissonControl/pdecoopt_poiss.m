clear all
hold off

def_setup = set_def_setup;

def_setup.pow  =5;
def_setup.beta = 1e-4;
def_setup.bc = 'dirichlet';
def_setup.ob = 1;
def_setup.plots = 1;
% def_setup.type = 'bound2d';
 def_setup.type = 'dist2d';
%def_setup.type = 'dist2d';

if strcmp(def_setup.type,'dist3d') == 1
    def_setup.plots = 0;
end
 
fprintf('Setting up the system...\n')
switch def_setup.type
    case 'dist2d'
        fprintf('Setting up the 2D distributed control problem...\n')
        [A,b,bdy_set,ubdy,uhat,def_setup,prob_setup] = pdecoDriver(def_setup);
    case 'dist3d'
        fprintf('Setting up the 3D distributed control problem...\n')
        [A,b,bdy_set,ubdy,uhat,def_setup,prob_setup] = pdecoDriver_3d(def_setup);
    case 'bound2d'
        fprintf('Setting up the 2D boundary control problem...\n')
        [A,b,bdy_set,ubdy,uhat,def_setup,prob_setup] = pdecoDriver_bc(def_setup);
        
end
        

pow = def_setup.pow;
N = 2^pow; h = 1/N;
beta = def_setup.beta;

%% solve
[def_soln] = set_def_soln(prob_setup);
def_soln.method = 'minres';
def_soln.kmethod = 'gmg';      % method for approx stiff matrix
def_soln.mmethod = 'chebit';      % method for approx mass matrix
def_soln.ucheb = 5;
def_soln.ycheb = 5;
def_soln.vcyc = 1;
def_soln.tol = 1e-6;
def_soln.loadprol = 1;
def_soln.conmethod = 'pblktri';  % diag or pblktri
cd ..
W = consolve(A,b,prob_setup,def_soln,def_setup,bdy_set);
cd PoissonControl

%% get optimal control, state and adjoint
nu = prob_setup.nu;
ny = prob_setup.ny;

F = W(1:nu);
U = W(nu+1:nu+ny);
lam = W(nu+ny+1:end);

if prob_setup.dim == 2
    v = 1:(N+1)^2;
    Un = sparse(1,(N+1)^2);
    Fn = sparse(1,(N+1)^2);
    vb = 1:4*N;
elseif prob_setup.dim == 3
    v = 1:(N+1)^3;
    Un = sparse(1,(N+1)^3);
    Fn = sparse(1,(N+1)^3);
    vb = 1:6*N^2+2;
end

v(bdy_set.dirn) = [];
vb(bdy_set.dbynodes) = [];
Un(bdy_set.dirn) = ubdy;
Un(v) =U;
if strcmp(def_setup.type,'dist2d') == 1
 Fn = sparse(1,(N+1)^2);
 Fn(bdy_set.dirn) = 0;
 Fn(v) =F;
 f = full(reshape(Fn,N+1,N+1));
end

 %make full so as not to confuse the plotting tool




%% plot the state
if def_setup.plots == 1;
    u = full(reshape(Un,N+1,N+1));
    x=0:h:1;
    y=x;
    subplot(2,2,2)
    surf(x,y,u);
    title({['Plot of state u, N = 2^' int2str(pow) ', \beta = ' num2str(beta) ];['Distributed control, ' def_setup.bc ' BCs']});
end

%% calculate the difference in u to \hat(u) and the value of J(u,f)
uhat(bdy_set.dirn) = [];

Mu = A(1:nu,1:nu);
My = A(nu+1:nu+ny,nu+1:nu+ny);

if def_setup.ob == 3;
    if strcmp(def_setup.type,'dist3d')==1
        excon = h^3/8;
    elseif strcmp(def_setup.type,'dist2d')==1
        excon = h^2/4;
    end
elseif def_setup.ob == 1;
    if strcmp(def_setup.type,'dist2d')==1
        excon = (h/6*(h^2-3*h+3))^2;
    elseif strcmp(def_setup.type,'dist3d')==1
        excon = (h/6*(h^2-3*h+3))^3;
    end
else
    excon = uhat*My*uhat';
end
        

erroru = 0.5*(U'*My*U-2*U'*b(1:nu)+excon);
cfunc = erroru+0.5*beta*F'*Mu*F;
fprintf('||u-uhat||^2 is %d \n',erroru);
fprintf('Cost function is %d \n',cfunc);

%% Plot the control
if (def_setup.plots == 1)&&(strcmp(def_setup.type,'dist2d') == 1);
    subplot(2,2,3)
    surf(x,y,f)
    title({['Plot of control, f, N = 2^' int2str(pow) ', \beta = ' num2str(beta)];['Distributed control, ' def_setup.bc ' BCs']});
end


%% Test right hand side gives the same solution
% switch method 
%     case 1
%         ufor = FPoissSolve(F,h,(bdy_set.dirn),ubdy);
%         norm(ufor-U)
% end
if strcmp(def_setup.type,'dist2d') == 1
ufor = FPoissSolve(Fn,h,(bdy_set.dirn),ubdy);
norm(ufor-Un')
end


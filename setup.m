% Example 2 - Distributed control, NBC
clear all
close all

cd PoissonControl

def_setup = set_def_setup;
def_setup.pow  =5;

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
[def_soln] = set_def_soln(prob_setup);

def_soln.method = 'bpcg';           % bpcg, minres, backslash, gmres, ppcg  
def_soln.kmethod = 'gmg';           % stiffness matrix approximation:
                                    % gmg, bslash (amg -- if hsl mi20 installed)
def_soln.mmethod = 'chebit';        % mass matrix approximation:
                                    % chebit, bslash
def_soln.gmgpre = 3;                % no of pre smoothing steps for mgrid (int >= 0)
def_soln.gmgpost = 0;               % no of post smoothing steps for mgrid (int >= 0)
def_soln.ucheb = 5;                 % no of cheb semi it steps for control mass matrix (int >= 0)
def_soln.ycheb = 5;                 % no of cheb semi it steps for control mass matrix (int >= 0)
def_soln.scale = 0.9;               % scaling parameter for bpcg

cd ..

fprintf('GMG...\n')
[Wg,iter,t2] = consolve(A,b,prob_setup,def_soln,def_setup);
F = Wg(1:prob_setup.nu);
U = Wg(prob_setup.nu+1:prob_setup.nu+prob_setup.ny);

nu = prob_setup.nu;
ny = prob_setup.ny;
N = 2^def_setup.pow; h = 1/N;

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

if def_setup.plots == 1;
    u = full(reshape(Un,N+1,N+1));
    x=0:h:1;
    y=x;
    subplot(2,2,2)
    surf(x,y,u);
    title({['Plot of state u, N = 2^' int2str(def_setup.pow) ', \beta = ' num2str(def_setup.beta) ];['Distributed control, ' def_setup.bc ' BCs']});
end
%% Plot the control
if (def_setup.plots == 1)&&(strcmp(def_setup.type,'dist2d') == 1);
    subplot(2,2,3)
    surf(x,y,f)
    title({['Plot of control, f, N = 2^' int2str(def_setup.pow) ', \beta = ' num2str(def_setup.beta)];['Distributed control, ' def_setup.bc ' BCs']});
end



cd ..
cd ..
cd 'Thesis code'
cd PoissonControl

def_setup = set_def_setup;
def_setup.plots = 0;

maxn = 6;
rslts = zeros(maxn-1,9);
for i = 2:maxn
    i
    rslts(i-1,1) = i;
    def_setup.pow  =i;
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
    def_soln.method = 'bpcg';
    def_soln.mmethod = 'bslash';      % method for approx mass matrix
    
    cd ..
    def_soln.kmethod = 'gmg';      % method for approx stiff matrix
    def_soln.vcyc = 1;
    fprintf('1 V cycle, 1 pre 1 post...\n')
    fprintf('...\n')
    def_soln.gmgpre = 1;
    def_soln.gmgpost = 1;
    [W,iter,t2] = consolve(A,b,prob_setup,def_soln,def_setup);
    
    rslts(i-1,2) = t2;
    rslts(i-1,3) = iter;
    
    fprintf('1 V cycle, 2 pre 2 post...\n')
    def_soln.gmgpre = 2;
    def_soln.gmgpost = 2;
    [W,iter,t2] = consolve(A,b,prob_setup,def_soln,def_setup);
    
    rslts(i-1,4) = t2;
    rslts(i-1,5) = iter;
    
    fprintf('1 V cycle, 2 pre 0 post...\n')
    fprintf('...\n')
    def_soln.gmgpre = 2;
    def_soln.gmgpost = 0;
    [W,iter,t2] = consolve(A,b,prob_setup,def_soln,def_setup);
       
    rslts(i-1,6) = t2;
    rslts(i-1,7) = iter;
    
    fprintf('1 V cycle, 1 pre 0 post...\n')
    fprintf('...\n')
    def_soln.gmgpre = 1;
    def_soln.gmgpost = 0;
    [W,iter,t2] = consolve(A,b,prob_setup,def_soln,def_setup);
    
    rslts(i-1,8) = t2;
    rslts(i-1,9) = iter;
    cd PoissonControl
end

cd ..
cd ..
cd Poisson/matlab

rslts
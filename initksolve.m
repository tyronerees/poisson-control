% an m file to initialize the things needed to approx the K solve

switch lower(def_soln.kmethod)
    case 'gmg'
        if def_soln.loadprol == 0
            switch prob_setup.dim
                case 2
                    multdata = genmultdata(def_setup.pow,def_setup.bc);
                case 3
                    multdata = genmultdata3d(def_setup.pow,bdy_set.bdn);
            end
        else
            switch prob_setup.dim
                case 2
                    switch def_setup.bc
                        case 'dirichlet'
                            load multdataQ1_dir
                        case 'neumann'
                            load multdataQ1_neu
                        case 'mixed'
                            load multdataQ1_mix
                    end
                case 3
                    switch def_setup.bc
                        case 'dirichlet'
                            load multdataQ1_dir_3d
%                         case 'neumann'
%                             load multdataQ1_neu_3d
%                         case 'mixed'
%                             load multdataQ1_mix_3d
                    end
            end
        end
        multdata(1).pow = def_setup.pow;
        multdata(1).dim = prob_setup.dim;
        multdata(1).spre = def_soln.gmgpre;
        multdata(1).spost = def_soln.gmgpost;
    case 'amg'
        control = mi20_control;
        control.v_iterations = def_soln.vcyc;
        inform = mi20_setup(A(nu+ny+1:2*ny+nu,nu+1:nu+ny),control);
end
if strcmp(def_soln.kmethod,'gmg') == 0
    multdata = [];
end

if strcmp(prob_setup.eqn,'stokes')==1
    nvel = prob_setup.nvel;
    nvel2 = nvel/2;
    switch def_soln.kmethod_stokes
        case 'amg'
            control = mi20_control;
            control.aggressive = 1;
            control.pre_smoothing = 3;
            control.post_smoothing = 3;
            control.v_iterations = def_soln.vcyc;
            inform = mi20_setup(A(nu+ny+1:ny+nu+nvel2,nu+1:nu+nvel2),control);
    end
end

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
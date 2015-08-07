function multdata = genA(def_setup,multdata,smoother,K)

% generates matrices once for geometric multigrid for the advection
% diffusion equation

for i = def_setup.pow:-1:1            % loop over all the levels
    h = 2^(-i+1);           % define h
    
    %% get the sd constant, delta
    hk = min(abs(h/def_setup.c), abs(h/def_setup.s));
    
    P_h = h*(abs(def_setup.gamma))/((2)*def_setup.epsilon);
    
    % delta (ESW pg 136)
    if P_h > 1
        def_setup.delta = hk*(1-1/P_h)/((2)*abs(def_setup.gamma));
    else
        def_setup.delta = 0;
    end
    
    %% assemble the stiffness matrix
    cd AdvDiffControl
    Abar = stiffsupg_co(h,def_setup,def_setup.delta);
    cd ..
    Abar(:,multdata(i).bdn) =[];
    Abar(multdata(i).bdn,:) = [];
    multdata(i).A1 = Abar;
    
    %% assemble the smoothing matrices
    if multdata(1).smoother == 'gs'
        multdata(i).GSfwd1 = tril(multdata(i).A,0); % (D-L);
        multdata(i).GSbwd1 = triu(multdata(i).A,0); % (D-U);
    end
end

if multdata(1).smoother == 'gs'
        multdata(pow+1).GSfwd1 = tril(K,0); % (D-L);
        multdata(pow+1).GSbwd1 = triu(K,0); % (D-U);
end



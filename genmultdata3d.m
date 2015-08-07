function multdata = genmultdata3d(pow,bdn)
N = 2^pow;
    [multdata(pow).Pro multdata(pow).bdn] = prol3d(N,linspace(0,1,N+1)',bdn);
    
    for level = pow-1:-1:1
       [multdata(level).Pro multdata(level).bdn] = prol3d(2^level,linspace(0,1,2^level+1)',multdata(level+1).bdn);
    end
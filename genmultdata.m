function multdata = genmultdata(pow,bc)
N = 2^pow;
    [multdata(pow).Pro] = prol(N,linspace(0,1,N+1)',bc);
    
    for level = pow-1:-1:1
       [multdata(level).Pro] = prol(2^level,linspace(0,1,2^level+1)',bc);
    end
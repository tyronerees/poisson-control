h = 2^(-4);
M=massmatrix1(h);
b = rand(length(M),1);

i=20;
elts = 'q1';
dim = 2;

[yold] = chebsemiit(M,b,i,dim,elts);

[ynew,resvec] = chebsemiit2(M,b,i,dim,elts);

yact = M\b;

norm(yold-ynew)
norm(ynew-yact)
norm(yold-yact)

semilogy(resvec)

alp = 2/5;

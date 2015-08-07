function [u,ebar,rbar] = mgvcyc_sd2(A,f,multdata,u,level,npre,npost,dim)
% V-cycle for the adjoint problem
M = sparse(diag(diag(A)));
n = length(A);

P = multdata(level).Pro;

Abar = multdata(level).Aadj;

ebar = 1;
J = speye(n)-M\A;
Mf = M\f;

switch dim
    case 2
        omega = 8/9;
    case 3
        omega = 1;
end

switch multdata(1).smoother
    case 'jacobi'
        M = sparse(diag(diag(A)));
        J = speye(n)-M\A;
        Mf = M\f;
        omega = 8/9;

        for j = 1:npre
            u = (1-omega)*u+omega*(J*u + Mf);  
        end
    case 'gs'
        for j =1:npre
            res = f - A*u;
            u = u +  multdata(level+1).GSfwdadj\res;
        end
    case 'blkgs'
        u = blkgs(A,f,u,'fwd',npre);
end


if level == 2
    u = A\f;
else
	rbar = ((1)*P')*(f-A*u);
	ebar = zeros(length(rbar),1);
	ebar = mgvcyc_sd2(Abar,rbar,multdata,ebar,level-1,npre,npost,dim);
	u = u+P*ebar;   
end 

switch multdata(1).smoother
    case 'jacobi'
        for j = 1:npost
            u = (1-omega)*u+omega*(J*u + Mf);  
        end
    case 'gs'
        for j =1:npost
            res = f - A*u;
            u = u +  multdata(level+1).GSbwdadj\res;
        end
    case 'blkgs'
        u = blkgs(A,f,u,'bwd',postits);
end



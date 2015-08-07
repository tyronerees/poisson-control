function [u,ebar,rbar] = mgvcyc(A,f,multdata,u,level,npre,npost,dim)
% M = sparse(diag(diag(A)));
% n = length(A);

P = multdata(level).Pro;

Abar = P'*A*P;
ebar = 1;
% J = speye(n)-M\A;
% Mf = M\f;

switch dim
    case 2
        om = 9/8; % 1/Relax Param
    case 3
        om = 1;   % 1/Relax Param
end

M = om*sparse(diag(diag(A)));

for j = 1:npre
    r = f-A*u;
    u = u + M\r;
%     u = (1-omega)*u+omega*(J*u + Mf);  % k steps? 
end

if level == 1
    u = A\f;
else
	rbar = ((1)*P')*(f-A*u);
	ebar = zeros(length(rbar),1);
	ebar = mgvcyc(Abar,rbar,multdata,ebar,level-1,npre,npost,dim);
	u = u+P*ebar;   
end 

for j = 1:npost
    r = f-A*u;
    u = u + M\r;
% 	u = (1-omega)*u+omega*(J*u + Mf);  % k steps? 
end


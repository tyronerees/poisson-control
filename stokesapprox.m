function x = stokesapprox(A,b,Mp,def_soln,prob_setup)

m = prob_setup.nvel;
n = prob_setup.np;

x0 = zeros(n+m,1);
% S = B*(K\B');
% P = [K sparse(m,n); B -S];

b(m+1:end) = b(m+1:end)-sum(b(m+1:end))*ones(n,1)./n;  %project again


%resvec = zeros(k+1,1);
% resvec(1) = norm(b);

for i = 1:def_soln.Uits;
    r0 = b - A*x0;
%     x1 =  x0+P\r0;
    x1 = x0 + PP(A,Mp,r0,n,m,def_soln,prob_setup);
    x0 = x1;    
%     % correct pressure
%     xcor = x1;
%     xcor(m+1:end) = xcor(m+1:end)-sum(xcor(m+1:end))*ones(n,1)./n;
%     resvec(i+1) = norm(b - [K B'; B sparse(n,n)]*xcor);
end

u = x1(1:m);
p = x1(m+1:end)-sum(x1(m+1:end))*ones(n,1)./n; % project into right space
x = [u;p];
% figure(14)
% semilogy(resvec)
% title('Residual')
% pause

%%
function x = PP(A,M,r,n,m,def_soln,prob_setup)

if strcmp(def_soln.conmethod,'pblktri') ==1
    switch def_soln.kmethod_stokes
        case 'bslash'
            xa = A(1:m,1:m)\r(1:m);
        case 'amg'
            xa1 = mi20_precondition(r(1:m/2));
            xa2 = mi20_precondition(r(m/2+1:m));
            xa = [xa1;xa2];
           % xa = mi20_precondition(r(1:m));
    end
    
    switch def_soln.mmethod_stokes
        case 'bslash'
            xb = -M\(r(m+1:end)-A(m+1:end,1:m)*xa);
        case 'chebit'
            xb = -chebsemiit2(M,r(m+1:end)-A(m+1:end,1:m)*xa,def_soln.sprecheb,prob_setup.dim,prob_setup.pelt);
    end
elseif strcmp(def_soln.conmethod,'diag') ==1
    xa  = A(1:m,1:m)\(r(1:m));
    xb = M\r(m+1:end);
end
x = [xa;xb];
%         zed5 = Mp\zed5b*(1/delta);
%         zed4 = K\(zed4b - B'*zed5);
%u = A(1:prob_setup.nvel,1:prob_setup.nvel)\b(1:prob_setup.nvel);
%p = -Mp\(A(prob_setup.nvel+1:end,1:prob_setup.nvel)*u - b(prob_setup.nvel+1:end));

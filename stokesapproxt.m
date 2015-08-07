function x = stokesapproxt(A,b,Mp,def_soln,prob_setup)

m = prob_setup.nvel;
n = prob_setup.np;

x0 = zeros(n+m,1);

% S = B*(K\B');
% [n,m] = size(B);
% P = [K B'; sparse(n,m) -S];

b(m+1:end) = b(m+1:end)-sum(b(m+1:end))*ones(n,1)./n;  %project again


%resvec = zeros(k+1,1);
%resvec(1) = norm(b);

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
    switch def_soln.mmethod_stokes
        case 'bslash'
            xb = -M\r(m+1:end);
        case 'chebit'
            xb = -chebsemiit2(M,r(m+1:end),def_soln.sprecheb,prob_setup.dim,prob_setup.pelt);
    end
    
    switch def_soln.kmethod_stokes
        case 'bslash'
            xa = A(1:m,1:m)\(r(1:m) - A(1:m,m+1:end)*xb);
        case 'amg'
            rhs = r(1:m) - A(1:m,m+1:end)*xb;
            xa1 = mi20_precondition(rhs(1:m/2));
            xa2 = mi20_precondition(rhs(m/2+1:m));
            xa = [xa1;xa2];
            % xa = mi20_precondition(rhs(1:m));
    end
elseif strcmp(def_soln.conmethod,'diag') ==1
    xa  = A(1:m,1:m)\(r(1:m));
    xb = M\r(m+1:end);
end
        
x = [xa;xb];


%p = Mp\b(prob_setup.nvel+1:end);
%u = A(1:prob_setup.nvel,1:prob_setup.nvel)\(b(1:prob_setup.nvel) -A(1:prob_setup.nvel,prob_setup.nvel+1:end)*p);

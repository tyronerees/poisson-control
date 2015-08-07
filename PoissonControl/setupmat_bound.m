% setupmat.m
%
% The function 16*(x-0.5)^2*(y-0.5)^2
% Boundary values set in the matrix

function [B,K,Mu,My,Mhat,ubdy,uhat] = setupmat_bound(h,bdy_set,objective,plots)

N=1/h;
% objective = 1;          % 1 for biquadratic peak

% method = 2;               % 1 for insert zeros, 2 for delete rows

m=(N+1)^2;
[x y quad bdy] = mesh_bdy(N); 

K = stiffmatrixbdy3(h);
My= massmatrix1(h);
C = bdymassmatrix(h);

C2 = sparse(m,m);
Ck = h/(6)*[2 1; 1 2];

for k=1:length(bdy(:,1))
    if (bdy(k,2)-bdy(k,1)) >= 1
        C2(bdy(k,1:2),bdy(k,1:2))=C2(bdy(k,1:2),bdy(k,1:2))+Ck(1:2,1:2); 
    elseif (bdy(k,2) - bdy(k,1)) <= 1
        C2(bdy(k,1:2),bdy(k,1:2))=C2(bdy(k,1:2),bdy(k,1:2))+Ck(1:2,1:2); 
    end
end

C(:,setxor(bdy_set.bdn,1:m)) = [];        % restrict to just the elements on the bdy
                                  % in the one dimension
Mhat = C;

% define the function u hat
switch objective
    case 1
        uhat = 16*(x-0.5).^2.*(y-0.5).^2;

        uhat((x>0.5)|(y>0.5))=0;

        uh = @(x,y) 16*(x-0.5).^2*(y-0.5).^2;                         %%%%%%%%%%
    
    case 2
        sig = 0.125;
%         uhat = (pi/sig^2)*exp(-((x-0.5).^2+(y-0.5).^2)/sig^2);
        uhat = exp(-((x-0.5).^2+(y-0.5).^2)/sig^2);
    case 3
        uhat = ones(1,length(x));
        uhat((x>0.5)|(y>0.5))=0;
    case 4
        uhat = x+y;
        uh = @(x,y) x+y;
    case 5
        uhat = exp(x+y);
        uh = @(x,y) exp(x+y);
    case 6
        sig = 0.125;
        rhs = exp(-((x-0.5).^2+(y-0.5).^2)/sig^2);
%         rhs = ones(length(M),1);
        uhat = FPoissSolve(rhs',h,bdy_set.bdn,zeros(4*N,1));
        uhat = uhat';
    case 7
        rhs = zeros(length(M),1);
%         bv = @(x,y) x+y;
        bv = @(x,y) exp(x^2+y^2);
        uhatbv = zeros(4*N,1);
        for i = 1:length(bdy_set.bdn)
            uhatbv(i)=bv(x(bdy_set.bdn(i)),y(bdy_set.bdn(i)));
        end
        uhat = FPoissSolve(rhs,h,bdy_set.bdn,uhatbv');
        uhat = uhat';
end
        

if plots == 1
    puhat = reshape(uhat,N+1,N+1);
    x1=0:h:1;
    y1=x1;
    figure(2)
    subplot(2,2,1)
    surf(x1,y1,puhat);
    title('Plot of \^u');

end

% define the rhs vector B = [a;b;c] ( b_i = \int( \hat{u} \phi_i ) )

a = zeros(4*N,1);
b = My*uhat';

f = zeros(m,1);               % define the Neumann bc
c = My*f;


% get the values of u on the boundary

ubdy = zeros(4*N,1);
switch objective
    case 1
        for i = 1:length(bdy_set.bdn)
            if x(bdy_set.bdn(i))<0.5 && y(bdy_set.bdn(i))<0.5
               ubdy(i)=uh(x(bdy_set.bdn(i)),y(bdy_set.bdn(i)));
            end
        end
    case 3
        for i = 1:length(bdy_set.bdn)
            if x(bdy_set.bdn(i))<0.5 && y(bdy_set.bdn(i))<0.5
               ubdy(i)=1;
            end
        end   
    case 4
        for i = 1:length(bdy_set.bdn)
            ubdy(i)=uh(x(bdy_set.bdn(i)),y(bdy_set.bdn(i)));
        end
    case 5
        for i = 1:length(bdy_set.bdn)
            ubdy(i)=uh(x(bdy_set.bdn(i)),y(bdy_set.bdn(i)));
        end
    case 7
        ubdy = uhat(bdy_set.bdn)';
end

ubdy = ubdy(bdy_set.dbynodes);             % pick out the values that are needed (better way to do this?!)

Mu = C;
Mu(setxor(bdy_set.bdn,1:m),:) = [];   



% create right hand side

B = [a;b;c];
 

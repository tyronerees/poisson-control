% setupmat.m
%
% The function 16*(x-0.5)^2*(y-0.5)^2
% Boundary values set in the matrix

function [B,K,M,ubdy,uhat] = setupmat(h,bdy_set,objective,plots)

N=1/h;
% objective = 1;          % 1 for biquadratic peak

% method = 2;               % 1 for insert zeros, 2 for delete rows


[x y quad] = mesh_bdy(N); 

K = stiffmatrixbdy3(h);
M=massmatrix1(h);

m=(N+1)^2;

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
if plots == 2
    puhat = reshape(uhat,N+1,N+1);
    x1=0:h:1;
    y1=x1;
    figure(1)
    surf(x1,y1,puhat);
end

% define the rhs vector B = [a;b;c] ( b_i = \int( \hat{u} \phi_i ) )

a = zeros(m,1);
% canonical basis functions

if objective == (1||3)
    b = rhsb(h,objective);
else
    b = M*uhat';
end


c = zeros(m,1);


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


% update right hand side

c = c-K(:,bdy_set.dirn)*ubdy;           
b = b;%M(:,bdy_set.dirn)*ubdy;


% switch method
%     case 1
%         % set boundary values
%         b(bdn) = 0;
%         c(bdn) = ubdy;
% 
%         % Put zeros on bdy of M
%         M2 = M;
% 
%         M2(:,bdn) = 0;
%         M2(bdn,:) = 0;
% 
%         % Put ones on diag of K
%         K(:,bdn) = 0;
%         K(bdn,:) = 0;
%         K(bdn,bdn) = speye(length(ubdy));
%     case 2
        % set boundary values
        a(bdy_set.dirn) = [];
        b(bdy_set.dirn) = [];
        c(bdy_set.dirn) = [];
        
       	% delete bdy of M
        M2 = M;

        M2(:,bdy_set.dirn) = [];
        M2(bdy_set.dirn,:) = [];
        M = M2;
        
        % Delete bdy of K
        K(:,bdy_set.dirn) = [];
        K(bdy_set.dirn,:) = [];
        
% end

% create right hand side

B = [a;b;c];
 

function u = FPoissSolve(f,h,bdn,ubdy,subfig)
if nargin == 4
    subfig =1;
end

% create matrices
M = massmatrix1(h);
K = stiffmatrixbdy3(h);


% get right hand side
F = M*f';

% apply bcs
F(bdn) = ubdy;
K(bdn,:) = 0;
K(bdn,bdn) = speye(length(bdn));

u = K\F;

U = reshape(u,1/h+1,1/h+1);

x=0:h:1;
y=x;

if subfig==2
    figure(4)
else
    subplot(2,2,4)
end
surf(x,y,full(U));

if subfig==1
title('Plot of soln to forward prob');
end
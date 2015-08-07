function [y2] = chebsemiit(A,b,i,dim,elts)

% Generates the ith chebyshev semi iteration vector using Jacobi iteration.
%
% Input : A - Matrix to be solved
%       : b - Right hand side
%       : i - number of iterations of the Cheyshev semi iteration
%
% Output : y2 - approximation to y satisfying A*y = b;

%dim = 2;                     % Switch for the dimesion of the problem 

% tic,
l = length(A);

I = speye(l);                     % define identity
invM = spdiags((1./diag(A)),0,l,l);   % Compute inverse once



% w is the relaxation parameter
% mu is the largest eig for this relaxation parameter
% using the mass matrix as A (see notes 30-11-07)

switch dim
    case 2
        if strcmp(elts,'q1') == 1 
            w = 4/5;
            mu = 5/4; 
        elseif strcmp(elts,'minres') == 1
                                       
        end
    case 3
        if strcmp(elts,'q1') == 1 
            w = 4/7;
            mu = 14/13;   
        elseif strcmp(elts,'minres') == 1
        end
end

g = w*invM*b;                           % define relaxed Jacobi g
S = (I-w*invM*A);                       % define the iteration Matrix

y1 = zeros(length(b),1);
y2 = S*y1 + g;                          % x_(m+1) = Sx_m + g

ck0 = 2;
ck1 = mu;

for j = 2:i
   ck2 = mu*ck1-0.25*ck0;

     w = mu*ck1/ck2;
    y3 = w*(S*y2+g-y1)+y1;  
    
    y1 = y2;    y2 = y3;
   ck0 = ck1;
   ck1 = ck2;
end
% t = toc;

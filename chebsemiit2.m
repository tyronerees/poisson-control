function [y] = chebsemiit2(A,b,i,dim,elts)

% Generates the ith chebyshev semi iteration vector using Jacobi iteration.
%
% Input : A - Matrix to be solved
%       : b - Right hand side
%       : i - number of iterations of the Cheyshev semi iteration
%       : dim - dimension of problem
%       : elts - elements used
%
% Output : y - approximation to y_e satisfying A*y_e = b;
%
% T. Rees, December 2010

y0 = zeros(length(b),1);
y00 = y0;
w = 1;


switch dim
    case 2
        if strcmp(elts,'q1') == 1 
            lmin = 1/4; % Min eig
            lmax = 9/4; % Max eig
        elseif strcmp(elts,'q2') == 1
            lmin = 1/4;
            lmax = 25/16;
        end
    case 3
        if strcmp(elts,'q1') == 1 
            lmin = 1/8;
            lmax = 27/8;
        elseif strcmp(elts,'q2') == 1
            lmin = 0.3924;  
            lmax = 2.0598; 
        end
end

alp = (lmin+lmax)/2; %1/Relaxation parameter           
rho = (lmax-lmin)/(lmax+lmin); % Max eig
Mdiag=alp*diag(A); 

for k = 1:i
   w = 1/(1-(w*rho^2)/4);
   r = b - A*y0;
   z = Mdiag.\r;
   y = w*(z + y0 - y00) + y00;
   y00 = y0; y0 = y;
end
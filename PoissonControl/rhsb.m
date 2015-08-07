% stiffmatrixbdy.m
%
% The stiffness matrix for the standard (Laplace) problem including the
% boundary

function [b] = rhsb(h,obj)
N=1/h;
NN = (N+1)^2;

b = zeros(NN,1);
%Ko = sparse(NN,NN);

[x y quad] = mesh2d(N);


% Global Load vector Matrix
if obj == 1
    % Matrix CHI, each row is CHI_i evaluated at quad points +/- \sqrt(3) 
    CHI = 1/12*[(1-sqrt(3))^2 2 2 (1+sqrt(3))^2;
        2 (1+sqrt(3))^2 2 (1-sqrt(3))^2;
        (1+sqrt(3))^2 2 2 (1-sqrt(3))^2;
        2 (1-sqrt(3))^2 2 (1+sqrt(3))^2];
    for k=1:length(quad(:,1))
        qq = 1/sqrt(3);
        UHAT = [(2*x(quad(k,1)) + h/2*(qq+1)-1)^2*(2*y(quad(k,1)) + h/2*(qq+1)-1)^2;
            (2*x(quad(k,1)) + h/2*(qq+1)-1)^2*(2*y(quad(k,1)) + h/2*(-qq+1)-1)^2;
            (2*x(quad(k,1)) + h/2*(-qq+1)-1)^2*(2*y(quad(k,1)) + h/2*(-qq+1)-1)^2;
            (2*x(quad(k,1)) + h/2*(-qq+1)-1)^2*(2*y(quad(k,1)) + h/2*(qq+1)-1)^2]; 
        if (x(quad(k,1))< 0.5)&&(y(quad(k,1))<0.5)
            bk = (h^2/4)*CHI*UHAT;
        else
            bk = [0; 0; 0; 0];
        end
    b(quad(k,1:4))=b(quad(k,1:4))+bk(1:4); 
    end
end

if obj == 3
    for k=1:length(quad(:,1))
        if (x(quad(k,1))< 0.5)&&(y(quad(k,1))<0.5)
            bk = (h^2/4)*[1; 1 ;1; 1];
        else
            bk = [0; 0; 0; 0];
        end
    b(quad(k,1:4))=b(quad(k,1:4))+bk(1:4); 
    end
end

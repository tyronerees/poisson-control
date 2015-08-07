% stiffmatrixbdy.m
%
% The stiffness matrix for the standard (Laplace) problem including the
% boundary

function [b] = rhsb3d(h,obj)
N=1/h;
NN = (N+1)^3;

b = zeros(NN,1);
%Ko = sparse(NN,NN);

[x y z cub bdy] = mesh3d_bdy(N);

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
    for k=1:length(cub(:,1))
        if (x(cub(k,1))< 0.5)&&(y(cub(k,1))<0.5&&(z(cub(k,1))<0.5))
            bk = (h^2/4)*ones(8,1);
        else
            bk = zeros(8,1);
        end
    b(cub(k,1:8))=b(cub(k,1:8))+bk(1:8); 
    end
end

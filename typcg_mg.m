function [x,i] = typcg_mg(M,B,multdata,pow,dim)

r=B; 
nr0 = norm(r);
if nr0==0   % added to stop division by zero
    x = zeros(length(B),1);   %
else
    z=mgvcyc(M,r,multdata,zeros(size(B)),pow,1,1,dim);
    p=z;
    w=M*p; 
    gam=p'*w;
    del=z'*r;
    alp=del/gam; 
    iterp=1;
    x=alp*p;
    r=r-alp*w;
    nr = norm(r);
    i = 1;
    
    while nr>1e-6*nr0 && i<100
         z=mgvcyc(M,r,multdata,zeros(size(B)),pow,2,2,dim);
         alp=z'*r;
%             if alp<1.e-15, 
%                 return, 
%             end
        bet=alp/del;
        p=z+bet*p;
        w=M*p;
        gam=p'*w;
        del=alp;
        alp=alp/gam;
        iterp=iterp+1;
        x=x+alp*p;
        r=r-alp*w;
        nr = norm(r);
        i = i+1;
    end
end

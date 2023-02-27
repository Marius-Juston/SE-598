function Ex02Sol_P1()
clear all; clc
% define the problem for x, g, dg/dx
ux = [3.5, 3.8]; stdx = [0.5, 0.5];
gx=@(x1, x2) 1- (x1+x2-5).^2/30-(x1-x2-12).^2/120;
dg1 = @(x1, x2) - (x1+x2-5)/15-(x1-x2-12)/60;
dg2 = @(x1, x2) - (x1+x2-5)/15+(x1-x2-12)/60;

% MCS method with N = 1000000 for PoF
Nmcs = 1000000;
xm1 = normrnd(ux(1), stdx(1), Nmcs, 1);
xm2 = normrnd(ux(2), stdx(2), Nmcs, 1);
gmcs = gx(xm1, xm2); 
PoF_mcs = sum(gmcs>0)/Nmcs


% FORM method with HL_RF for PoF
u=zeros(1,2); iter=0;  Dif=1; sign = 1;        
while Dif >= 1d-5 && iter < 20
    iter=iter + 1;
    x = ux + u.*stdx;
    g = gx(x(1), x(2)); 
    DG=[dg1(x(1), x(2)), dg2(x(1), x(2))].*stdx;                       
    u=(DG*u'-g)/norm(DG)^2*DG;
    U(iter,:)=u/norm(u);
    if iter ==1
        sign = -g/abs(g);
    elseif iter>1
        Dif=abs(U(iter-1,:)*U(iter,:)' - 1);
    end
end
beta = sign*norm(u);
PoF_HL_RF = normcdf(-beta)


function Ex02Sol_P2()
clear all; clc
% define the problem for x, g
w =2; t = 4; L = 100; E = 2.9e7; S0 = 35000; D0 = 2.5;
ux = [500, 1000]; stdx = [100, 100];

sqroot = @(Px, Py) ((Px/w^2).^2+(Py / t^2).^2).^(0.5);
gx =@(Px, Py) (4 *L^3) *  sqroot(Px, Py) / (E * t * w)- D0;

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
    [g,DG]=cons(u,ux,stdx);                     
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


function [g,DG]=cons(u, ux, stdx)
xtemp = u.* stdx +ux;
Px = xtemp(1); Py = xtemp(2);
w =2; t = 4; L = 100; E = 2.9e7; S0 = 35000; D0 = 2.5;

sqroot = ((Px/w^2)^2+(Py / t^2)^2)^(0.5);
g=(4 *L^3) *  sqroot / (E * t * w) - D0;
DG(1) = 4 * L^3 * Px / ( E * t * w^5 * sqroot)*stdx(1);
DG(2) = 4 * L^3 * Py / (E * t^5 * w * sqroot)*stdx(2);

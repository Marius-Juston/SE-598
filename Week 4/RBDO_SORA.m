%%%%%%%%%%  RBDO with FORM using SORA  %%%%%%%%
function RBDO_SORA()
clear all; close all; clc;
global nc nd bt stdx mpp0
nc=3;  nd=2; bt=norminv(0.99,0,1); stdx=[0.6,0.6];
mpp0 = [0, 0; 0 0; 0 0];
x0=[5,5]; lb=[0,0]; ub=[10,10];
options = optimset('GradConstr','on','GradObj','on','LargeScale','off');
dif = 1; ki = 0;
while dif>1e-5 && ki <50
    ki = ki+1;
    [x, ~]=fmincon(@Costfun,x0,[],[],[],[],lb,ub,@nonlcon,options)
    MPP(1,:) = x-AMV(x,1);
    MPP(2,:) = x-AMV(x,2);
    MPP(3,:) = x-AMV(x,3);
    dif = sum(sum((mpp0-MPP).^2));
    MPP_Val{ki} =mpp0;
    mpp0 = MPP
end

%====================      Obj. Function   ===============================%
function [f,g]= Costfun(x)
f=x(1)+x(2);
g=[1 1];

%============================ Constraint Fun. ============================%
function [c,ceq,GC,GCeq]=nonlcon(x)
global mpp0
ceq=[]; GCeq=[];

% for constraint 1/2/3
x1 = x-mpp0(1,:); 
x2 = x-mpp0(2,:);
x3 = x-mpp0(3,:);
[c1,GC1]=cons(x1,1);
[c2,GC2]=cons(x2,2);
[c3,GC3]=cons(x3,3);
c = [c1, c2, c3];
GC = [GC1, GC2, GC3]; 

%================== RIA Approach with AMV Algorithm ====================%
function [xmpp]=AMV(x,kc)
global bt stdx nd
u=zeros(1,nd); iter = 0; Dif=1;
while Dif>1d-5 & iter<20
    iter=iter+1;
    if iter>1
        u=DG*bt/norm(DG);
    end
    [~,DG]=cons(x+u.*stdx,kc);
    U(iter,:)=u/bt;
    if iter>1
        Dif=abs(U(iter,:)*U(iter-1,:)'-1);
    end
end
xmpp = x+u'.*stdx;

function [c,GC]=cons(x,kc)
if kc == 1
    c=1-x(1)^2*x(2)/20;
    GC(1,1)=-x(1)*x(2)/10;
    GC(2,1)=-x(1)^2/20;
elseif kc == 2
    c=1-(x(1)+x(2)-5)^2/30-(x(1)-x(2)-12)^2/120;
    GC(1,1)=(-(x(1)+x(2)-5)/15-(x(1)-x(2)-12)/60);
    GC(2,1)=(-(x(1)+x(2)-5)/15+(x(1)-x(2)-12)/60);
elseif kc == 3
    c=1-80/(x(1)^2+8*x(2)+5);
    GC(1,1)=x(1)*160/((x(1)^2+8*x(2)+5))^2;
    GC(2,1)=80*8/((x(1)^2+8*x(2)+5))^2;
end
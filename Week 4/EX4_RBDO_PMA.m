%%%%%%%%%%  RBDO with FORM using PMA Approach  %%%%%%%%
function RBDO_PMA()
clear; close all; clc;
global nc nd bt stdx Iters Cost G1Constraints G2Constraints Pmu Pstdx np

G1Constraints = [];
G2Constraints = [];

nc=2;  nd=2;  bt=norminv(0.99,0,1); 
lb=[0,0]; ub=[10,10]; 

x0=[4,  6];                                % intial design point
stdx=[0.2,0.2];                        % standard deviation of the random variable

Pmu = [500, 1000];
Pstdx = [100, 100];
np = 2;

xp=x0; Iters=0;
options = optimset('GradConstr','on','GradObj','on','LargeScale','off');
[x_opt,~] = fmincon(@Costfun,x0,[],[],[],[],lb,ub,@frelcon,options)

disp(G1Constraints);
disp(G2Constraints);
 
size(G1Constraints)
size(G2Constraints)

%====================      Obj. Function   ===============================%
function [f,g]= Costfun(x)
    f=x(1)*x(2);
    g=[x(2) x(1)];
    Cost=f;
end
%====================  Define Constraints and Gradiants  =================%
function [c,ceq,GC,GCeq] = frelcon(x) 
    ceq=[]; GCeq=[];
    for j = 1:nc
            [G,DG]=AMV(x,j);
            beta(j)=G;
            dbeta(:,j)=DG./stdx;
    end

    c=beta; GC=dbeta; 
    dx=norm(x-xp);
    if  dx>1d-5  || Iters == 0
        Iters=Iters+1;
        SHOW(Iters,x,c,GC);
    end
    xp = x;
end
%=================  PMA Approach with AMV Algorithm  =====================%
function [G,DG]=AMV(x,kc)
    u=zeros(1,nd + np); iter = 0; Dif=1;
    while Dif>1d-5 && iter<20                      
        iter=iter+1; 
        if iter>1
            u=DG*bt/norm(DG);
        end        
        [G,DG]=cons(u,x,kc);
        U(iter,:)=u/bt;        
        if iter>1
            Dif=abs(U(iter,:)*U(iter-1,:)'-1);
        end
    end

    DG = DG(1:2);
end
%============================ Constraint Fun. ============================%
function [g,DG]=cons(u,d,kc)
    x = u(1:2).*stdx+d;
    xP = u(3:4).*Pstdx + Pmu;

    w = x(1);
    t = x(2);
    Px = xP(1);
    Py = xP(2);

    L = 100;
    E = 2.9e7;
    S0 = 35000;
    D0 = 2.5;

    if kc == 1
        g= (6 * L) * (Px / w + Py / t) / (t * w) - S0;
        DG(1)=-(6 * L * Px) / (t * w^3) - ( (6 * L) * (Py / t + Px / w))/ (t * w^2);
        DG(2)=-(6 * L * Py) / (t^3 * w) - ( (6 * L) * (Py / t + Px / w))/ (t^2 * w);
        
        DG(3) = (6 * L) / (t * w^2);
        DG(4) = (6 *L) / (t^2 * w);

        G1Constraints = [G1Constraints; x, g, DG];
    elseif kc == 2
        sqr = sqrt((Py / t^2)^2 + (Px/w^2)^2);

        g=(4 *L^3) *  sqr / (E * t * w)- D0;
        DG(1)=-(8 * L^3 * Px^2)/(E * t * sqr * w^6) - (4 *L^3 * sqr)/(E *t * w^2);
        DG(2)=-(8 * L^3 * Py^2)/(E * t^6 * sqr * w) - (4 *L^3 * sqr)/(E *t^2 * w);
        
        DG(3) = 4 * L^3 * Px / ( E * t * w^5 * sqr);
        DG(4) = 4 * L^3 * Py / (E * t^5 * w * sqr);

        G2Constraints = [G2Constraints; x, g, DG];
    end

    DG(1:2) = DG(1:2).*stdx;
    DG(3:4) = DG(3:4).*Pstdx;
end
function  SHOW(Iters,x,c,GC)%====== Display the Iteration Information=====%
    fprintf(1,'\n********** Iter.%d ***********\n' ,Iters);
    disp(['Des.: ' sprintf('%6.4f  ',x)]);
    disp(['Obj.: ' sprintf('%6.4f',Cost)]);  
    disp(['Cons.: ' sprintf('%6.4f  ',c)]);
    disp(['Sens.: ' sprintf('%6.4f  ',GC)]);
    fprintf('\n\n')
end
end

%%%%%%%%%%  RBDO with FORM using SORA  %%%%%%%%
function EX4_RBDO_SORA()
    clear all; close all; clc;
    global nc nd bt stdx mpp0 Pmu Pstdx np
    nc=2;  nd=2; bt=norminv(0.99,0,1); stdx=[0.2,0.2];

    Pmu = [500, 1000];
    Pstdx = [100, 100];
    np = 2;

    mpp0 = [0 0; 0 0];
    x0=[4,6]; lb=[0,0]; ub=[10,10];
    options = optimset('GradConstr','on','GradObj','on','LargeScale','off');
    dif = 1; ki = 0;
    while dif>1e-5 && ki <50
        ki = ki+1;
        [x, ~] = fmincon(@Costfun,x0,[],[],[],[],lb,ub,@nonlcon,options)
        MPP(1,:) = x-AMV(x,1);
        MPP(2,:) = x-AMV(x,2);
        dif = sum(sum((mpp0-MPP).^2));
        MPP_Val{ki} =mpp0;
        mpp0 = MPP;
    end
end

%====================      Obj. Function   ===============================%
function [f,g]= Costfun(x)
    f=x(1)*x(2);
    g=[x(2) x(1)];
end

%============================ Constraint Fun. ============================%
function [c,ceq,GC,GCeq]=nonlcon(x)
    global mpp0 Pmu
    ceq=[]; GCeq=[];
    
    % for constraint 1/2/3
    x1 = x-mpp0(1,:); 
    x2 = x-mpp0(2,:);

    x1 = [x1 Pmu];
    x2 = [x2 Pmu];

    [c1,GC1]=cons(x1,1);
    [c2,GC2]=cons(x2,2);

    GC1 = GC1(1:2);
    GC2 = GC2(1:2);

    c = [c1, c2];
    GC = [GC1', GC2']; 
end

%================== RIA Approach with AMV Algorithm ====================%
function [xmpp]=AMV(x,kc)
    global bt stdx nd np Pstdx Pmu
    u=zeros(1,nd + np); iter = 0; Dif=1;
    while Dif>1d-5 && iter<20
        iter=iter+1;
        if iter>1
            u=DG*bt/norm(DG);
        end
         xi = u.*[stdx, Pstdx]+ [x, Pmu];

        [~,DG]=cons(xi, kc);
        U(iter,:)=u/bt;
        if iter>1
            Dif=abs(U(iter,:)*U(iter-1,:)'-1);
        end
    end
    u = u(1:2);
    xmpp = x+u.*stdx;
end

function [c,GC]=cons(x,kc)
    w = x(1);
    t = x(2);
    Px = x(3);
    Py = x(4);

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
    elseif kc == 2
        sqr = sqrt((Py / t^2)^2 + (Px/w^2)^2);

        g=(4 *L^3) *  sqr / (E * t * w)- D0;
        DG(1)=-(8 * L^3 * Px^2)/(E * t * sqr * w^6) - (4 *L^3 * sqr)/(E *t * w^2);
        DG(2)=-(8 * L^3 * Py^2)/(E * t^6 * sqr * w) - (4 *L^3 * sqr)/(E *t^2 * w);
        
        DG(3) = 4 * L^3 * Px / ( E * t * w^5 * sqr);
        DG(4) = 4 * L^3 * Py / (E * t^5 * w * sqr);
    end

    c = g;
    GC = DG;
end
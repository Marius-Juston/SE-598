% RBDO CODE with System PoF constraint using FORM_HLRF method
function EX5_RBDO_SysRel()
close all; 
clear; 
clc 
global nc nd Rt_sys stdx Iters Cost Pmu Pstdx np    % define global variables that can be used in different subroutines 
nc=2;                                     % # of constraints   
nd=2;                                     % # of design variables 
Rt_sys=0.90;                            % target system reliability
x0=[4,  6];                                % intial design point
stdx=[0.2,0.2];                        % standard deviation of the random variable
lb=[0,0]; ub=[10,10];              % low bound and upper bound

Pmu = [500, 1000];
Pstdx = [100, 100];
np = 2;


xp=x0;                                    % xp is used to store the design point of previous iteration - for convergence check
Iters=0;                                   % iteration index
options = optimoptions('fmincon','Display', 'iter-detailed', 'SpecifyConstraintGradient', false, 'SpecifyObjectiveGradient',true);
[xopt,~]=fmincon(@Costfun,x0,[],[],[],[],lb,ub,@frelcon,options)
%====================      Obj. Function   ===============================%
    function [f, g]= Costfun(x)
        f=x(1) * x(2);
        g = [x(2) x(1)];
        Cost=f;
    end
%==========  Define Outer Loop Constraints and Gradiants  =================%
    function [c, ceq] = frelcon(x)
        ceq=[];
        for j = 1:nc
            [G]=HL_RF(x,j);
            beta(j)=G;
        end
        Rel = normcdf(beta, 0,1);
        PoF_sys =  sum(1-Rel)-(1-Rel(1))*(1-Rel(2)); 
        c=Rt_sys-(1-PoF_sys);
        dx=norm(x-xp);
%         if  dx>1d-5  || Iters == 0
%             Iters=Iters+1;
%             SHOW(Iters,x,c);
%         end
        xp = x;
    end

%======== Inner Loop using FORM with HL_RF Algorithm ==============%
    function [beta]=HL_RF(x,kc)
        u=zeros(1,nd + np); iter=0;  Dif=1; sign = 1;
        while Dif >= 1d-5 && iter < 20
            iter=iter + 1;
            [g,DG]=cons(u,x,kc);
            u=(DG*u'-g)/norm(DG)^2*DG;
            U(iter,:)=u/norm(u);
            if iter ==1
                sign = -g/abs(g);
            elseif iter>1
                Dif=abs(U(iter-1,:)*U(iter,:)' - 1);
            end
        end
        beta = sign*norm(u);
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
        g=6 * L * (Px * t + Py * w) / (S0 * t^2 * w^2) - 1;
        DG(1)=-6 * L * (2 * Px * t + Py * w) / (S0 * t^2 * w^3);
        DG(2)=-6 * L * (Px * t + 2 * Py * w) / (S0 * t^3 * w^2);
        
        DG(3) = 6 * L / (S0 * t * w^2);
        DG(4) = 6 *L / (S0 * t^2 * w);
    elseif kc == 2
        g=4 *L^3 * sqrt(Py^2 / t^4 + Px^2/w^4) / (D0 * E * t * w)-1;
        DG(1)=-4 * L^3 * (Px^2 *t^4 + 3 * Py^2 * w^4)/ (D0 * E * t^5 * w^6 * sqrt(Py^2/t^4 + Px^2/w^4));
        DG(2)=-4 * L^3 * (3 * Px^2 * t^4 + Py^2 * w^4) / (D0 * E * t^6 * w^5 * sqrt(Py^2/t^4 + Px^2/w^4));
        
        DG(3) = 4 * L^3 * Px / (D0 * E * t * w^5 * sqrt(Py^2/t^4 + Px^2/w^4));
        DG(4) = 4 * L^3 * Py / (D0 *E * t^5 * w * sqrt(Py^2 / t^4 + Px^2 / w^4));
    end

    DG(1:2) = DG(1:2).*stdx;
    DG(3:4) = DG(3:4).*Pstdx;
end
end

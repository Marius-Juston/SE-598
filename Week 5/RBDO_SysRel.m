% RBDO CODE with System PoF constraint using FORM_HLRF method
function RBDO_SysRel()
close all; 
clear; 
clc 
global nc nd Rt_sys stdx Iters Cost    % define global variables that can be used in different subroutines 
nc=3;                                     % # of constraints   
nd=2;                                     % # of design variables 
Rt_sys=0.90;                            % target system reliability
x0=[5,5];                                % intial design point
stdx=[0.6,0.6];                        % standard deviation of the random variable
lb=[0,0]; ub=[10,10];              % low bound and upper bound
xp=x0;                                    % xp is used to store the design point of previous iteration - for convergence check
Iters=0;                                   % iteration index
options = optimoptions('fmincon','Display', 'iter-detailed', 'SpecifyConstraintGradient', false, 'SpecifyObjectiveGradient',true);
[xopt,~]=fmincon(@Costfun,x0,[],[],[],[],lb,ub,@frelcon,options)
%====================      Obj. Function   ===============================%
    function [f, g]= Costfun(x)
        f=x(1)+x(2);
        g = [1 1];
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
        PoF_sys =  sum(1-Rel)-(1-Rel(1))*(1-Rel(2))-max((1-Rel(1))*(1-Rel(3)),(1-Rel(2))*(1-Rel(3))); 
        c=Rt_sys-(1-PoF_sys);
        dx=norm(x-xp);
%         if  dx>1d-5  || Iters == 0
%             Iters=Iters+1;
%             SHOW(Iters,x,c);
%         end
        xp = x;
    end

%======== Inner Loop using FORM with HL_RF Algorithm ==============%
    function [beta,dbeta]=HL_RF(x,kc)
        u=zeros(1,nd); iter=0;  Dif=1; sign = 1;
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
        x = u.*stdx+d;
        if kc == 1
            g=1-x(1)^2*x(2)/20;
            DG(1)=-x(1)*x(2)/10*stdx(1);
            DG(2)=-x(1)^2/20*stdx(2);
        elseif kc == 2
            g=1-(x(1)+x(2)-5)^2/30-(x(1)-x(2)-12)^2/120;
            DG(1)=(-(x(1)+x(2)-5)/15-(x(1)-x(2)-12)/60)*stdx(1);
            DG(2)=(-(x(1)+x(2)-5)/15+(x(1)-x(2)-12)/60)*stdx(2);
        elseif kc == 3
            g=1-80/(x(1)^2+8*x(2)+5);
            DG(1)=x(1)*160*stdx(1)/((x(1)^2+8*x(2)+5))^2;
            DG(2)=80*8*stdx(2)/((x(1)^2+8*x(2)+5))^2;
        end
    end
end

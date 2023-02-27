% RBDO CODE with Double Loop formulation using FORM_HLRF method
function Ex03SOL_RBDO_HLRF()
close all; 
clear; 
clc 
global nc nd bt stdx uP stdP Iters Cost x_all_f x_all_g1 x_all_g2  % define global variables that can be used in different subroutines 
nc=2;                                     % # of constraints   
nd=2;                                     % # of design variables 
bt=norminv(0.99,0,1);               % target reliability index
x0=[4,6];                                 % intial design point
stdx=[0.2,0.2];                        % standard deviation of the design variable
uP = [500, 1000];                     % mean value of the Px and Py
stdP = [100, 100];                    % standard deviation of the random variable Px and Py
lb=[0,0]; ub=[10,10];                % low bound and upper bound
xp=x0;                                    % xp is used to store the design point of previous iteration - for convergence check
Iters=0;                                  % iteration index
x_all_f =[]; x_all_g1 =[]; x_all_g2 =[];   % matrix to save the function evaluations for f, G1, G2

options = optimoptions('fmincon', 'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient',true);
[xopt,~]=fmincon(@Costfun,x0,[],[],[],[],lb,ub,@frelcon,options)
x_all_f 
x_all_g1
x_all_g2
%====================      Obj. Function   ===============================%
    function [f,g]= Costfun(x)
        f=x(1)*x(2);
        g=[x(2) x(1)];
        Cost=f;
        x_all_f = [x_all_f; [x, f]]; 
    end
%==========  Define Outer Loop Constraints and Gradiants  =================%
    function [c,ceq,GC,GCeq] = frelcon(x)
        ceq=[]; GCeq=[];
        for j = 1:nc
            [G,DG]=HL_RF(x, j);
            c(j)=bt-G;
            GC(:,j)=DG;
        end
        dx=norm(x-xp);
        if  dx>1d-5  || Iters == 0
            Iters=Iters+1;
            SHOW(Iters,x,c,GC);

        end
        xp = x;
    end
%======== Inner Loop using FORM with HL_RF Algorithm ==============%
    function [beta,dbeta]=HL_RF(x,kc)
        u=zeros(1,nd+2); iter=0;  Dif=1; sign = 1;  % here we have 2 additional random variables, Px, Py. Thus, u will be a 4 dimenional vector
        while Dif >= 1d-5 && iter < 20
            iter=iter + 1;
            [g,DG]=cons(u,x,kc);                             % call the constraint function to find out the g and DG
            u=(DG*u'-g)/norm(DG)^2*DG;
            U(iter,:)=u/norm(u);
            if iter ==1
                sign = -g/abs(g);
            elseif iter>1
                Dif=abs(U(iter-1,:)*U(iter,:)' - 1);
            end
        end
        beta = sign*norm(u);
        dbeta_All = u./(beta*[stdx, stdP]);
        dbeta =dbeta_All(1:2); 
    end
%============================ Constraint Fun. ============================%
    function [g,DG]=cons(u,d,kc)
    xtemp = u.*[stdx, stdP] +[d, uP];      
    w = xtemp(1);  t = xtemp(2); Px = xtemp(3); Py = xtemp(4);
    L = 100; E = 2.9e7; S0 = 35000; D0 = 2.5;

    if kc == 1
        g= (6 * L) * (Px / w + Py / t) / (t * w) - S0;
        DG(1)=-(6 * L * Px) / (t * w^3) - ( (6 * L) * (Py / t + Px / w))/ (t * w^2);
        DG(2)=-(6 * L * Py) / (t^3 * w) - ( (6 * L) * (Py / t + Px / w))/ (t^2 * w);   
        DG(3) = (6 * L) / (t * w^2);
        DG(4) = (6 *L) / (t^2 * w); 
        x_all_g1 = [x_all_g1; [xtemp, g, DG]]; 
    elseif kc == 2
        sqroot = sqrt((Px/w^2)^2+(Py / t^2)^2);
        g=(4 *L^3) *  sqroot / (E * t * w)- D0;
        DG(1)=-(8 * L^3 * Px^2)/(E * t * sqroot * w^6) - (4 *L^3 * sqroot)/(E *t * w^2);
        DG(2)=-(8 * L^3 * Py^2)/(E * t^6 * sqroot * w) - (4 *L^3 * sqroot)/(E *t^2 * w);
        DG(3) = 4 * L^3 * Px / ( E * t * w^5 * sqroot);
        DG(4) = 4 * L^3 * Py / (E * t^5 * w * sqroot);
        x_all_g2 = [x_all_g2; [xtemp, g, DG]];
    end
    
    end
%====== Display the Iteration Information=====%
    function  SHOW(Iters,x,c,GC)
        fprintf(1,'\n********** Iter.%d ***********\n' ,Iters);
        disp(['Des.: ' sprintf('%6.4f  ',x)]);
        disp(['Obj.: ' sprintf('%6.4f',Cost)]);
        disp(['Index.: ' sprintf('%6.4f ',bt-c)]);
        disp(['Grad.: ' sprintf('%6.4f ',GC)]);
        fprintf('\n\n')
    end
end

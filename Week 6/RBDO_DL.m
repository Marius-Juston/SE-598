% RBDO CODE with Double Loop formulation using FORM_HLRF method
function RBDO_DL()
close all; 
clear; 
clc 
global nc nd bt stdx Iters Cost x_all   % define global variables that can be used in different subroutines 
nc=3;                                     % # of constraints   
nd=2;                                     % # of design variables 
bt=norminv(0.99,0,1);              % target reliability index
x0=[5,5];                                % intial design point
stdx=[0.6,0.6];                        % standard deviation of the random variable
lb=[0,0]; ub=[10,10];              % low bound and upper bound
xp=x0;                                    % xp is used to store the design point of previous iteration - for convergence check
Iters=0;                                   % iteration index
x_all = [];
%%%% Plot the constraint contours %%%%%%%%
[xp1, xp2] = meshgrid([0:0.1:10; 0:0.1:10]);

g1=1-xp1.^2.*xp2/20;
g2=1-(xp1+xp2-5).^2/30-(xp1-xp2-12).^2/120;
g3=1-80./(xp1.^2+8*xp2+5);

figure(1) 
hold on
contour(xp1, xp2, g1, [0, 0],'-b');
contour(xp1, xp2, g2, [0, 0],'-b');
contour(xp1, xp2, g3, [0, 0],'-b');
%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimoptions('fmincon','Display', 'iter-detailed', 'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient',true);
[xopt,~]=fmincon(@Costfun,x0,[],[],[],[],lb,ub,@frelcon,options)
plot(x_all(:,1), x_all(:,2), '*b');
size(x_all)
%====================      Obj. Function   ===============================%
    function [f,g]= Costfun(x)
        f=x(1)+x(2);
        g=[1 1];
        Cost=f;
    end
%==========  Define Outer Loop Constraints and Gradiants  =================%
    function [c,ceq,GC,GCeq] = frelcon(x)
        ceq=[]; GCeq=[];
        for j = 1:nc
            [G,DG]=HL_RF(x,j);
            beta(j)=bt-G;
            dbeta(:,j)=DG;
        end
        c=beta; GC=dbeta;
        dx=norm(x-xp);
        if  dx>1d-5  || Iters == 0
            Iters=Iters+1;
%             SHOW(Iters,x,c,GC);
            plot(x(1),x(2),'or');
            pause(0.1);
        end
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
        dbeta = u./(beta*stdx);
    end
%============================ Constraint Fun. ============================%
    function [g,DG]=cons(u,d,kc)
        x = u.*stdx+d;
        x_all = [x_all; x];
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
%====== Display the Iteration Information=====%
    function  SHOW(Iters,x,c,GC)
        fprintf(1,'\n********** Iter.%d ***********\n' ,Iters);
        disp(['Des.: ' sprintf('%6.4f  ',x)]);
        disp(['Obj.: ' sprintf('%6.4f',Cost)]);
        disp(['Index.: ' sprintf('%6.4f ',bt-c)]);
        fprintf('\n\n')
    end
end

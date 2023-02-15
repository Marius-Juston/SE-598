%%%%%%%%%%  RBDO with FORM using PMA Approach  %%%%%%%%
function RBDO_PMA()
clear; close all; clc;
global nc nd bt stdx Iters Cost G1Constraints G2Constraints G3Constraints 

G1Constraints = [];
G2Constraints = [];
G3Constraints = [];


G1 = @(x1, x2 ) 1-x1.^2.*x2/20;
G2 = @(x1, x2) 1-(x1+x2-5).^2/30-(x1-x2-12).^2/120;
G3 = @(x1, x2) 1-80./(x1.^2+8*x2+5);

x1p = 0:0.1:10; x2p = 0:0.1:10;

[x1, x2] = meshgrid(x1p, x2p);

G1p = G1(x1, x2);
G2p = G2(x1, x2);
G3p = G3(x1, x2);

contour(x1, x2, G1p, [0,0])
hold on
contour(x1, x2, G2p, [0,0])
hold on
contour(x1, x2, G3p, [0,0])

nc=3;  nd=2;  bt=norminv(0.99,0,1); 
x0=[5,5]; stdx=[0.6,0.6]; lb=[0,0]; ub=[10,10]; 
xp=x0; Iters=0;
options = optimset('GradConstr','on','GradObj','on','LargeScale','off');
[x_opt,~] = fmincon(@Costfun,x0,[],[],[],[],lb,ub,@frelcon,options)

G1Constraints
G2Constraints
G3Constraints

size(G1Constraints)
size(G2Constraints)
size(G3Constraints)

%====================      Obj. Function   ===============================%
function [f,g]= Costfun(x)
    f=x(1)+x(2);
    g=[1 1];
    Cost=f;

    plot(x(1), x(2), 'or')
    pause(.5);
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
    u=zeros(1,nd); iter = 0; Dif=1;
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
end
%============================ Constraint Fun. ============================%
function [ceq,GCeq]=cons(u,d,kc)
    x = u.*stdx+d;
    if kc == 1
       ceq=1-x(1)^2*x(2)/20;
       GCeq(1)=-x(1)*x(2)/10*stdx(1);
       GCeq(2)=-x(1)^2/20*stdx(2);

       G1Constraints = [G1Constraints; x, ceq, GCeq];
    elseif kc == 2
       ceq=1-(x(1)+x(2)-5)^2/30-(x(1)-x(2)-12)^2/120;
       GCeq(1)=(-(x(1)+x(2)-5)/15-(x(1)-x(2)-12)/60)*stdx(1);
       GCeq(2)=(-(x(1)+x(2)-5)/15+(x(1)-x(2)-12)/60)*stdx(2);

       G2Constraints = [G2Constraints; x, ceq, GCeq];
      
    elseif kc == 3
       ceq=1-80/(x(1)^2+8*x(2)+5);
       GCeq(1)=x(1)*160*stdx(1)/((x(1)^2+8*x(2)+5))^2;
       GCeq(2)=80*8*stdx(2)/((x(1)^2+8*x(2)+5))^2;

       G3Constraints = [G3Constraints; x, ceq, GCeq];
    end
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

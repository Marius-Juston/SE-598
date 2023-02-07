
G1 = @(x) (x(1).^2) .*(x(2)) / 20 - 1;
G2  = @(x) (x(1) + x(2) - 5).^2 / 30 + (x(1) - x(2) - 12).^2 / 120 - 1;
G3 = @(x) 80 / (x(1).^2 + 8 * x(2) + 5) - 1;

dmin = [0,0];
dmax = [10, 10];

d0 = [5 , 5];

Outputs = [];
n = 1;

options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp', 'SpecifyObjectiveGradient', True);
[x, fval, ~, outputs] = fmincon(@Lecture3_x01, d0, [], [], [], [], dmin, dmax, @conditions ,options);

function [obj] = Lecture3_x01(x)
    global Outputs n

    obj =  x(1) + x(2); 

    Outputs(n, :) = [n, x, obj];
    n = n + 1;
end

function [c, ceq, GC, GCeq] = conditions(x)

    G1 = @(x) (x(1).^2) .*(x(2)) / 20 - 1;
    G2  = @(x) (x(1) + x(2) - 5).^2 / 30 + (x(1) - x(2) - 12).^2 / 120 - 1;
    G3 = @(x) 80 / (x(1).^2 + 8 * x(2) + 5) - 1;

    G = {G1, G2, G3};



    FORM_Normal()
    
   
   c = [g1, g2, g3];
   ceq = [];
   GCeq = [];
end

function [beta] = FORM_Normal(mu, stdx,G, DG_F)
    % define variables
    [~,nd] = size(mu); u=zeros(1,nd); iter=0;  Dif=1;
    U = [];
    
    % start the HL_RF loop
    while Dif >= 1d-5 && iter < 20
        iter=iter + 1;
        x = mu+u.*stdx;
        Gx= G(x);

        DG = zeros(1, nd);

        for i = 1:nd
            DG(i) = DG_F{i}(x);
        end

        DG = DG.*stdx;
        u=(DG*u'-Gx)/norm(DG)^2*DG;
        U(iter,:)=u/norm(u);
        if iter>1
            Dif=abs(U(iter-1,:)*U(iter,:)' - 1);
        end
    end
    beta = norm(u);  % reliability index
end
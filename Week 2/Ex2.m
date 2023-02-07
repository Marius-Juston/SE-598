%% Exercise P1
clear all; format long; clc

%% Monte Carlo Simimulation
% STEP 1: Generate N number of uniform random sample 
N = 1000000;

% STEP 2: Converting the sample to Normal Distribution
X1 = normrnd(3.5, 0.5, N, 1);
X2 = normrnd(3.8, 0.5, N, 1);

% STEP 3: Calculating the performance function G (X)
G = 1 - ((X1 + X2 - 5).^2)/30 - ((X1 - X2 - 12).^2)/120;

L = G(G < 0);
U = G(G>= 0);
histogram(L, 'FaceColor','red', 'EdgeColor','red');
hold on;
histogram(U, 'FaceColor','blue', 'EdgeColor','blue');

% STEP 4: Estimate Reliability
PoF = sum(G>0)/N;
%disp(PoF);
fprintf("P1 MCS: PoF: %f\n", PoF);

%% FORM Computation
mu = [3.5, 3.8];
stdx = [0.5, 0.5];

G = @(x)1 - ((x(1) + x(2) - 5).^2)/30 - ((x(1) - x(2) - 12).^2)/120;
DG1 = @(x) 1/60 * (32 - 5 * x(1) - 3 * x(2));
DG2 = @(x) 1/60 * (8 - 3 * x(1) - 5 * x(2));

DG_Fs = {DG1, DG2};

[beta, PoF, U] = FORM_Normal(mu, stdx, G, DG_Fs);

%disp(beta);
%disp(PoF);
% disp(U);

fprintf("P2 FORM: Beta value: %f, PoF: %f\n", beta, PoF);

%% Exercise P2 
%% Monte Carlo Simulation
% STEP 1: Generate N number of uniform random sample 
N = 1000000;
D0 = 2.5; % inches
E = 2.9e7; % psi
L = 100; % inches
w = 2; % inches
t = 4; % inches

% STEP 2: Converting the sample to Normal Distribution
X1 = normrnd(500, 100, N, 1);
X2 = normrnd(1000, 100, N, 1);

% STEP 3: Calculating the performance function G (X)
G = 4 * L^3/(E * w * t) * sqrt((X1 / w^2).^2 + (X2 / t^2).^2) - D0;

L = G(G < 0);
U = G(G>= 0);
histogram(L, 'FaceColor','red', 'EdgeColor','red')
hold on;
histogram(U, 'FaceColor','blue', 'EdgeColor','blue')

% STEP 4: Estimate Reliability
PoF = sum(G>0)/N;
%disp(PoF);

fprintf("P2 MCS: PoF: %f\n", PoF);

%% FORM Computation
D0 = 2.5; % inches
E = 2.9e7; % psi
L = 100; % inches
w = 2; % inches
t = 4; % inches
mu = [500, 1000];
stdx = [100, 100];

G = @(x)4 * L^3/(E * w * t) * sqrt((x(1) / w^2).^2 + (x(2) / t^2).^2) - D0;
DG1 = @(x) 4 * L^3 * x(1)/ (E * t * w^5 * sqrt(x(1).^2/w^4 + x(2).^2/t^4));
DG2 = @(x) 4 * L^3 * x(2)/ (E * t^5 * w * sqrt(x(1).^2/w^4 + x(2).^2/t^4));

DG_Fs = {DG1, DG2};

[beta, PoF, U] = FORM_Normal(mu, stdx, G, DG_Fs);

fprintf("P2 FORM: Beta value: %f, PoF: %f\n", beta, PoF);

%% Functions
function [beta, PoF, U] = FORM_Normal(mu, stdx,G, DG_F)
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
    PoF = 1- normcdf(beta,0,1);  % PoF value
end
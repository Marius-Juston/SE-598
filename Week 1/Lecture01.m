%% Lecture 
%% Problem 1
clc, close all, clear all

% Minimize f(x) = exp^x - 5x -1
% Solve unconstrained optimization problem to see the usage of data in 
% optiomization

% define functions
fx = @(x) exp(x) - 5 * x - 1; % optimization function
dfx = @(x) exp(x) - 5;        % first derivative
d2fx = @(x) exp(x);           % second derivative

x = .5;
tol = 1e-5;
Nmax = 20;
n = 1;

outputs = [];

while (n < Nmax)
    fx_val = fx(x);
    dfx_val = dfx(x) + normrnd(0, 1e-4, 1, 1);
    d2fx_val = d2fx(x);

    outputs(n, :) = [n, x, fx_val, dfx_val, d2fx_val];

    x = x - dfx_val / d2fx_val;

    if abs(dfx_val) < tol
        break;
    end

    n = n + 1;
end

disp(outputs);

%% Problem 2
clc, close all, clear all

% Minimize f(x) = exp^x - 5x -1
% Solve unconstrained optimization problem to see the usage of data in 
% optiomization

% define functions
fx = @(x, O) x(1)^2 + 3 * x(2)^2 + 6 * x(1) + 18 * x(2) + 22 * sin(0.1 * x(1) * x(2) + 1.5) - 20; % optimization function

x0 = [0 , 0];

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
[x, fval, exitflag, outputs] = fminunc(fx, x0,options);

x
fval
disp(outputs);

%% HW Exercises
%% Ex 1

% Clear previous data
clc, close all, clear all

%
x0 = [0 , 0];

global Outputs n;
Outputs = [];
n = 1;

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
[x, fval, exitflag, outputs] = fminunc(@Lecture1_x01, x0,options);

x
fval
disp(outputs);
Outputs


%% Exercise 2

clc, close all, clear all

% f(x) = x1 + x2
% g1 =  1- x1^2 * x2 / 20
% g2 = 1 - (x1 + x2 - 5) ^2/30 - (x1 - x2 - 12)^2/120
% g3 = 1-80 / (x1^2 + 8*x2 +5)
% x = [x1 x2]
% x1 in [0, 10]
% x2 in [0, 10]

global Outputs n;
Outputs = [];
n = 1;

x0 = [5  5];

options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp');
[x, fval, exitflag, outputs] = fmincon(@Lecture1_x02, x0, [], [], [], [], [0 ,0], [10, 10], @conditions ,options);

x
fval
disp(outputs);
Outputs

%% Functions

%% Ex 1 Function
% Objective function for ex1 as well as recorder method
function [obj] = Lecture1_x01(x)
    % Uses global variables for recording
    global Outputs n

    % Calculates the objective function
    obj = x(1)^2 + 3 * x(2)^2 + 6 * x(1) + 18 * x(2) + 22 * sin(0.1 * x(1) * x(2) + 1.5) - 20; 

    % Appends to the output matrix the current input information
    Outputs(n, :) = [x, obj];
    n = n + 1;
end

%% Ex 2 Functions

function [obj] = Lecture1_x02(x)
    global Outputs n

    obj =  x(1) + x(2); 

    Outputs(n, :) = [x, obj];
    n = n + 1;
end

function [c, ceq] = conditions(x)
   g1 =  1- x(1)^2 * x(2) / 20;
   g2 = 1 - (x(1) + x(2) - 5) ^2/30 - (x(1) - x(2) - 12)^2/120;
   g3 = 1-80 / (x(1)^2 + 8*x(2) +5);
   
   c = [g1, g2, g3];
   ceq = [];
end
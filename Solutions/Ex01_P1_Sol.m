function Ex01_P1_Sol()
close all; clear;  clc
global objVal xVal     % Define two global variables to record obj and X iteratively 
objVal = [];             % initialize the global variable
xVal=[]; 
x0 = [0,0];              % set the starting point x0
options = optimoptions(@fminunc,'Algorithm','quasi-newton'); 
[x,fval,exitflag,output] = fminunc(@objfun,x0,options)     % using the Quasi-Newton method to find the minimum

figure(1)       % plot the objective function values iteratively 
plot(1:length(objVal), objVal, '-*')
title('History of Ojective Function Value'); 
xlabel('Iterations')
ylabel('Objective')

figure(2)       % plot the sample points iteratively 
plot(xVal(:,1), xVal(:,2), '-o' )
title('History of Sample Points'); 
xlabel('X_1')
ylabel('X_2')

function obj = objfun(x)   % here we define the objective function with input x and output obj
global xVal objVal
obj =  x(1).^2+3*x(2).^2+6*x(1)+18*x(2)+22*sin(0.1*x(1).*x(2)+1.5)-20;
xVal =[xVal; x];
objVal = [objVal; obj]; 

% Note: you can also define the obj function using "@", without using "function", such as 
% obj = @(x) x(1).^2+3*x(2).^2+6*x(1)+18*x(2)+22*sin(0.1*x(1).*x(2)+1.5)-20;



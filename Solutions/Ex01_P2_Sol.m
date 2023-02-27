function Ex01_P2_Sol()
close all; clear; clc
x0=[5, 5];                       % define starting point x0
lb = [0, 0];                      % define low bound of X
ub = [10, 10];                  % define upper bound of X
objfun = @(x) x(1)+x(2);  % define objective fun

% using the SQP  to solve nonlinear constrained optimization problem. 
% see fmincon in help document for more details 
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
x = fmincon(objfun,x0,[],[],[],[],lb,ub,@nonlcon,options)


function [c, ceq]=nonlcon(x)   % here we define the nonlinear constraints  
ceq = [];                               % no equality constraint  
c(1)=1-x(1)^2*x(2)/20;
c(2)=1-(x(1)+x(2)-5)^2/30-(x(1)-x(2)-12)^2/120;
c(3)=1-80/(x(1)^2+8*x(2)+5);

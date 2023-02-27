% Initial trial with surrogate modeling 
close all; clear; clc

% we first plot the true constraint contours for G = 0; 
g=@(x1, x2) 1-(x1+x2-5).^2/30-(x1-x2-12).^2/120;
[xp1, xp2] = meshgrid([0:0.1:10; 0:0.1:10]);
figure(1) 
hold on;
contour(xp1, xp2, g(xp1,xp2),'-b', 'ShowText','on');

% Step 1: 
% model form: 3rd order polynomials  
% y = w1+w2*x1 + w3*x2 + w4*x1^2+x5*x2^2+w6*x1*x2
% prepare the data
xsamp = 1:1:10; 
[xs1, xs2] = meshgrid(xsamp, xsamp);  % generating sample points
xs=[reshape(xs1, [], 1), reshape(xs2, [],1)]; 
ys = g(xs(:,1), xs(:,2));

% fit the model parameters using least square
fun0 = @(w, x1, x2) w(1)+w(2)*x1+w(3)*x2 +w(4)*x1.^2+ ... 
             w(5)*x2.^2 +w(6)*x1.*x2;
fun = @(w) sum((fun0(w, xs(:,1), xs(:,2))-ys).^2);

x0 = ones(1,6);
[w1,~] = fminunc(fun,x0)   %min MSE

% plot contours using predicted models 

hold on;
contour(xp1, xp2, fun0(w1,xp1,xp2),'--r');
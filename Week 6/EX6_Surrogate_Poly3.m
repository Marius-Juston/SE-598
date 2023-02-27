% Initial trial with surrogate modeling 
close all; clear; clc

% we first plot the true constraint contours for G = 0; 
L = 100; E = 2.9e7; S0 = 35000; D0 = 2.5;
Px = 500; Py = 1000;

g1 = @(w,t)(6 * L) * (Px ./ w + Py ./ t) ./ (t .* w) - S0;
g2 =@(w,t)(4 *L^3) *  sqrt((Px./w.^2).^2+(Py ./ t.^2).^2) ./ (E * t .* w)- D0;

[w, t] = meshgrid([1:0.1:10; 1:0.1:10]);
figure(1) 
hold on;
y1 = g1(w,t);
y2 = g2(w,t);
contour(w, t, y1,[0,0], '-b', 'ShowText','on' );
hold on
contour(w, t, y2,[0,0],'-r', 'ShowText','on');


%% Polynomial 3rd order weight calcualtion using w = X^+ y solution 
xsamp = linspace(1, 10, ceil(sqrt(300))); 
[xs1, xs2] = meshgrid(xsamp, xsamp);  % generating sample points
xs=[reshape(xs1, [], 1), reshape(xs2, [],1)]; 
ys1 = g1(xs(:,1), xs(:,2));
ys2 = g2(xs(:,1), xs(:,2));

xE = [ones(size(xs, 1),1), xs, xs.^2, xs(:,1).*xs(:,2)];

w1 = pinv(xE) * ys1;
w2 = pinv(xE) * ys2;

% plot contours using predicted models 

ys1_pred = xE * w1;
ys1_pred_z = reshape(ys1_pred, size(xs1));

ys2_pred = xE * w2;
ys2_pred_z = reshape(ys2_pred, size(xs2));

hold on;
contour(xs1, xs2, ys1_pred_z,[0,0], '--b');
contour(xs1, xs2, ys2_pred_z,[0,0], '--r');


%% RMSE Calculation

[w, t] = meshgrid([1:2:10; 1:2:10]);
y1 = g1(w,t);
y2 = g2(w,t);
xs=[reshape(w, [], 1), reshape(t, [],1)]; 
xE = [ones(size(xs, 1),1), xs, xs.^2, xs(:,1).*xs(:,2)];

ys1_pred = xE * w1;
y1 = reshape(y1, size(ys1_pred));

ys2_pred = xE * w2;
y2 = reshape(y2, size(ys2_pred));

rmse1 = sqrt(mean((y1 - ys1_pred).^2));
rmse2 = sqrt(mean((y2 - ys2_pred).^2));

disp(rmse1);
disp(rmse2);


% Comment on the RMSE:
% The accuracy of both surrogate models using the 3rd order polynomial is
% not good. For the constraint G1 the error is 1.255e+5 which is very
% large, informing us that the system is not being properly modeled. The
% constraint G2 is being better represented using the 3d order polynomial
% with a much smaller error of 22.8; however, this is still very large
% error. As such the exploration of either a higher model or a different
% model function should be explored to achieve more accurate results.


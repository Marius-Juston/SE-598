clear; close all; clc
Pi = [0.87, 0.75, 0.92, 0.90, 0.95, 0.93, 0.88, 0.79, 0.99];

invPi = 1 - Pi;

lower_bd = max(invPi);
upper_bound = min(sum(invPi), 1);

disp(lower_bd);
disp(upper_bound);
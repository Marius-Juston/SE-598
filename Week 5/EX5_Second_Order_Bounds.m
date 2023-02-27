clear; close all; clc

Pi = [0.87, 0.75, 0.92, 0.90, 0.95, 0.93, 0.88, 0.79, 0.99];
% Pi = [0.82, 0.78, 0.9, 0.72, 0.95];
n = size(Pi, 2);

invPi = 1 - Pi;

invPi = sort(invPi, "descend");


upper_bound = min(sum(invPi), 1);

lower_bd = invPi(1);

for i=2:n
    sub = invPi(i);

    temp = 0;

    for j=1:i-1
        temp = temp + (sub * invPi(j));
    end

    lower_bd = lower_bd + max(sub - temp, 0);
end

upper_bd = sum(invPi);

max_sec = invPi(1).*invPi;
max_sec = max_sec(2:n);
upper_bd = min(upper_bd -  sum(max_sec), 1);

disp(lower_bd);
disp(upper_bd);
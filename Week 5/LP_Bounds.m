clear; close all; clc
Pi = [0.72, 0.78, 0.82, 0.90, 0.95];
Pij = zeros(5,5);
 for i = 2:5
     for j = 1:i-1
         Pij(i,j) = Pi(i)*Pi(j);
     end
 end 
% binary decomposition with 5 components
bd_mat = [1 1 1 1 1; 0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0;
    0 0 1 1 1; 0 1 0 1 1; 0 1 1 0 1; 0 1 1 1 0; 1 0 0 1 1; 1 0 1 0 1;
    1 0 1 1 0; 1 1 0 0 1; 1 1 0 1 0; 1 1 1 0 0; 0 0 0 1 1; 0 0 1 0 1;
    0 0 1 1 0; 0 1 0 0 1; 0 1 0 1 0; 0 1 1 0 0; 1 0 0 0 1; 1 0 0 1 0;
    1 0 1 0 0; 1 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0;
    1 0 0 0 0; 0 0 0 0 0];
f = [1, zeros(1,31)];
A = -eye(32); b = zeros(32,1);
Aeq = []; beq = [];
for ki = 1:5
    for kj = 1:ki
        Aeq=[Aeq; (bd_mat(:, ki)>0 & bd_mat(:, kj)>0)'];
        if ki == kj
            beq = [beq, Pi(ki)];
        else
            beq = [beq, Pij(ki, kj)];
        end
    end
end
x_low = linprog(f,A,b,Aeq,beq);
x_up = linprog(-f,A,b,Aeq,beq);
LPB = [1-x_up(1), 1-x_low(1)]
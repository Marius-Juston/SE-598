clear; close all; clc
n = 9;
Pi = [0.87, 0.75, 0.92, 0.90, 0.95, 0.93, 0.88, 0.79, 0.99];
Pij = zeros(n,n);
 for i = 2:n
     for j = 1:i-1
         Pij(i,j) = Pi(i)*Pi(j);
     end
 end 
% binary decomposition with n components
bd_mat = [1 1 1 1 1; 0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0;
    0 0 1 1 1; 0 1 0 1 1; 0 1 1 0 1; 0 1 1 1 0; 1 0 0 1 1; 1 0 1 0 1;
    1 0 1 1 0; 1 1 0 0 1; 1 1 0 1 0; 1 1 1 0 0; 0 0 0 1 1; 0 0 1 0 1;
    0 0 1 1 0; 0 1 0 0 1; 0 1 0 1 0; 0 1 1 0 0; 1 0 0 0 1; 1 0 0 1 0;
    1 0 1 0 0; 1 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0;
    1 0 0 0 0; 0 0 0 0 0];

bd_mat = [ones(1, n)];

n_2 = 2^n;

for i=1:n
    temp = ones(1, n);
    temp(i) = 0;

    bd_mat = [bd_mat; temp];
end

for i=1:n
    for j=i + 1:n
        temp = ones(1, n);
        temp(i) = 0;
        temp(j) = 0;
    
        bd_mat = [bd_mat; temp];
    end
end

for i=1:n
    for j=i + 1:n
        for e=j + 1:n
            temp = ones(1, n);
            temp(i) = 0;
            temp(j) = 0;
            temp(e) = 0;
        
            bd_mat = [bd_mat; temp];
        end
    end
end

for i=1:n
    for j=i + 1:n
        for e=j + 1:n
            for c=e + 1:n
                temp = ones(1, n);
                temp(i) = 0;
                temp(j) = 0;
                temp(e) = 0;
                temp(c) = 0;
            
                bd_mat = [bd_mat; temp];
            end
        end
    end
end

for i=1:n
    for j=i + 1:n
        for e=j + 1:n
            for c=e + 1:n

                for r=c + 1:n
                temp = ones(1, n);
                temp(i) = 0;
                temp(j) = 0;
                temp(e) = 0;
                temp(c) = 0;
                temp(r) = 0;
            
                bd_mat = [bd_mat; temp];
                end
            end
        end
    end
end

for i=1:n
    for j=i + 1:n
        for e=j + 1:n
            for c=e + 1:n
                for r=c + 1:n
                     for t=r + 1:n
                        temp = ones(1, n);
                        temp(i) = 0;
                        temp(j) = 0;
                        temp(e) = 0;
                        temp(c) = 0;
                        temp(r) = 0;
                        temp(t) = 0;
                    
                        bd_mat = [bd_mat; temp];
                     end
                end
            end
        end
    end
end

for i=1:n
    for j=i + 1:n
        for e=j + 1:n
            for c=e + 1:n
                for r=c + 1:n
                     for t=r + 1:n
                         for y=t + 1:n
                            temp = ones(1, n);
                            temp(i) = 0;
                            temp(j) = 0;
                            temp(e) = 0;
                            temp(c) = 0;
                            temp(r) = 0;
                            temp(t) = 0;
                            temp(y) = 0;
                        
                            bd_mat = [bd_mat; temp];
                         end                 
                     end
                end
            end
        end
    end
end

for i=1:n
    for j=i + 1:n
        for e=j + 1:n
            for c=e + 1:n
                for r=c + 1:n
                     for t=r + 1:n
                         for y=t + 1:n
                             for l=y + 1:n
                                temp = ones(1, n);
                                temp(i) = 0;
                                temp(j) = 0;
                                temp(e) = 0;
                                temp(c) = 0;
                                temp(r) = 0;
                                temp(t) = 0;
                                temp(y) = 0;
                                temp(l) = 0;
                            
                                bd_mat = [bd_mat; temp];
                             end    
                         end
                     end
                end
            end
        end
    end
end


bd_mat = [bd_mat; zeros(1, n)];

f = [1, zeros(1,n_2 - 1)];
A = -eye(n_2); b = zeros(n_2,1);
Aeq = []; beq = [];
for ki = 1:n
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
clc;
clear;
close all;

%Number of Initial Population & Dimensions
pop_size = 100;
dim = 2;

%Iteration Condition
max_iter = 200;

%Mutation Rate
%m_rate = 0.05;

%Domain of Benchmarks
from = -5.12;
to = -1*from;

%Results for n Times Execution
num_of_result = 5;
%Columns of total Result : dim,(gbest_fitness),(time)
total_result = zeros(num_of_result,dim+2);

for n=1:num_of_result
    tic;
    
    %Nfe Condition
    max_nfe = 20000;
    
    %Initialize Best Fitness and Position by a Large Value
    g_best = zeros(1, dim);
    g_best(1,:) = to;
    g_best_fitness = F1(g_best(1,:)); %nfe++
    nfe = 1;
    
    %Initialize Population
    X = unifrnd(from,to,[pop_size dim]);
    %Initialize Personal Bests (Equal to First Positions)
    p_best = X;
    %Initialize Velocity (Equal to First Position)
    V = X;
    %V = zeros(pop_size,dim);
    
    F_result = zeros(1, pop_size);
    
    %c_max = 2.5;
    %c_min = 0.5;
    
    %Calculate Fitness of X and Save the Best Fitness and its Position
    for i = 1:pop_size
        F_result(1,i) = F1(X(i,:));
        nfe = nfe + 1;
        if (F_result(1,i) <= g_best_fitness)
            g_best_fitness = F_result(1,i);
            g_best(1,:) = X(i,:);
        end
    end
    
    %Main Loop
    for m=1:max_iter
        F1_result = zeros(1, pop_size);
        r1 = rand;
        r2 = rand;
        
        %Cognitive Coefficient
        %c1 = (c_min - c_max) * (j/max_iter) + c_max;
        %Global Coefficient
        %c2 = (c_max - c_min) * (j/max_iter) + c_max;
        %Inertia Coefficient for Standard PSO
        %w = (c1*r1) + (c2*r2);
        
        r_w = rand;
        w0 = 0 + (0.4-0)*rand;
        alpha0 = 0.5 + (1-0.5)*rand;
        
        %Interia Coefficient for AWPSO
        w = w0 + r_w * (1 - w0);
        %Acceleration Factor
        alpha = alpha0 + m/max_iter;
        
        %Convex (for Breeding PSO)
        %lambda1 = rand;
        %lambda2 = 1 - lambda1;
        %Average
        %lambda1 = 0.5;
        %lambda2 = lambda1;
        %Affine
        %lambda1 = 1.5;
        %lambda2 = -0.5;
        %Linear
        %lambda1 = rand;
        %lambda2 = rand;
        
        for j=1:pop_size
            %Equation of Velocity (Update Velocity of each Particle)
            %V(j,:) = (w * V(j,:)) + (c1*r1*(p_best(j,:) - X(j,:)))...
            %    + (c2*r2*(g_best(1,:) - X(j,:) )); %Standard PSO
            V(j,:) = (w * V(j,:)) + alpha*((r1*(p_best(j,:) - X(j,:)))...
                + (r2*(g_best(1,:) - X(j,:) ))); %AWPSO
            
            %Update Velocity of each Particle for Breeding PSO
            %for z=1:2:pop_size
            %    V(z,:) = ((V(z,:) + V(z+1,:))/(norm(V(z,:) + V(z+1,:))))...
            %        * norm(V(z,:));
            %end
            
            %Control the Domain of the new Velocity
            for p=i:dim
                r_v = rand;
                if (V(j,p) < from)
                    V(j,p) = from + r_v;
                elseif (V(j,p) > to)
                    V(j,p) = to - r_v;
                end
            end
            
            %Equation of new Position (Update Position of each Particle)
            X(j,:) = X(j,:) + V(j,:);
            
            %Update Position of each Particle for Breeding PSO
            %for z=1:2:pop_size
            %    x_temp_a = X(z,:);
            %    x_temp_b = X(z+1,:);
            %    X(z,:) = lambda1*(x_temp_a) + lambda2*(x_temp_b);
            %    X(z+1,:) = lambda1*(x_temp_b) + lambda2*(x_temp_a);
            %end
            
            %Mutation
            %y = unifrnd(from,to,[floor(m_rate*pop_size) dim]);
            %for i = 1:size(y,1)
            %rand_Mutation = randi(pop_size);
            %X(rand_Mutation,1:dim) = X(rand_Mutation,1:dim) + y(i,:);
            %end
            
            %Control the Domain of the new Positions
            for t=1:dim
                r_x = rand;
                if (X(j,t) < from)
                    X(j,t) = from + r_x;
                elseif (X(j,t) > to)
                    X(j,t) = to - r_x;
                end
            end
        end
        
        %Update the Personal Bests and Global Best
        for k = 1:pop_size
            
            if (nfe >= max_nfe)
                break;
            end
            
            F1_result(1,k) = F1(X(k,:));
            nfe = nfe + 1;
            if (F1_result(1,k) <= F1(p_best(k,:)))
                nfe = nfe +1;
                p_best(k,:) = X(k,:);
            end
            
            if (F1_result(1,k) <= g_best_fitness)
                g_best_fitness = F1_result(1,k);
                nfe = nfe +1;
                g_best(1,:) = X(k,:);
            end
        end
    end
    
    total_result(n,1) = toc;
    total_result(n,2) = g_best_fitness;
    total_result(n,3:end) = g_best;
end

min_fitness = min(total_result(:,2));
max_fitness = max(total_result(:,2));
mean_fitness = mean(total_result(:,2));
std_fitness = std(total_result(:,2));
mean_time = mean(total_result(:,1));

disp(strcat('Popsize:', num2str(pop_size), ', Dimension:', num2str(dim)));
disp(strcat('mean fitness: ', num2str(mean_fitness)));
disp(strcat('max fitness: ', num2str(max_fitness)));
disp(strcat('min fitness: ', num2str(min_fitness)));
disp(strcat('std fitness: ', num2str(std_fitness)));
disp(strcat('mean time: ', num2str(mean_time)));
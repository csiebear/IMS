clear;
% Problem Definition and Parameters Setting
PT = [10,10,13,4,9,4,8,15,7,1,9,3,15,9,11,6,5,14,18,3];
DUE= [50,38,49,12,20,105,73,45,6,64,15,6,92,43,78,21,15,50,150,99];
WEIGHT = [10,5,1,5,10,1,5,10,5,1,5,10,10,5,1,10,5,5,1,5];
ma_num= 1;% Single_Machine
j_num = size(PT,2);

% Setting Genetic Algorithm parameter
% 母體大小(population_size),交配率(crossover_rate),突變率(mutation_rate)
% 重複執行的次數(Num_Iteration)
population_size = 30;
crossover_rate = 0.9;
mutation_rate = 0.1;
Num_Iteration = 50000;

tic % Start stopwatch timer

% the population_list is a population_size-by-j_num array
% It is used to put the random job sequence (In GA , we called it chromosome)
population_list = zeros(population_size, j_num); % Record population (includes all chromosomes).

% Initialize the Solution
% Pick one popultaion_list and put row vector containing a random
% permutation of integers from 1 to j_num
for m = 1:population_size
    population_list(m,1:j_num) = randperm(j_num);
    % disp(population_list(m,1:j_num));
    for j = 1:j_num 
        for k = 0:(j_num-1)
            if population_list(m,j) > k && population_list(m,j) <= (k+1)
                population_list(m,j) = k+1;
            end
        end
    end
end

Makespan_best = 9999999999; 
Scheduling_best = zeros(1, j_num);  

% Crossover->Mutation->Repairment->Evaluation->Selection loop and we will
% run (Num_iteration) times
%for i = 1:Num_Iteration
    %record the original population_list in population_list_tmp
    population_list_tmp = population_list;
    % Crossover
    S = randperm(population_size); 
    
    for m = 1:(population_size/2)
        crossover_prob = rand();
        if (crossover_rate >= crossover_prob) 
            % pick two chromosomes from popultaion_list become a pair of parent 
            parent_1 = population_list(S(-1+2*m),1:j_num*ma_num); 
            parent_2 = population_list(S(2*m),1:j_num*ma_num); 
            %assume the child_1 and child_2 are same as their parent
            child_1 = parent_1; 
            child_2 = parent_2; 
            % randomly decide the cutpoint
            cutpoint = randperm(j_num); 
            % choose 2 numbers as start and end point
            cutpoint = sort(cutpoint(1:2));
            % Show the cut range
            % disp(cutpoint);
            
            % You can see the change if you remove the comment mark %
            % replace child string with the other parent string in the
            % specified range
            for k = cutpoint(1):cutpoint(2) 
                % disp('before');
                % fprintf('The k =%d \n',k);
                % disp(child_1);
                % disp(parent_2(k));
                child_1(k) = parent_2(k);
                % disp('After');
                % disp(child_1);
                child_2(k) = parent_1(k);
            end
            % store the child string(chromosome)
            population_list(S(-1+2*m),1:j_num*ma_num) = child_1;
            population_list(S(2*m),1:j_num*ma_num) = child_2;
        end%end if
    end%end for-loop
    
    % Mutation
    for m = 1:population_size
        for j = 1:j_num*ma_num
            mutation_prob = rand();  
            if mutation_rate >= mutation_prob  
                ran_num = rand(); 
                for k = 0:(j_num-1)
                    if ran_num > k*(1/j_num) && ran_num <= (k+1)*(1/j_num)
                        population_list(m,j) = k+1;
                    end
                end
            end
        end
    end
    
    % Repairment
    for m = 1:population_size 
        repair_or_not = zeros(1, j_num);
        for j = 1:j_num*ma_num
            for k = 1:j_num
                if population_list(m,j) == k
                    repair_or_not(k) = repair_or_not(k) + 1;
                end
            end
        end
        
        for k = 1:j_num
            if repair_or_not(k) > ma_num
                r_ran_num = randperm(repair_or_not(k)); 
                r_ran_num = sort(r_ran_num(1:(repair_or_not(k)-ma_num)));
                appeartime = 0;  
                for j = 1:j_num*ma_num
                    if population_list(m,j) == k
                        appeartime = appeartime + 1;
                        for n = 1:(repair_or_not(k)-ma_num)
                            if appeartime == r_ran_num(n)
                                population_list(m,j) = 0;   
                            end
                        end
                    end
                end
            end
        end
        
        for k = 1:j_num 
            if repair_or_not(k) < ma_num
                zeroappeartime = 0;
                appeartime = 0;      
                for j = 1:j_num*ma_num  
                    if population_list(m,j) == 0
                        zeroappeartime = zeroappeartime + 1;
                    end
                end
                r_ran_num = randperm(zeroappeartime);  
                r_ran_num = sort(r_ran_num(1:ma_num-(repair_or_not(k))));
                for j = 1:j_num*ma_num
                    if population_list(m,j) == 0
                        appeartime = appeartime + 1;
                        for n = 1:(ma_num-repair_or_not(k))
                            if appeartime == r_ran_num(n)
                                population_list(m,j) = k;
                            end
                        end
                    end
                end
            end
        end
    end
%end

% Report the Results
disp('--- Final Report ---');
fprintf('Optimal_Value : %d\n',Makespan_best);

toc % Read elapsed time from stopwatch

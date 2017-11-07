clear;
% Problem Definition and Parameters Setting
PT = [10;10;13;4;9;4;8;15;7;1;9;3;15;9;11;6;5;14;18;3];
DUE= [50;38;49;12;20;105;73;45;6;64;15;6;92;43;78;21;15;50;150;99];
WEIGHT = [10;5;1;5;10;1;5;10;5;1;5;10;10;5;1;10;5;5;1;5];
ma_num= 1;% Single_Machine
Ma=ones(size(PT,1),1);
j_num = size(PT,1);

% Setting Genetic Algorithm parameter
% 母體大小(population_size),交配率(crossover_rate),突變率(mutation_rate)
% 重複執行的次數(Num_Iteration)
population_size = 30;
crossover_rate = 0.8;
mutation_rate = 0.1;
Num_Iteration = 10000;
Tbest=99999999;%Objective value default setting to max
num_tardy=99999999;
Best_position=0;

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

% Crossover->Mutation->Repairment->Evaluation->Selection loop and we will
% run (Num_iteration) times
for i = 1:Num_Iteration
    
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
    % 計算tardiness
    tardiness_list=zeros(1,population_size);
    fitness_list=zeros(1,population_size);
    pk_list=zeros(1,population_size);
    q=zeros(1,population_size);
    R=zeros(1,population_size);
    total = 0 ;
    
    for m = 1:population_size
        J=population_list(m,:);
        Ptime = 0 ;
        Tardiness = 0;
        for i = 1:j_num;
            Ptime = Ptime + PT(J(1,i),1);
            Tardiness = Tardiness + WEIGHT(J(1,i),1)*max(Ptime-DUE(J(1,i),1),0);
            % fprintf('The Ptime: %d The weight: %d The Due: %d The Tardiness:%d \n',Ptime,WEIGHT(J(1,i),1),DUE(J(1,i),1),Tardiness);
        end
        %計算完的Tardiness是要最小為最佳，因此這邊我們將數值取倒數方便找出最小值
        tardiness_list(m) = Tardiness;
        fitness_list(m)= 1/Tardiness;
    end    

    for i = 1:population_size   % total時間計算
        total = total+fitness_list(i);
    end
    
    for i = 1:population_size   % 機率
        pk_list(i) = fitness_list(i) / total ;
    end

    % The summation of the pk_list always equals to 1
    %disp(sum(pk_list));
    
    % 計算u值
    for  i = 1:population_size
        u=0 ;
        for j = 1:i
            u=u+pk_list(i);
        end
        % Step3: Calculate cumulative probability qk for each chromosome vk
        % The array q is the cummlative probability value
        q(i)= u ;
        % Step4: Generate a random number r from the range [0, 1].
        % The array R is the sequence random number generated by the uniformly
        % distributed random numbers
        R(i)=rand;
    end
    
    clone_population_list= population_list;

    % roulette wheel approach Step5:
    %If r<=q1,then select the first chromosome v1;
    %otherwise, select the kth chromosome vk (2<=k<=popSize) such that qk-1<r<= qk
    for i = 1 :population_size
        if (R(i)<= q(1))
            population_list(i,:) = clone_population_list(i ,:);
        else 
            for j = 1:(population_size-1);
                if (R(i)>q(j) && R(i)<= q(j+1))
                    population_list(i,:) = clone_population_list(j+1,:) ;
                break;
                end
            end
        end
    end 
    %Best result record
    min_tardiness=sort(tardiness_list());
    ThisIterationBest=min_tardiness(1);
    if ThisIterationBest<Tbest
        Tbest=ThisIterationBest;
        for i=1:population_size
            if Tbest==tardiness_list(i)
                Best_position=i;
            end
        end
        Best_job_sequence=population_list(Best_position,: );
        % Calculate the num of tardy
        num_tardy=0;
        for b = 1:j_num;
            Ptime = Ptime + PT(Best_job_sequence(1,b),1);
            if (Ptime > DUE(Best_job_sequence(1,b)) )
                num_tardy = num_tardy + 1;
            end
        end
    end
end%end the iteration

% Report the Results
disp('--- Final Report ---');
disp('Optimal Solution ( i.e., job sequence) = '); 
    disp(Best_job_sequence);
fprintf('Optimal function ( i.e., fitness) Value : %d\n',Tbest);
fprintf('Running time : %.10f\n',toc);
fprintf('Average (Weighted) Tardiness : %.2f\n',Tbest/j_num);
fprintf('Number of Tardy : %d\n',num_tardy);
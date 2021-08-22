% S O M A - T 3 A
% Code accompanying the paper:
% Diep Q.B., Zelinka I., Das S., Senkerik R. (2020) SOMA T3A for Solving 
% the 100-Digit Challenge. In: Zamuda A., Das S., Suganthan P., Panigrahi B. (eds) 
% Swarm, Evolutionary, and Memetic Computing and Fuzzy and Neural Computing. 
% SEMCCO 2019, FANCCO 2019. Communications in Computer and Information Science, 
% vol 1092. Springer, Cham. https://doi.org/10.1007/978-3-030-37838-7_14
function [Best , array_digit , FEs , Mig] = optimizerSOMA_T3A(Info , SOMApara , CostFunction, lBounds, uBounds, varargin)
    % -------------- Extract Information ----------------------------------
                    the_func            = Info.the_func;
                    FEs_Max             = Info.FEs_Max;
                    dimension           = Info.dimension;
                    f_star              = Info.f_star;
    % -------------- The domain of the function ---------------------------
                    VarMin              = Info.Search_Range(1);
                    VarMax              = Info.Search_Range(2);
    % -------------- Initial Parameters of SOMA ---------------------------
                    PopSize             = SOMApara.PopSize;
                    N_jump              = SOMApara.N_jump;
					m                   = SOMApara.m;
					n                   = SOMApara.n;
					k                   = SOMApara.k;
    % --------------------- Create Initial Population ---------------------
    pop               = repmat(lBounds, 1, PopSize) + rand(dimension,PopSize).*repmat(uBounds-lBounds, 1, PopSize);
    fitness = evaluate_fitness(CostFunction, pop, varargin{:});
    
    FEs_count         = PopSize;
    [global_cost, id] = min(fitness);
    global_leader     = pop(:,id);
    % ---------------- SOMA MIGRATIONS ------------------------------------
    FEs.reach_digit(10) = 0;
    flag_digit(10)      = 0;
    error               = Inf;
    Mig                 = 0;
    m = 50;
    k = 100;
    while (FEs_count < FEs_Max)
        Mig   = Mig  + 1;
		% ------------ Control parameters ---------------------------------
       %PRT   = 0.33 + 0.25*cos(2*pi*T*1e-7*FEs_count + pi);
        Step  = 0.02 + 0.005*cos(0.5*pi * 1e-7 * FEs_count);
        PRT   = 0.08 + 0.500*(FEs_count / FEs_Max);
       %Step  = 0.025 - 0.01*(FEs_count / FEs_Max);
        % ------------ Migrant selection: m -------------------------------
        A               = randi([1 PopSize],1,m);
        sub_pop_sortA   = [fitness(A) ; pop(:,A) ; A];
        sub_pop_sortA   = sortrows(sub_pop_sortA')';
        sub_popA        = sub_pop_sortA(2:end-1,:);
        pop_A           = sub_pop_sortA(end,:);
        % ------------ movement of n Migrants -----------------------------
        for j  = 1 : n
            indi_moving     = sub_popA(:,j);
            % ------------ Target selection: k ----------------------------
            B               = randi([1 PopSize],1,k);
            sub_pop_sortB   = [fitness(B) ; pop(:,B) ; B];
            sub_pop_sortB   = sortrows(sub_pop_sortB')';
            sub_popB        = sub_pop_sortB(2:end-1,:);
            pop_B           = sub_pop_sortB(end,:);
            Target          = sub_popB(:,1);
            % ------------- Moving process --------------------------------
            coincide        = pop_A(j) - pop_B(1);
            if  coincide
                indi_journey = [];
                for move     = 1 : N_jump
                    nstep    = move * Step;
                    %----- SOMA Mutation ----------------------------------
                    PRTVector      = rand(dimension,1) < PRT;
                    %----- SOMA Crossover ---------------------------------
                    indi_offspring =  indi_moving + (Target - indi_moving)*nstep.*PRTVector;
                    indi_journey   = [indi_journey   indi_offspring];
                end % END JUMPING
                %----- Checking Boundary and Replaced Outsize Individuals -
                number_indi_journey = size(indi_journey,2);
                for cl = 1 : number_indi_journey
                    for rw = 1 : dimension
                        if  (indi_journey(rw,cl) < lBounds(rw)) || (indi_journey(rw,cl) > uBounds(rw))
                             indi_journey(rw,cl) = lBounds(rw)  +  rand*(uBounds(rw) - lBounds(rw));
                        end
                    end
                end
                %----- SOMA Re-Evaluate Fitness Fuction -------------------
                new_cost               = evaluate_fitness(the_func, indi_journey, varargin{:});
                FEs_count              = FEs_count + number_indi_journey;
                %----- SOMA Accepting: Place Best Individual to Population-
                [min_new_cost, idzz]   = min(new_cost);
                if  min_new_cost       < fitness(pop_A(j))
                    indi_replace       = indi_journey(:,idzz);
                    pop(:,pop_A(j))    = indi_replace;
                    fitness(pop_A(j))  = min_new_cost;
                    %----- SOMA Update Global_Leader ----------------------
                    if  min_new_cost   < global_cost
                        global_cost    = min_new_cost;
                        global_leader  = indi_replace;
                        error          = global_cost - f_star;
                    end
                end
            end % END indi_moving ~= Target
            %----- SOMA Record FEs value ----------------------------------
            if error <= 1
                for i = 1:10
                    if error < 10^(-i+1) && flag_digit(i) == 0
                        FEs.reach_digit(i)  = FEs_count;
                        flag_digit(i)       = 1;
                    end
                end
            end
            %----- SOMA Checking Stop Conditon ----------------------------
            if flag_digit(10), break; end
        end  % END PopSize   (For Loop)
        if flag_digit(10),  break; end
        %if FEs_count > 3e8, break; end
    end   % END MIGRATIONS (While Loop)
    Best.Value      = global_cost;
    Best.Positon    = global_leader;
    FEs.stop        = FEs_count;
    array_digit     = flag_digit;
end
%--------------------------------------------------------------------------


function fitness = evaluate_fitness(CostFunction, pop, varargin)
[~, PopSize] = size(pop);
 fitness           = NaN(1, PopSize);
    for i =1:PopSize
       fitness(i) = feval(CostFunction, pop(:, i), varargin{:});
    end
end
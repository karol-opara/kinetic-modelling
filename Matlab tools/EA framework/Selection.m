classdef Selection
    methods (Static)
        function offspringPopulation = Proportional...
                (parentPopulation, ... % parent population matrix - one individual is one column
                parentFitness, ...     % parent fitness
                offspringNumber)       % (optional) number of offsprings
            
            [~, popSize] = size(parentPopulation);
            if (nargin == 2)
                offspringNumber = popSize;
            end
            
            % first, let's make the fitness nonnegative
            % (with this transformation the least-fit individual cannot reproduce)
            fx = parentFitness - min(parentFitness);
            
            sumOfFitness = sum(fx);
            if sumOfFitness == 0 % in case all individuals have the same, zero parentFitness
                fx = fx + 1/length(fx);
            else
                fx = fx ./ sum(fx); % norming to 1
            end
            for i=2:popSize      % building "inverse cdf"
                fx(i)=fx(i)+fx(i-1);
            end
            ind = zeros(1,offspringNumber);
            for i=1:offspringNumber % proportional (roulette) succession
                rnd=rand();
                b=1;
                e=popSize;
                while e-b>1 % searching in ordered array using biscetion
                    k=ceil(b+(e-b)/2);
                    if fx(k)<rnd
                        b=k;
                    else
                        e=k;
                    end
                end
                k=e;
                % we need to ensure that the first item
                % in the array can be selected as well
                if k<=2 && rnd < fx(1)
                    k=1;
                end
                ind(i)=k;
            end
            
            offspringPopulation=parentPopulation(:,ind);
            
        end
    end
end
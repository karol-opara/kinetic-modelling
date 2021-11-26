classdef Constraints
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function pop = Reflection(pop, lBounds, uBounds)
            % REFLECTION hanles constraint breach by reflection.
            % If single reflection is not enough, variable violating the
            % constraint is set randomly within feasible set.
            
            if(numel(lBounds) == 1 && numel(uBounds) == 1 ...
                    && lBounds == -Inf && uBounds == Inf)
                return;
            end
            
            % If bounds are infinite, return.
            if any(~isfinite(lBounds) | ~isfinite(uBounds))
                if any(isfinite(lBounds) | isfinite(uBounds))
                    warning('ConstraintsReflection:constrReflection',...
                        'Mix of finite and infinite bounds is not supported, no bounds assumed');
                end
                return;
            end
            
            [dim, popSize] = size(pop);
            L = repmat(lBounds, 1, popSize);
            U = repmat(uBounds, 1, popSize);
            %             lowerConstr = pop < L;
            %             pop = pop + lowerConstr.*(2*L - 2*pop);
            %             upperConstr = pop > U;
            %             pop = pop + upperConstr .* (2*U - 2*pop);
            lowerConstr = pop < L;
            pop(lowerConstr) = pop(lowerConstr) + (2*L(lowerConstr) - 2*pop(lowerConstr));
            upperConstr = pop > U;
            pop(upperConstr) = pop(upperConstr) + (2*U(upperConstr) - 2*pop(upperConstr));
            
            % checking in case of large violations of constraints
            largeViolation = pop < L | pop > U;
            if (any(any(largeViolation)))
                nonFeasibleRatio = sum(sum(largeViolation))/(dim*popSize);
                if(nonFeasibleRatio > 0.01)
                    warning('ConstraintsReflection:LargeViolation', ...
                        'Large violations of constraints, uniform random resampling used instead');
                end
                randomPopulation = InitializePopulation.Uniform(lBounds, uBounds, popSize);
                pop(largeViolation)= randomPopulation(largeViolation);
            end
        end % Reflection
        
    end
end
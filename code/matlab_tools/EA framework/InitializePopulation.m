classdef InitializePopulation
    methods (Static)
        function initialPopulation = Uniform(lBounds, uBounds, popSize, randState)
            % Uniform random initialization within given limits.
            % Returns population matrix with each column denoting a
            % single individual.
            
            if nargin == 4
                RandStream.getDefaultStream.State = randState;
            elseif nargin ~= 3
                error('Wrong number of arguments');
            end
            
            MinimizationProblem.CheckBounds(lBounds, uBounds);
            [dim, ~] = size(lBounds);
            
            % calculate and return initial population matrix
            initialPopulation = (repmat(lBounds, 1, popSize) + ...
                repmat((uBounds - lBounds), 1, popSize).*rand(dim, popSize));
            
        end % initPopulation
               
    end
end
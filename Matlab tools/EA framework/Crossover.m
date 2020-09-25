classdef Crossover
    methods (Static)
        function pop = Binomial(pop, pop2, crossoverProbability)
            % CROSSBINOMIAL Crossover two populations, substitute genes of individuals
            % from base population with individuals
            [dim, popSize] = size(pop);
            cross = rand(dim,popSize) < crossoverProbability;
            pop(cross) = pop2(cross);
        end % crossBinomial
    end
end
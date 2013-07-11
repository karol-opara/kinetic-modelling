classdef Mutation
    methods (Static)
        function pop = DeRandInfty(pop, F)
            [~,popSize] = size(pop);
            donors = pop(:,randi(popSize,1,popSize));
            covMtx = 2*F^2*cov(pop.');
            pop = Mutation.Gaussian(donors,covMtx);
        end
        
        function pop = Gaussian(pop, covMtx)
            % GAUSSIANMUTATION mutates each individual (column) of population pop
            % by adding to it a realization of a normally distributed variable with
            % zero mean and covariance matrix covMtx.
            % Alternatively sigma^2 (square of std) can be given as input.
            
            [dim, popSize] = size(pop);
            [r, c] = size(covMtx);
            % If covMtx is scalar make it a vector
            if (r == 1 && c == 1)
                covMtx = covMtx * ones(dim,1);
            end
            mu = zeros(dim,1);
            pop = pop + mvnrnd(mu, covMtx, popSize).';
        end % Gaussian
        
        function ui = DERand1(pop, F)
            % DERAND1 conduct DE/rand/1 mutation on population pop with scaling
            % factor F
            [dim popSize] = size(pop);
            if popSize < 3
                error('Too little number of vectors to conduct mutation in DE');
            end
            
            % Draw indices of points taken to calculate difference vectors
            ind2=randi(popSize,1,popSize);
            ind3=randi(popSize,1,popSize);
            
            % Ensure indices are pairwise different
            ind1 = 1:popSize;
            eq1=(ind1==ind2); % first ensure ind1 ~= ind 2
            while any(eq1)
                randInd = randi(popSize,sum(eq1),1);
                ind2(eq1)=randInd;
                eq1=(ind1==ind2);
            end
            eq2=(ind1==ind3 | ind2==ind3); % next ensure ind3 ~= ind2 and ind3~= ind1
            while any(eq2)
                randInd = randi(popSize,sum(eq2),1);
                ind3(eq2)=randInd;
                eq2=(ind1==ind3 | ind2==ind3);
            end
            
            % Perform differential mutation
            if(dim == 1 || popSize == 1)
                % One-dimensional matrices must be handled separately
                ui = pop + F*(pop(ind2) - pop(ind3));
            else
                ui = pop + F*(pop(:,ind2) - pop(:,ind3));
            end
        end % DERand1
    end
end
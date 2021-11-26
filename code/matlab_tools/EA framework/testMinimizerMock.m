classdef testMinimizerMock < Minimizer
    methods (Static)
        function [algorithmOpts] = PrintAlgorithmDefaultOptions()
            algorithmOpts = struct('Name','''Test algorithm''','Version','''0.0.0.0''');
        end
    end
    methods
        function InitializeMinimizer(this)
            this.MinProb = MinimizationProblem.GetTestProblem();
            this.Pop = [1.2 2 3 2
                        5   2 3 4.5];
        end
        
        function PerformIteration(this, fVal)
        end
    end
end
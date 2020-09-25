classdef DEClass < Minimizer
    
    properties
        IterationNo = 0;
    end
    
    methods (Static)
        function [algorithmOpts] = PrintAlgorithmDefaultOptions()
            algorithmOpts.Name =     '''DE''';
            algorithmOpts.Version =  '''0.0.0.1''';
            algorithmOpts.LastModificationDate = '''17-Apr-2012''';
            algorithmOpts.PopSize =           '1e3  % aka Np or mu';
            algorithmOpts.ScalingFactor =     '0.4    % aka F';
            algorithmOpts.RandState =         'RandStream.getDefaultStream.State % default rand state';
            algorithmOpts.ConstraintHandling ='''reflection'' % constraint handling technique';
        end
        
        function RunDERand1Bin(fitnessFunction, initLBounds, initUBounds)
            minProb = MinimizationProblem(fitnessFunction, initLBounds, initUBounds);
            DERand1 = DE(minProb,struct());
            DERand1.RunMinimization();
        end
    end
    
    methods
        function obj = DEClass(minProblem, inOpts)
            obj = obj@Minimizer(minProblem, inOpts);
            obj.InitializeMinimizer();
        end
        
        function InitializeMinimizer(this)
            this.AlgOpts = DE.PrintAlgorithmDefaultOptions();
            this.SetInitialPopulation(InitializePopulation.Uniform...
                (this.MinProb.InitLBounds, this.MinProb.InitUBounds, eval(this.AlgOpts.PopSize)));
            end
        
        function PerformIteration(this, fVal)
            this.IterationNo = this.IterationNo+1;
            
            mutatnts = Mutation.DERand1(this.Pop, AlgOpts.ScalingFactor);
            [fValmut] = this.MinProb.EvaluateFitnessFunction(mutants);
            
            % parent - offspring selection
            this.Pop(fValmut < fVal,:) = mutatnts(fValmut < fVal,:);
            
            this.Pop = Constraints.Reflection(this.Pop, this.MinProb.LBounds, this.MinProb.UBounds);
        end
        
    end
end
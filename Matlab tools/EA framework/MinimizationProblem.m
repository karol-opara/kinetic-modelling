classdef MinimizationProblem < handle
    properties
        Name;
        FitnessFunction
        EvalParallel;
        Dimension
        InitLBounds
        InitUBounds
        LBounds;
        UBounds;
    end
    
    methods
        function obj = MinimizationProblem(fitnessFunction, initLBounds, initUBounds, minProbOpts)
            if (nargin == 4)
                obj.InitializeOptions(minProbOpts);
            else
                obj.InitializeOptions();
            end
            obj.FitnessFunction = fitnessFunction;
            
            MinimizationProblem.CheckBounds(initLBounds, initUBounds)
            if (any(any(isfinite([initLBounds initUBounds])== false)))
                error('Initialization range bounds must be finite');
            end
            obj.InitLBounds = initLBounds;
            obj.InitUBounds = initUBounds;
            obj.Dimension = length(initLBounds);
            
            MinimizationProblem.CheckBounds(obj.LBounds, obj.UBounds);
            if (islogical(obj.EvalParallel) == false)
                error('Parameter EvalParallel must be logical (boolean)');
            end
        end
        
        function InitializeOptions(this, varargin)
            if (nargin == 2)
                minProbOpts = varargin{1};
            else
                minProbOpts = struct();
            end
            % First get default options
            optStr = MinimizationProblem.PrintDefaultOpts();
            % Convert default opts from strings to values
            opts = Utils.ParseDefaultOptions(optStr);
            % Add user options
            opts = Utils.SetUserOptions(opts, minProbOpts);
            % Save options as class properties
            fields = fieldnames(opts);
            for i = 1:length(fields)
                this.(fields{i}) = opts.(fields{i});
            end
        end
        
        function [fVal fMin xMin funEvals] = EvaluateFitnessFunction(this, pop, varargin)
            [~, popSize] = size(pop);
            if (this.EvalParallel)
                fVal = feval(this.FitnessFunction, pop, varargin{:});
            else
                fVal = zeros(1, popSize)*NaN;
                if (numel(varargin)>0)
                    additionalArgs = varargin{:};
                    for i=1:popSize
                        fVal(i) = feval(this.FitnessFunction, pop(:,i), additionalArgs{:,i});
                    end
                else
                    for i=1:popSize
                        fVal(i) = feval(this.FitnessFunction, pop(:,i));
                    end
                end
            end
            funEvals = popSize;
            
            % finding the best individual
            [fMin, indMin] = min(fVal);
            xMin = pop(:,indMin);
        end % evaluateFitnessFunction
        
       
    end
    
    methods (Static)         
        function problemOpts = PrintDefaultOpts()
            problemOpts.Name =              '''''   % (empty string)';
            problemOpts.EvalParallel =      'true   % objective function accepts NxM matrices, with M>1';
            problemOpts.LBounds =           '-Inf   % lower bounds, scalar or Nx1-vector';
            problemOpts.UBounds =           'Inf    % upper bounds, scalar or Nx1-vector';
        end

        function CheckBounds(lBounds, uBounds)
            [rl cl] = size(lBounds);
            [ru cu] = size(uBounds);
            if (rl ~= ru)
                error('Lower and upper bounds must have the same dimension');
            end
            if(cl ~= 1 || cu ~= 1)
                error('Upper and lower bounds must be given as column vectors');
            end
            if (any(any(lBounds > uBounds)))
                error('Lower bounds cannot be greater than upper ones');
            end
        end
        
        function minProb = GetTestProblem()
           A = [9 4
                4 3];
            b = [2 3];
            pOpts.Name = 'MinimizationProblem.GetTestProblem (ellipse)';
            pOpts.EvalParallel = false;
            minProb = MinimizationProblem(@(x) x.'*A*x+b*x, [1 3].', [1 5].', pOpts); 
        end
    end
    
    
end
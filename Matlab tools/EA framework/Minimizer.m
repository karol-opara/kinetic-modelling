classdef Minimizer < handle
    properties (Access = public)
        AlgOpts;
        MinProb; % minimization problem
        Output;
        StopCond;
        Log;
    end
    properties %(Access = protected)
        Pop;        % current population matrix, each individual encoded as one column
    end
    properties (Access = private)
        StartTime;    % timer used to measure computation time
    end
    
    methods (Static)
        function [stoppingCriteria output] = PrintDefaultOptions()
            % Returns a set of default options
            
            % stopping criteria
             stoppingCriteria.MaxFunEvals =       '1e4*dim  % maximal number of fitness function evaluations';
             stoppingCriteria.MaxIter =           '1e4*dim  % maximal iteration number';
            % stoppingCriteria.TolX =              'NaN'; % Not implemented yet
            % stoppingCriteria.TolFun =            'NaN'; % Not implemented yet
            
            % Implementation options
            
            output.Iter =                '0;	% number of iterations done within current run';
            output.FunEvals =            '0;    % number of fitness function evaluations';
            output.XMin =                'NaN;  % coordinates of the minimum';
            output.FMin =                'Inf;  % value of fitness function in the minimum';
            output.AlgOpts =    'NaN;  % algorithm options structure';
            output.MinProb = 'NaN;  % miminimzation problem';
            output.StopCond =    'NaN;  % stopping criteria structure';
            output.StopReason =          '''(unknown)'';        % stopping reason';
            
            
            %logging.NormPopCov
        end % printDefaultOptions
    end
    
    methods
        function obj = Minimizer(minProblem, inOpts)
            if (nargin==0)
                return;
            end
            obj.MinProb = minProblem;
            obj.SetDefaultOptions(minProblem);
            if (nargin > 1)
                obj.SetUserOptions(inOpts);
            end
        end
        
        function SetInitialPopulation(this,initPop)
            % sets few first initial population vectors
            
            % first, check arguments
            if (nargin ~= 2)
                error('Wrong number of arguments');
            end
            if (~isnumeric(initPop))
                error('Initial population must be numeric');
            end
            
            [initDim initSize] = size(initPop);
            if(initDim ~= this.MinProb.Dimension)
                error('Wrong dimension of initial population');
            end
            if(~isfinite(all(this.MinProb.LBounds)))
                lowerBoundMatrix = repmat(this.MinProb.LBounds(:,1),1,initSize);
                if (any(any((initPop < lowerBoundMatrix))))
                    error('Initial population is out of the feasible set');
                end
            end
            if(~isfinite(all(this.MinProb.UBounds)))
                upperBoundMatrix = repmat(this.MinProb.UBounds(:,2),1,initSize);
                if (any(any(initPop > upperBoundMatrix)))
                    error('Initial population is out of the feasible set');
                end
            end
            
            % second, set initial population
            this.Pop(:,1:initSize) = initPop;
        end
        
        function [xMin fMin output] = RunMinimization(this)
            this.InitializeOutput();
            
            this.InitializeMinimizer();
            while (this.IsStopConditionMet() == false)
                % Evaluate fitness function
                [fVal fMin xMin funEvals] = ...
                    this.MinProb.EvaluateFitnessFunction(this.Pop);
                
                % Update output for further study
                this.UpdateOutput(funEvals, fMin, xMin);
                
                % Perform algorithm-specific iteration
                this.PerformIteration(fVal);
            end
            
            [xMin fMin output] = this.FinishOutput();
        end
        
        function [xMin fMin output] = FinishOutput(this)
            computationTimeSeconds = toc(this.StartTime);
            [hour minute seconds] = sec2hms(computationTimeSeconds);
            seconds = floor(seconds);
            miliseconds = round(1e3*mod(computationTimeSeconds,1));
            this.Output.ComputationTime = [num2str(hour) ':'...
                num2str(minute) ':' num2str(seconds) '.'...
                num2str(miliseconds)];
            %this.Output.ComputationTimeSeconds = computationTimeSeconds;
            
            xMin = this.Output.XMin;
            fMin = this.Output.FMin;
            output = this.Output;
            
            
        end
        
        function InitializeOutput(this)
            this.Output.AlgOpts = this.AlgOpts;
            this.Output.MinProb = this.MinProb;
            this.Output.StopCond = this.StopCond;
            this.Output.StartTime = datestr(this.StartTime,'yyyy.mm.dd-HHMMSS');
            this.StartTime = tic;
        end
        
        function UpdateOutput(this, funEvals, fMin, xMin)
            if(fMin < this.Output.FMin)
                this.Output.XMin = xMin;
                this.Output.FMin = fMin;
            end
            this.Output.FunEvals = this.Output.FunEvals + funEvals;
            this.Output.Iter = this.Output.Iter + 1;
        end
        
        
        function [stop] = IsStopConditionMet(this)
            stop = false;
            if(this.Output.FunEvals >=this.StopCond.MaxFunEvals)
                stop = true;
                this.Output.StopReason = 'MaxFunEvals exceeded';
            end
            if (this.Output.Iter >= this.StopCond.MaxIter)
                stop = true;
                this.Output.StopReason = 'MaxIter exceeded';
            end
        end % IsStopConditionMet
        
        function SetDefaultOptions(this, minimizationProblem)
            [stopCondStr outputStr] = Minimizer.PrintDefaultOptions();
            this.StopCond = Utils.ParseDefaultOptions(stopCondStr, minimizationProblem);
            this.Output = Utils.ParseDefaultOptions(outputStr);
            [algorithmOptsStr] = this.PrintAlgorithmDefaultOptions();
            this.AlgOpts = Utils.ParseDefaultOptions(algorithmOptsStr, minimizationProblem);
        end
        
        
        function warnings = SetUserOptions(this, inOpts)
            % SETUSEROPTIONS overwrites default option values with values
            % provided by user.
            % Structure inOpts contains up to 3 substructures, namely:
            % AlgOpts, Output and StopCond,
            % which are documented in
            % PrintDefault[Algorithm]Options function
            
            warnings = false;
            optionStruct = fieldnames(inOpts);
            for i=1:length(optionStruct) % first find substructure names
                switch(optionStruct{i})
                    case {'AlgOpts', 'Output', 'StopCond'}
                        [this.(optionStruct{i}) warnings] = ...
                            Utils.SetUserOptions...
                            (this.(optionStruct{i}), inOpts.(optionStruct{i}));
                    otherwise
                        warning('Minimizer:SetUserOptions:Structure',...
                            ['Unrecognized structure name: '...
                            (optionStruct{i}) ', option discarded']);
                        warnings = true;
                end
            end
        end
        
        function value = getPopSize(this)
            [~, value] = size(this.Pop);
        end
        
        function value = getPop(this)
            value = this.Pop;
        end
    end
    
    methods (Abstract, Static)
        [algorithmOpts] = PrintAlgorithmDefaultOptions()
        
    end
    methods (Abstract)
        InitializeMinimizer(this)
        PerformIteration(this, fVal)
    end
    
    
end
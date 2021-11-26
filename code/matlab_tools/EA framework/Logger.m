classdef Logger < handle
    properties
        Iter
        FunEvals
        Pop
        PopFVal
        XMin
        FMin  
        Output
        Filename
    end
    
    properties (Access = private)
       IterToLog
       NextIterToLogInd
       NumberOfIterLogged 
       SaveFolder
       SaveFileName
       StartTime
    end
    
    methods (Static)
        function logging = PrintDefaultOpts()
            logging.IterToLog = ['[1 2 5 1e1 2e1 5e1 1e2 2e2 5e2 1e3 2e3 5e3 '...
                '1e4 2e4 5e4 1e5 2e5 5e5 1e6 2e6 5e6 1e7 2e7 5e7 1e8 2e8 5e8]; '...
                '% logging after these iterations'];
            logging.NextIterToLogInd = '1; % index of the next iteration to log';
            logging.NumberOfIterLogged = '0';
            logging.SaveFolder = '''.''; % default save folder is current folder';
            logging.SaveFileName = '''alg''; % default file name';
            logging.StartTime = 'datestr(this.StartTime,''yyyy.mm.dd-HHMMSS'')';
        end
    end
    
    methods
        
        function obj = Logger(logOpts)
            % First get default and user options
            optStr = Logger.PrintDefaultOpts();
            defOpts = Utils.ParseDefaultOptions(optStr);
            opts = Utils.SetUserOptions(defOpts, logOpts);
            fields = fieldnames(opts);
            for i = 1:length(fields)
                obj.(fields{i}) = opts.(fields{i});
            end
            
            obj.InitializeLog();
        end
        
        function InitializeLog(this)
            this.Iter = NaN*zeros(1,length(this.IterToLog));
            this.NextIterToLogInd = 1;
            
            maxLogLength = length(this.IterToLog);
            this.Pop{ind} = cell(1,maxLogLength);
            this.PopFVal{ind} = cell(1,maxLogLength);
            this.XMin{ind} = cell(1,maxLogLength);
            this.FMin(ind) = NaN*zeros(2,maxLogLength);
            this.FunEvals(ind) = NaN*zeros(2,maxLogLength);
            
            this.Filename = [this.SaveFolder '/logMinimizerRun-' this.SaveFileName...
                '-' this.StartTime];
        end
        
        function LogIteration(this, minimizer, fVal)
            if (this.Output.Iter == this.Log.IterToLog(this.Log.NextIterToLogInd))
                this.Log.NextIterToLogInd = this.Log.NextIterToLogInd + 1;
                this.LogCurrentIteration(minimizer, fVal);
            end
        end
        
        function FinishLogging(this, output)
            this.LogCurrentIteration(minimizer, fVal);
            range = 1:this.NumberOfIterLogged;
            this.Pop = this.Pop{range};
            this.PopFVal = this.PopFVal{range};
            this.XMin = this.XMin{range};
            this.FMin = this.FMin(range);
            this.FunEvals = this.FunEvals(range);
            
            this.Output = output;
            
            save(this.Filename);
            
        end
    end
    
    methods (Access = private)
        function LogCurrentIteration(this, minimizer)
            ind = this.NumberOfIterLogged + 1;
            
            this.Pop{ind} = minimizer.Pop;
            this.PopFVal{ind} = fVal;
            this.XMin{ind} = minimizer.Output.XMin;
            this.FMin(ind) = minimizer.Output.FMin;
            this.FunEvals(ind) = minimizer.Output.FunEvals;
            this.IterLogged(ind) = minimizer.Outpt.Iter;
            
            this.NumberOfIterLogged = ind;
        end
    end
    
end
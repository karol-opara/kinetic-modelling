classdef Utils
    methods (Static)
        function [outOpts warnings] = SetUserOptions(defaultOpts, userOpts)
            if(isstruct(defaultOpts) == false || isstruct(userOpts) == false)
                error('Utils:SetUserOptions',...
                    'Both defaul and user options must be structures');
            end
            warnings = false;
            
            outOpts = defaultOpts;
            
            userFields = fieldnames(userOpts);
            for j = 1:length(userFields) % then find field names
                if(isfield(outOpts, userFields{j}))
                    % if name exists overwrite it with user
                    % provided setting
                    outOpts.(userFields{j}) = ...
                        eval(['userOpts.' (userFields{j})]);
                else
                    warning('Utils:SetUserOptions:Field',...
                        ['Unrecognized field name: '...
                        'userOpts.' userFields{j}...
                        ', setting discarded']);
                    warnings = true;
                end
            end
        end
        
        function opts = ParseDefaultOptions(optStr, minProblem)
            if (nargin == 2)
                % dimension is used as a parameter to some default values
                dim = minProblem.Dimension;
            end
            
            % Convert default opts from strings to values
            fields = fieldnames(optStr);
            try
                for i = 1:length(fields)
                    opts.(fields{i}) = eval(optStr.(fields{i}));
                end
            catch err
                warning('Utils:ParseDefaultOptions',['Failed to evaluate line: ' ...
                    'opts.' fields{i} ' = eval(optStr.' fields{i} ')']);
                str = [' Possible reason: string values should be given with ' ...
                    'triple apostroph \ne.g. algorithmOptions.Name = ''''''Simple GA'''''';\n'];
                fprintf(str);
            end
            
        end
    end
end
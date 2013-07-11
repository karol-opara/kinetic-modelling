classdef ViewPlot
    methods (Static)
        function [hl hiqr hl2 hiqr2 horig] = QuantileOverIterations(iter, fVal, col, iter2, fVal2, col2, originTau)
      
            if (nargin == 0)
                warning('ViewPlot:QuantileOverIterations', ...
                    'No input arguments, test arguments used instead');
                iter= (1:100).';
                a = hilb(100)+hilb(100).*rand(100);
                fVal = a(1:25,1:100);
            end
            if (nargin < 3)
                col = [0.1 0.1 0.1];
            end
            
            [hl hiqr]=ViewPlot.PlotSingleQuantileOverIterations(iter,fVal,col);
                        grid('on');
            %xlabel('Iterations');
           % ylabel('Fitness');
            %legend([hl hiqr qmin], 'Runs', 'IQR', 'Min, median, max');
            
            if nargin>3
                [hl2 hiqr2]=ViewPlot.PlotSingleQuantileOverIterations(iter2,fVal2,col2);
                %legend([hl hiqr qmin hl2 hiqr2 qmin2], 'Runs', 'IQR', 'Min, median, max');
            end
            
           horig = line(get(gca,'XLim'),[originTau, originTau]);
           set(horig,'Color', 'k', 'LineWidth', 2);
           uistack(horig,'bottom');

            %ViewPlot.Save(gcf, 'Clayton');
        end
        
        function [hiqr qmedian]=PlotSingleQuantileOverIterations(iter,fVal,col)
            [popSize, iterNo] = size(fVal);
            if(iscolumn(iter)== false)
                error('ViewPlot:QuantileOverIterations', 'Iter must be a column vector');
            end
            if (iterNo ~= length(iter))
                error('ViewPlot:QuantileOverIterations',...
                    'Number of rows in fVal must be the same as the number of loged iterations');
            end
            q = quantile(fVal,[0.0 .05 .25 .50 .75 .95 1.0]);
            
            lightCol = (col+0.8*[1 1 1])/2;
            lighterCol = (col+3*[1 1 1])/4;
            qCol = 0.8*col;
            
            %set(ax,'XScale','log','YScale','log');
            if(any(isnan(iter)))
                len = 1:length(iter);
                ind = len(isnan(iter));
                i1 = 1:(ind-1);
                i2 = (ind+1):length(iter);
                %             hiqr = jbfill(iter.', q(3,:), q(5,:),lightCol, lightCol, true,0.5);
                hiqr = jbfill(iter(i1).', q(2,i1), q(6,i1),lightCol, lightCol, true,0.5);
                hiqr = jbfill(iter(i2).', q(2,i2), q(6,i2),lightCol, lightCol, true,0.5);
            end
            
            lIter = repmat([iter; NaN], popSize, 1);
            lFVal = reshape([fVal NaN*zeros(popSize,1)].', popSize*(iterNo+1),1);
            %hl = line(lIter,lFVal);
            %set(hl, 'Color', lightCol);
            %set(gca, 'XScale', 'log');
            
            width = 2;
            %jbfill(iter.', q(2,:), q(6,:),col, col, true,0.1);
            
%             qmin = line(iter, q(1,:));
%             set(qmin, 'Color', qCol, 'LineWidth', width);
            qmedian = line(iter, q(4,:));
            set(qmedian, 'Color', qCol, 'LineWidth', width);
%             qmax = line(iter, q(7,:));
%             set(qmax, 'Color', qCol, 'LineWidth', width);
%             q25q75 = line([iter; NaN; iter], [q(3,:).'; NaN; q(5,:).']);
%             set(q25q75,'Color',qCol);
            
            %             qmean = line(iter, mean(fVal));
%             set(qmean, 'Color', qCol, 'LineStyle', '--', 'LineWidth', width);
        end
        
        function Save(figureNumber, filename)
            figure(figureNumber);
            set(figureNumber, 'PaperPositionMode', 'auto');
            set(figureNumber, 'color', 'none');
            export_fig([filename '.pdf']);
            eval(['print -dpng -r1200 ' filename '.png']);
            eval(['print ' filename '.fig']);
            eval(['print ' filename ' -dmeta']);
        end
    end
end
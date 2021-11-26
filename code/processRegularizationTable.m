function processRegularizationTable()
wd = pwd();
[lambdas, ~] = getRegularizatonCoefficients(wd);


pnorms = {'log', 0.5, 1, 2};
qnorms = {'identity','log', 0.5, 1, 2};
plen = length(pnorms);
qlen = length(qnorms);

tab = ['\\hline \n Inter-component loss $L_\\text{comp}$ & Single-component loss $L_\\text{time}$ & Minimal sum of losses $M$ & Best relative reg. coefficient $\\lambda_\\text{rel}$ & Regularization coefficient $\\lambda$ \\\\ \\hline \n'];
for j = 1:qlen
    fr = true;
    for i = 1:plen
        if (fr)
            firstRow = ['\\hline \n' '\\multirow{4}{=}{' getReadablePnorms(qnorms{j}) '}'];
            fr = false;
        else
            firstRow = '';
        end
        lambda = lambdas(i,j);
        [~, lossMult] = getRegularizatonCoefficients(NaN, i, j);
        tab = [tab firstRow ' & ' getReadablePnorms(pnorms{i}) ' & ' sprintf('%1.1f', lossMult) ' & ' sprintf('%1.2f', lambda) ' & ' sprintf('%1.2f', lambda * lossMult)  ' \\\\ \n'];
    end
end
tab = [tab '\\hline\n'];
sprintf(tab)
end

function pn = getReadablePnorms(pnorm)
switch(pnorm)
    case 0.5
        pn = 'Root';
    case 1
        pn = 'Absolute';
    case 2
        pn = 'Square';
    case 'log'
        pn = 'Log';
    case 'rel'
        pn = 'Relative';
    case 'identity'
        pn = 'Identity (one-step loss)';
    otherwise
        error('processBenchmarkingExperiment:getReadablePnorms', 'Unsupported pnorm');
end
end
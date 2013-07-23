function err = KineticsModelLandscapeConstrainedObjectiveFunctionK(k, zml, timeZ, z0,p,q,type,lambda,sumForward, sumBackward)
k = projectK(k, sumForward, sumBackward);

err = ObjectiveFunction(k, zml, timeZ, z0, p, q, type, lambda);
kPenalty = log1p(abs(sum(k)-sumForward-sumBackward));
if (kPenalty/err > 1e2)
    warning('KineticsModelLandscape:ConstrainedObjectiveFunction',['Large k penalty kPenalty/err = ' num2str(kPenalty/err)]);
end
err = err + 1e-4*kPenalty;
end


function k = projectK(k, sumForward, sumBackward)
ifw = [1 3 5];
ibk = [2 4 6];
k(ifw) = sumForward .* k(ifw)./sum(k(ifw));
k(ibk) = sumBackward .* k(ibk)./sum(k(ibk));
end
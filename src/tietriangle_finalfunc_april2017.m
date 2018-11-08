function [x, gamma_out] = tietriangle_finalfunc_april2017(z,tietri_gamma, guess)
gamma_out = zeros(size(guess));
gamma_out(z(:,2) == 0) = 0;
gamma_out(z(:,3) == 0) = 1;
calcvals = (z(:,2) ~= 0) & (z(:,3) ~= 0);
if any(calcvals)
    e = ones(sum(calcvals),1);
    A = spdiags(e, 0, sum(calcvals), sum(calcvals));
    options = optimoptions(@fsolve,'display','off','tolfun',1e-8,'tolx',1e-8,...
    'jacobpattern',A);
    guess(calcvals) = max(0,min(1,z(calcvals,2)./(z(calcvals,2) + z(calcvals,3))));
    gamma_fit = fsolve(@(y)crossw_tietri_april2017(z(calcvals,:),y,tietri_gamma),guess(calcvals),options);
    gamma_out(calcvals) = gamma_fit;
end
x = permute(tietri_gamma(gamma_out),[2 3 1]);

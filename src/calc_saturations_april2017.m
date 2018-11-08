function [xij, S, gamma_3p] = calc_saturations_april2017(z, fofs, rho, gamma_3p)
[xij, alpha, gamma_3p] = vectorize_flash_2017(z, fofs, gamma_3p);
error = 1e5;
alpha_squeeze = squeeze(alpha)';
S = alpha_squeeze;

while error>1e-10
    Snew = bsxfun(@times,bsxfun(@times,alpha_squeeze,S*rho),1./rho');
    Snew = bsxfun(@times,Snew,1./sum(Snew,2));
    diff = abs((Snew-S));
    error = max(diff(:));
    S = Snew;
end

end

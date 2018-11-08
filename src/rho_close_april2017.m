function [Residue, x_new, S_new, Z_new, gamma_3p] = rho_close_april2017(vd, G_old, H_old, S_old, Z_old, x_old, rho_list, s, rho_func, dt, dx, gamma_3p)
vd = [1;vd];
S_new = S_old;
x_new = x_old;
G_new = (G_old(2:end,:) - ...
    (dt/dx).*(repmat(vd(2:end),[1 4]).*H_old(2:end,:) - ...
    repmat(vd(1:end-1),[1 4]).*H_old(1:end-1,:)));

vd_original = vd;
G_out = G_new;
G_new = max([G_old(1,:);G_out],0);
Z_new = G_new./repmat(sum(G_new,2),[1 4]);

[is_in_old,Loc_in_old] = ismembertol(Z_new,Z_old,1e-20,'byrows',true);
[x_new(:,:,~is_in_old), S_new(~is_in_old,:), gamma_3p] = calc_saturations_april2017(Z_new(~is_in_old,:),s,rho_list, gamma_3p);
x_new(:,:,is_in_old) = max(0,x_old(:,:,Loc_in_old(is_in_old)));
S_new(is_in_old,:) = max(0,S_old(Loc_in_old(is_in_old),:));
Residue_tmp = abs(1 - sum(G_new, 2)./ rho_func(S_new));
Residue = Residue_tmp(2:end) + 1e3*sum((G_out.*(G_out<0)).^2,2);
Residue(vd(2:end)<0) = Residue(vd(2:end)<0) + 1e2*vd_original(vd(2:end)<0).^2;
end
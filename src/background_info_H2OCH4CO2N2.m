% Relevant information for fractional flow problem

% Molecular weights for relevant gases
M_N2 = 28.0134;
M_CH4 = 16.04;
M_H20 = 18.01;
M_CO2 = 44.01;


% Densities for relevant phases
rho_hyd = 52132.517;
rho_water = 55707.94;
rho_gas = 3491.27*3;


% Non-dimensionalize densities
rho_Dim = rho_gas;
rho_hyd_D = rho_hyd/rho_Dim;
rho_gas_D = rho_gas/rho_Dim;
rho_water_D = rho_water/rho_Dim;
rho_list = [rho_water_D;rho_gas_D;rho_hyd_D];

% Viscosities for flowing phases
mu_g = 2e-5;
mu_l = 1.31e-3;
M = mu_l/mu_g;
phi = 0.5;

% Residual saturation levels
sgr = 0.02;
slr = 0.1;

k_rg = @(sg,sh) (((sg./(1 - sh) - sgr)).^2).*((sg./(1-sh))>=sgr).*((sg./(1-sh))<=(1 - slr)) + ...
    (0).*((sg./(1-sh))<sgr) + (1).*((sg./(1-sh))>(1-slr));
k_rw = @(sg,sh) ((((1 - sg - sh)./(1 - sh) - slr)).^4).*((sg./(1-sh))>=sgr).*((sg./(1-sh))<=(1 - slr)) + ...
    (1).*((sg./(1-sh))<sgr) + (0).*((sg./(1-sh))>(1-slr));
f_g = @(sg,sh) (k_rg(sg,sh)./(k_rg(sg,sh) + k_rw(sg,sh)./M));
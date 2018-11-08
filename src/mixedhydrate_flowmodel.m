%% Dynamic flow solution to gas injection process
% Author: K. Darnell
% Date: Nov. 6, 2016
% Modifying existing
close all
clearvars -except s

% Import functions necessary for simulation
if ~exist('s','var')
    prepare_lookup_functions
end

% Left (zL) and right (zR) state at initial condition for injection
% z = [H2O,CH4,CO2,N2] mole fraction
zL = [0;0;0.2;0.8]; %zW,zM,zC,zN
zR = [0.9321;0.0679;0;0]; %zW,zM,zC,zN sh50saq50
% zR = [0.475;0.525;0;0]; %zW,zM,zC,zN sh50sv50
% zR = [1;0;0;0]

Nc = length(zR);
save_str = '_MMDDYYY_Nx300_8020_saq50sh50';

% Options for plotting (leftover from previous script)
vd_BC = 1;
plot_int = 1e2; saveint = 1e2;

%Setup Grid and what not using Hesse formulation
xmin = 0; xmax = 1.0; Nx = 300;
Grid = initiate_Grid(xmin,xmax,Nx);
Grid = build_grid(Grid);
[D,Grad,I] = build_ops(Grid);

%% Set up fractional flow aspect
background_info_H2OCH4CO2N2;
% The output will be ordered such that [S_aq,S_v,S_h];
%% Define important functions

G_func = @(x,rho,S) permute(sum(bsxfun(@times,x,permute(bsxfun(@times,S,rho'),[3 2 1])),2),[3 1 2]);
Flux = @(S) repmat(rho_list',[Nx 1]).*[1 - f_g(S(:,2),S(:,3)),f_g(S(:,2),S(:,3)),zeros(Nx,1)];%Flux takes f_g(sg,sh)
Hfunc = @(x,S) [sum(squeeze(x(1,:,:))'.*Flux(S),2),...
    sum(squeeze(x(2,:,:))'.*Flux(S),2),...
    sum(squeeze(x(3,:,:))'.*Flux(S),2),...
    sum(squeeze(x(4,:,:))'.*Flux(S),2)];
% Hfunc = @(x,S) bsxfun(@times, Hfunc_p(x,S),(1 - S(:,3)).^2);
rho_func = @(S) S*rho_list;

%%
vd = ones(Nx,1); vd(1) = vd_BC;
%%
[x, S] = calc_saturations_april2017([zL';zR'], s, rho_list, 0.5);
[xL, SL, xR, SR] = deal(x(:,:,1),S(1,:), x(:,:,2), S(2,:));
%%
GL = G_func(xL,rho_list,SL);
GR = G_func(xR,rho_list,SR);

xall = ones(4,3,Nx);
xall(:,:,1) = xL;
xall(:,:,2:end) = repmat(xR,[1,1,Nx-1]);
Z = [zL';repmat(zR',[Nx-1,1])];
S = [SL;repmat(SR,[Nx-1,1])];
G = [GL;repmat(GR,[Nx-1,1])];
%%
H = Hfunc(xall,S);
rho = rho_func(S);
cfl = 0.025;
dt = Grid.dx*cfl;
dx = Grid.dx;
%%
if exist('sim_results') == 0
   mkdir 'sim_results' 
end
H_old = H;
G_old = G;
Z_old = Z;
S_old = S;
x_old = xall;
BIG = 1e6;
TOL = 1e-6;
do_tern = false;
Z_store = [];
x_store = [];
S_store = [];


Z_init = Z;
S_init = S;
gamma_3p = 0.5;
frontcalc = 3;

for tt = 1:round(1/dt)*0.5
    error = BIG;
    G_new = G_old;
    x_new = x_old;
    S_new = S_old;
    Z_new = Z_old;
    count = 0;
    
    if tt==1
        save(['sim_results/','0',save_str],...
            'G_new', 'S_new', 'Z_new','x_new', 'vd', 'Grid', 'cfl', ...
            'G_old', 'S_old', 'Z_old', 'x_old', 'rho_list')
    end
    
    frontcalc = min(Nx,max(frontcalc,find(max(abs(Z_new - Z_init),[],2)>1e-5,1,'last')+1));
    if isempty(frontcalc)
        frontcalc = 3;
    end
    e = ones(frontcalc-1,1);
    A = spdiags([e e], -1:0, frontcalc-1, frontcalc-1);
    vdguess = max(0,(sum(G_old(2:frontcalc,:),2) - rho(2:frontcalc) + ...
                cfl.*vd(1:frontcalc-1).*sum(H_old(1:frontcalc-1,:),2))./...
                (cfl.*sum(H_old(2:frontcalc,:),2)));
    rho_close_vd = @(vd) rho_close_april2017(...
        vd, G_old(1:frontcalc,:), H_old(1:frontcalc,:), S_old(1:frontcalc,:), ...
        Z_old(1:frontcalc,:), x_old(:,:,1:frontcalc), ...
        rho_list, s, rho_func, dt, dx, gamma_3p);
    options = optimoptions(@fsolve, 'display','off','findifftype',...
                            'central', 'tolfun', 1e-6, 'tolx', 1e-6,...
                             'scaleproblem', 'jacobian','maxiter', 100,...
                             'algorithm', 'trust-region-reflective',...
                             'Jacobpattern', A);
    vd(2:frontcalc) = fsolve(rho_close_vd, vdguess, options);
    G_new(2:frontcalc,:) = (G_old(2:frontcalc,:) - ...
        (dt/dx).*(repmat(vd(2:frontcalc),[1 4]).*H_old(2:frontcalc,:) - ...
        repmat(vd(1:frontcalc-1),[1 4]).*H_old(1:frontcalc-1,:)));
    [~, x_tmp, S_tmp, Z_tmp] = rho_close_april2017(...
        vd(2:frontcalc),  G_old(1:frontcalc,:), H_old(1:frontcalc,:), S_old(1:frontcalc,:), ...
        Z_old(1:frontcalc,:), x_old(:,:,1:frontcalc), ...
        rho_list, s, rho_func, dt, dx, gamma_3p);
    
    x_new(:,:,1:frontcalc) = x_tmp;
    S_new(1:frontcalc,:) = S_tmp;
    Z_new(1:frontcalc,:) = Z_tmp;
    
    rho = rho_func(S_new);
    H_new = Hfunc(x_new,S_new);
    G_old = G_new;
    H_old = H_new;
    S_old = S_new;
    Z_old = Z_new;
    x_old = x_new;
    
    if ~mod(tt,saveint) || tt==0
        save(['sim_results/',num2str(tt),save_str],...
            'G_new', 'S_new', 'Z_new','x_new', 'vd', 'Grid', 'cfl', ...
            'G_old', 'S_old', 'Z_old', 'x_old', 'rho_list')
        [tt*dt 0 0 mean(S_new)]
    end
    MOC_plot_4Nc
end

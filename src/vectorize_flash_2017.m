% function [xij,S,alpha] = get_H2OCH4CO2N2_sats_2017(z,rho,fofs)
function [xij, alpha, gamma_3p] = vectorize_flash_2017(z, fofs, gamma_3p)
    
    % Define tolerance level for the calculation error
    TOL = 1e-5;
    
    % Normalize z in the event it is not normalized
    z = bsxfun(@times,z,1./sum(z, 2));
    z_T = z';

    
    % Potentially useful information
    Nx = length(z(:,1));
    Np = 3;
    Nc = length(z(1,:));
    alpha = zeros(Np, 1, Nx);
    alpha_best = zeros(Np, 1, Nx);
    alpha_tmp = zeros(Np, 1, Nx);
    xij = zeros(Nc, Np, Nx);
    xij_best = zeros(Nc, Np, Nx);
    error_best = 1e6*ones(Nx,1);
    error_tmp = zeros(Nx,1);
    
    % Useful indices for future calculations
    Hgt_ind = z(:,1) > 0.86;
    adj_min = fofs.min_n2*ones(Nx,1);
    n2_fraction = bsxfun(@times, z(:,4), 1./sum(z(:,2:4),2));
    Aq_pure = z(:,1) > fofs.max_h2o_Aq;
    V_pure = z(:,1) < fofs.min_h2o_V;
    n2_gt = n2_fraction > fofs.mean_n2;
    n2_AqV = n2_fraction > fofs.max_n2;
    n2_notAqV = n2_fraction < fofs.min_n2;
    adj_min(Hgt_ind) = fofs.min_n2;
    adj_min(~Hgt_ind) = fofs.min_n2;
    n2_between = (n2_fraction >= adj_min) & (n2_fraction <= fofs.max_n2);
    % Initialized
    phases_found = false(Nx,1);
    checked_3phase = false(Nx,1);
    checked_AqV = false(Nx,1);
    checked_AqH = false(Nx,1);
    checked_HV = false(Nx,1);
    checked_all = false(Nx,1);
    checked_pureV = false(Nx,1);
    checked_pureAq = false(Nx,1);
    H_negative = false(Nx,1);
    V_negative = false(Nx,1);
    Aq_negative = false(Nx,1);

    alpha_find = true(size(alpha_tmp));
    x_include = true(size(xij));
    twophase_negatives = false(size(alpha_tmp));
    
    % Dummy counter for now
    iter = 0;
    while any(~phases_found & ~checked_all) && iter < 5
        
        if any(Aq_pure & ~checked_pureAq)
           error_best(Aq_pure) = 0; 
           xij(:,1,Aq_pure) = z(Aq_pure,:)';
           xij_best(:,1,Aq_pure) = xij(:,1,Aq_pure);
           alpha(:,:,Aq_pure) = 0;
           alpha_best(:,:,Aq_pure) = 0;
           alpha(1,1,Aq_pure) = 1;
           alpha_best(1,1,Aq_pure) = 1;
           phases_found(Aq_pure) = true; 
           checked_pureAq(Aq_pure) = true;
        end
        
        if any(V_pure & ~checked_pureV)
           error_best(V_pure) = 0; 
           xij(:,2,V_pure) = z(V_pure,:)';
           xij_best(:,2,V_pure) = xij(:,2,V_pure);
           alpha(:,:,V_pure) = 0;
           alpha_best(:,:,V_pure) = 0;
           alpha(2,1,V_pure) = 1;
           alpha_best(2,1,V_pure) = 1;
           phases_found(V_pure) = true;
           checked_pureV(V_pure) = true;
        end
        
        check_AqV = ~phases_found & ...
            ((n2_AqV & ~checked_AqV) | (checked_3phase & H_negative)) ...
            & ~checked_AqV;
        check_HV = ~phases_found & ...
            ((n2_notAqV & ~checked_HV) | (checked_3phase & Aq_negative)) ...
            & ~Hgt_ind & ~checked_HV;
        check_AqH = ~phases_found & ...
            ((n2_notAqV & ~checked_AqH) | (checked_3phase & V_negative)) ...
            & Hgt_ind & ~checked_AqH;
        check_3phase = ~phases_found & ~checked_3phase & ...
            ((n2_AqV & checked_AqV) | (~n2_notAqV & (checked_AqH | checked_HV) | (n2_between)));
        checked_all = checked_AqV & (checked_AqH | checked_HV) & checked_3phase;
        
        if any(check_AqV & ~phases_found)
           [Aqcells, Vcells] = cellfun(fofs.AqV_func,num2cell(z(check_AqV, :),2),...
               'UniformOutput',false);
           xij(:,1,check_AqV) = cell2mat(Aqcells)';
           xij(:,2,check_AqV) = cell2mat(Vcells)';
           xij(:,3,check_AqV) = 0;
           x_include(:,[1,2],check_AqV) = true;
           x_include(:,3,check_AqV) = false;
           alpha_find(3,:,check_AqV) = false;
           alpha_find([1,2],:,check_AqV) = true;
           checked_AqV(check_AqV) = true;
        end
        
        if any(check_AqH & ~phases_found)
           [Aqcells, Hcells] = cellfun(fofs.AqH_func,num2cell(z(check_AqH, :),2),...
               'UniformOutput',false);
           xij(:,1,check_AqH) = cell2mat(Aqcells)';
           xij(:,3,check_AqH) = cell2mat(Hcells)';
           xij(:,2,check_AqH) = 0;
           x_include(:,[1,3],check_AqH) = true;
           x_include(:,2,check_AqH) = false;
           alpha_find(2,:,check_AqH) = false;
           alpha_find([1,3],:,check_AqH) = true;
           checked_AqH(check_AqH) = true;
        end
        
        if any(check_HV & ~phases_found)
           [Hcells, Vcells] = cellfun(fofs.HV_func,num2cell(z(check_HV, :),2),...
               'UniformOutput',false);
           xij(:,3,check_HV) = cell2mat(Hcells)';
           xij(:,2,check_HV) = cell2mat(Vcells)';
           xij(:,1,check_HV) = 0;
           x_include(:,[2,3],check_HV) = true;
           x_include(:,1,check_HV) = false;
           alpha_find(1,:,check_HV) = false;
           alpha_find([2,3],:,check_HV) = true;
           checked_HV(check_HV) = true;
        end
        
        if any(check_3phase & ~phases_found)
            [x_out, gamma_3p] = fofs.AqHV_func(...
                 z(check_3phase, :),gamma_3p*ones(sum(check_3phase),1));
             gamma_3p = 0.5;
             xij(:, :, check_3phase) = x_out;
            checked_3phase(check_3phase) = true;
            x_include(:,:,check_3phase) = true;
            alpha_find(:,:,check_3phase) = true;
        end
        
        alpha_tmp(:,:,~phases_found) = ...
            solve_alpha(xij(:,:,~phases_found), z_T(:, ~phases_found));
        
        twophase_negatives(:,:,~check_3phase) = alpha_tmp(:,:,~check_3phase)<0;
        if any(twophase_negatives(:))
           alpha_tmp(twophase_negatives) = 0;
           set_pure = squeeze(any(twophase_negatives));
           alpha_tmp(:,:,set_pure) = min(1,alpha_tmp(:,:,set_pure)); 
           xij(:,twophase_negatives) = 0;
        end
        
        error_tmp(~phases_found) = ...
            calc_comp_error(xij(:,:,~phases_found), alpha_tmp(:,:,~phases_found), z(~phases_found,:));

        
        negatives = alpha_tmp(:,:,check_3phase)<0;
        
        H_negative(check_3phase) = negatives(3,:,:);
        V_negative(check_3phase) = negatives(2,:,:);
        Aq_negative(check_3phase) = negatives(1,:,:);
        
        
        
        
        xij_best(:,:,~phases_found & (error_tmp < error_best)) = ...
            xij(:,:,~phases_found & (error_tmp < error_best));
        alpha_best(:,:,~phases_found & (error_tmp < error_best)) = ...
            alpha_tmp(:,:,~phases_found & (error_tmp < error_best));
        error_best(~phases_found & (error_tmp < error_best)) = ...
            error_tmp(~phases_found & (error_tmp < error_best));
        
        phases_found(~phases_found & (error_tmp < TOL)) = true;
        iter = iter + 1;
    end
xij = xij_best;
alpha = alpha_best;
    % Assuming xij is Nc x Np x Nx
    % Assuming z is Nx x Nc
    % Assuming alpha is Np x 1 x Nx
    function error = calc_comp_error(xij, alpha, z)
        calc_z = sum(bsxfun(@times,xij,permute(alpha, [2 1 3])),2);
        error = sum((squeeze(calc_z)' - z).^2,2);
        error(squeeze(any(alpha<0,1))) = 1e6;
        error(squeeze(any(alpha>1,1))) = 1e6;
        if any(error == 1e6)
            error = error;
        end
    end

    function alpha_out = solve_alpha(xij, z_T)
        % Not sure if this is the best way to do this, but it's much faster
        % than any other option.
        reduced_alpha_inds = alpha_find(:,:,~phases_found);
        alpha_intermed = zeros(3*length(z_T(1,:)),1);
        c = num2cell(xij,[1 2]);
        c = cellfun(@(x) sparse(x),c,'uni',0);
        B = blkdiag(c{:});
        alpha_direct = B(:,reduced_alpha_inds)\z_T(:);
        alpha_intermed(reduced_alpha_inds) = alpha_direct;
        alpha_out = reshape(alpha_intermed,3,1,[]);
        alpha_out = bsxfun(@times,alpha_out,1./sum(alpha_out,1));
    end
end
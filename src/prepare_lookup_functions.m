%% Synthesize lookup table info
%%%%% ------Do all of this only once-------%%%%%%
if ~exist('AqHV_x','var')
    load ../Data/277.15K70Bar_LookupTable_April2017.mat
    x_mat_full = x_mat;
    alpha_mat_full = alpha_mat;
    C_full = C;
    k_inds = squeeze(sum((abs(sum(x_mat)-1)<1e-3)|(abs(sum(x_mat))<1e-3))==4);
    for ii=1:length(x_mat(1,1,:))
        tmp = x_mat(C(ii,:)>0,:,ii)>1e-9;
        k_inds2(ii,1) = all((sum(tmp)==length(tmp(:,1)))|(sum(tmp)==0));
        k_inds3(ii,1) = ~any(any(x_mat(:,:,ii)==1));
    end
    k_inds = (k_inds & k_inds2) & k_inds3 ;
    keep_num = sum(k_inds);
    C_keep = C(k_inds,:);
    x_mat_keep = x_mat(:,:,k_inds);
    alpha_mat_keep = alpha_mat(:,:,k_inds);
    
    
    %Condense x_mat_keep, so that each section only has compositions of stable phases
    %''''
    % This is new starting on 11/29/16
    % x_mat_keep is originally ordered as 'H20','CH4','N2','CO2'
    % I will switch ordering to match the incoming z
    % x_matrix is now ordered 'H20','CH4','CO2','N2' and the columns are
    % 'Aq','V','L','H' (the rows of alpha_mat_keep are also ordered this
    % way).

    x_matrix = [x_mat_keep(1,:,:);x_mat_keep(2,:,:);x_mat_keep(4,:,:);x_mat_keep(3,:,:)];
    x_matrix = x_matrix(:,:,squeeze(sum(all(x_matrix==0)))<3);
    alpha_mat_keep = alpha_mat_keep(:,:,squeeze(sum(all(x_matrix==0)))<3);


    % Pull out appropriate indices
    AqV_inds = squeeze((alpha_mat_keep(1,:,:)>0)&(alpha_mat_keep(2,:,:)>0)&(alpha_mat_keep(3,:,:)<1e-4));
    HV_inds = squeeze((alpha_mat_keep(2,:,:)>0)&(alpha_mat_keep(3,:,:)<1e-3)&(alpha_mat_keep(4,:,:)>0));
    AqH_inds = squeeze((alpha_mat_keep(1,:,:)>0)&(alpha_mat_keep(3,:,:)<1e-3)&(alpha_mat_keep(4,:,:)>0));
    LH_inds = squeeze((alpha_mat_keep(1,:,:)<1e-3)&(alpha_mat_keep(3,:,:)>0)&(alpha_mat_keep(4,:,:)>0));
    AqHV_inds = squeeze((alpha_mat_keep(1,:,:)>1e-3)&(alpha_mat_keep(2,:,:)>1e-3)&(alpha_mat_keep(3,:,:)==0)&(alpha_mat_keep(4,:,:)>1e-3))&squeeze(x_matrix(4,1,:)~=0)&squeeze(x_matrix(4,4,:)>1e-3);
    LHV_inds = squeeze((alpha_mat_keep(1,:,:)<1e-3)&(alpha_mat_keep(2,:,:)>0)&(alpha_mat_keep(3,:,:)>0)&(alpha_mat_keep(4,:,:)>0));
    % Now two columns where the first column is Aq and second column is V

    % Single-guest indices
    %     n2single = squeeze(all(x_matrix([2:3],2,:) == 0,1));
    co2single = squeeze(all(x_matrix([2,4],4,:) == 0,1)) & squeeze(all(x_matrix([1,3],1,:) ~= 0,1));
    ch4single = squeeze(all(x_matrix([3:4],4,:) == 0,1));
    co2liquid = squeeze(x_matrix(3,2,:))>0.85;

    % For the two-phase zones, just keep tie lines
    % Make tessellations of the space and setup function to check if point
    % lies in tessellation.
    AqV_tielines = x_matrix(:,1:2,[AqV_inds | AqHV_inds]);
    AqV_tielines = AqV_tielines(:,:,~any(AqV_tielines(1,:,:)==0));

    HV_tielines =  x_matrix(:,[2,4],[(HV_inds | ch4single | AqHV_inds) & ~co2liquid]);
    HV_tielines = HV_tielines(:,:,~any(HV_tielines(1,:,:)==0));

    AqH_tielines = x_matrix(:,[1,4],[AqH_inds | co2single | ch4single | AqHV_inds]);
    AqH_tielines = AqH_tielines(:,:,~any(AqH_tielines(1,:,:)==0));

    % Setup tessellation for single phase zones
    Aq_tess_inds = AqH_inds | AqV_inds | AqHV_inds;
    V_tess_inds = AqV_inds | HV_inds | AqHV_inds;

    AqHV_tietriangles = x_matrix(:,[1,2,4],AqHV_inds);
    AqHV_tietriangles = AqHV_tietriangles(:,:,~any(AqHV_tietriangles(1,:,:)==0));

    AqHV_Aq = squeeze(AqHV_tietriangles(:,1,:))'; %Aq composition of AqHV tietriangles
    AqHV_H = squeeze(AqHV_tietriangles(:,3,:))'; %H composition of AqHV tietriangles
    AqHV_V = squeeze(AqHV_tietriangles(:,2,:))'; %V composition of AqHV tietriangles
    [~,tietri_sort] = sort(AqHV_V(:,4));
    AqHV_Aq = AqHV_Aq(tietri_sort,:);
    AqHV_H = AqHV_H(tietri_sort,:);
    AqHV_V = AqHV_V(tietri_sort,:);
    AqHV_unq = uniquetol([AqHV_Aq AqHV_V AqHV_H],1e-5,'byrows',true);

    AqHV_Aq_unq = AqHV_unq(:,1:4);
    AqHV_V_unq = AqHV_unq(:,5:8);
    AqHV_H_unq = AqHV_unq(:,9:12);

    [~,~,fof_gammaAqHV_Aq] = interparc(0.5,AqHV_Aq_unq(:,1),AqHV_Aq_unq(:,2),AqHV_Aq_unq(:,3),AqHV_Aq_unq(:,4));
    [~,~,fof_gammaAqHV_H] = interparc(0.5,AqHV_H_unq(:,1),AqHV_H_unq(:,2),AqHV_H_unq(:,3),AqHV_H_unq(:,4));
    [~,~,fof_gammaAqHV_V] = interparc(0.5,AqHV_V_unq(:,1),AqHV_V_unq(:,2),AqHV_V_unq(:,3),AqHV_V_unq(:,4));


    q_points = linspace(0,1,1e3)';
    Aq_pts_out = abs(fof_gammaAqHV_Aq(q_points));
    V_pts_out = abs(fof_gammaAqHV_V(q_points));
    H_pts_out = abs(fof_gammaAqHV_H(q_points));

    pp_Aq1 = interp1(q_points,Aq_pts_out(:,1),'pchip','pp');
    pp_Aq2 = interp1(q_points,Aq_pts_out(:,2),'pchip','pp');
    pp_Aq3 = interp1(q_points,Aq_pts_out(:,3),'pchip','pp');
    pp_Aq4 = interp1(q_points,Aq_pts_out(:,4),'pchip','pp');
    new_Aqfunc = @(x)[ppval(pp_Aq1,x),...
                      ppval(pp_Aq2,x)...
                      ppval(pp_Aq3,x),...
                      ppval(pp_Aq4,x)];
    pp_V1 = interp1(q_points,V_pts_out(:,1),'pchip','pp');
    pp_V2 = interp1(q_points,V_pts_out(:,2),'pchip','pp');
    pp_V3 = interp1(q_points,V_pts_out(:,3),'pchip','pp');
    pp_V4 = interp1(q_points,V_pts_out(:,4),'pchip','pp');
    new_Vfunc = @(x)[ppval(pp_V1,x),...
                      ppval(pp_V2,x)...
                      ppval(pp_V3,x),...
                      ppval(pp_V4,x)];
    pp_H1 = interp1(q_points,H_pts_out(:,1),'pchip','pp');
    pp_H2 = interp1(q_points,H_pts_out(:,2),'pchip','pp');
    pp_H3 = interp1(q_points,H_pts_out(:,3),'pchip','pp');
    pp_H4 = interp1(q_points,H_pts_out(:,4),'pchip','pp');
    new_Hfunc = @(x)[ppval(pp_H1,x),...
                      ppval(pp_H2,x)...
                      ppval(pp_H3,x),...
                      ppval(pp_H4,x)];
    %Function returns the coordinates of all tie triangle at some
    %fractional distance along vertices in tie triangle phase path
    tietri_gamma = @(gamma) permute(reshape([new_Aqfunc(gamma);new_Vfunc(gamma);new_Hfunc(gamma)],[length(gamma) 3 4]),[1 3 2]);



    AqV_Aq = squeeze(AqV_tielines(:,1,:))'; %Aq composition of AqV tielines
    AqV_V = squeeze(AqV_tielines(:,2,:))'; %V composition of AqV tielines
    AqV_Aq = [AqV_Aq;Aq_pts_out];
    AqV_V = [AqV_V;V_pts_out];

    AqV_hullpts = [AqV_Aq;AqV_V];
    AqV_tess = convhulln(AqV_hullpts(:,1:3));
    [tol,AqV_aN,AqV_nrmls] = inhullpre(AqV_hullpts(:,1:3),AqV_tess);
    AqV_check = @(z) inhullpost(z(:,1:3),AqV_nrmls,AqV_aN,1e-8);

    HV_V = squeeze(HV_tielines(:,1,:))'; %V composition of HV tielines
    HV_H = squeeze(HV_tielines(:,2,:))'; %H composition of HV tielines
    HV_H = [HV_H;H_pts_out];
    HV_V = [HV_V;V_pts_out];

    HV_hullpts = [HV_H;HV_V];
    HV_tess = convhulln(HV_hullpts(:,1:3));
    [tol,HV_aN,HV_nrmls] = inhullpre(HV_hullpts(:,1:3),HV_tess);
    HV_check = @(z) inhullpost(z(:,1:3),HV_nrmls,HV_aN,1e-8);

    AqH_Aq = squeeze(AqH_tielines(:,1,:))'; %Aq composition of AqH tielines
    AqH_H = squeeze(AqH_tielines(:,2,:))'; %H composition of AqH tielines
    AqH_Aq = [AqH_Aq;Aq_pts_out];
    AqH_H = [AqH_H;H_pts_out];

    AqH_hullpts = [AqH_Aq;AqH_H];
    AqH_tess = convhulln(AqH_hullpts(:,1:3));
    [tol,AqH_aN,AqH_nrmls] = inhullpre(AqH_hullpts(:,1:3),AqH_tess);
    AqH_check = @(z) inhullpost(z(:,1:3),AqH_nrmls,AqH_aN,1e-8);

    % Each tie line can be uniquely defined by its midpoint
    AqV_midpt = squeeze(0.5*(AqV_V + AqV_Aq));
    HV_midpt = squeeze(0.5*(HV_V + HV_H));
    AqH_midpt = squeeze(0.5*(AqH_Aq + AqH_H));
    AqHV_mdpt = squeeze(mean(AqHV_tietriangles(:,:,tietri_sort),2))';
    AqHV_midpt_dist = @(z) sqrt(sum(bsxfun(@minus,AqHV_mdpt,z).^2,2));


    %Ordered such that the columns are 'K_{H2O}','K_{CH4}','K_'{CO2}','K_{N2}' and rows are each tie line
    % Not sure what to do with NaNs yet.....
    AqV_K_VoverAq = squeeze(AqV_V./AqV_Aq)';
    HV_K_VoverH = squeeze(HV_V./HV_H)';
    AqH_K_HoverAq = squeeze(AqH_H./AqH_Aq)';
    AqHV_K_HoverAq = squeeze(AqHV_H./AqHV_Aq)';
    AqHV_K_VoverAq = squeeze(AqHV_V./AqHV_Aq)';

   
    dist_toz = @(z,P1,P2) (z - P1);
    dist_f = @(z,P1,P2) (P2 - P1);

    factor = @(z,P1,P2) (sum(dist_f(z,P1,P2).*dist_toz(z,P1,P2),2)./...
        sum(dist_f(z,P1,P2).*dist_f(z,P1,P2),2));
    dist_single = @(z,P1,P2) sqrt(sum((repmat(z,[length(P1(:,1)) 1]) - P1 - ...
        repmat(factor(repmat(z,[length(P1(:,1)) 1]),P1,P2),[1 4]).*dist_f(repmat(z,[length(P1(:,1)) 1]),P1,P2)).^2,2));

    % Pass all tie lines into function
    AqV_dist_full = @(z) dist_single(z,AqV_Aq,AqV_V);
    HV_dist_full = @(z) dist_single(z,HV_H,HV_V);
    AqH_dist_full = @(z) dist_single(z,AqH_H,AqH_Aq);

    % Output functions
    AqV_func = @(z) tieline_finalfunc_april2017(AqV_Aq,AqV_V,z,dist_single);
    HV_func = @(z) tieline_finalfunc_april2017(HV_H,HV_V,z,dist_single);
    AqH_func = @(z) tieline_finalfunc_april2017(AqH_Aq,AqH_H,z,dist_single);
    AqHV_func = @(z, guess) tietriangle_finalfunc_april2017(z, tietri_gamma, guess);


    mean_n2 = mean(bsxfun(@times,AqHV_V(:, 4),1./sum(AqHV_V(:,2:4),2)));
    max_n2 = max(bsxfun(@times,AqHV_V(:,4),1./sum(AqHV_V(:,2:4),2)));
    min_n2 = min(bsxfun(@times,AqHV_H(:,4),1./sum(AqHV_H(:,2:4),2)));
    mean_h2o_Aq = mean([AqV_Aq(:,1);AqH_Aq(:,1)]);
    max_h2o_Aq = max([AqV_Aq(:,1);AqH_Aq(:,1)]);
    mean_h2o_V = mean([AqV_V(:,1);HV_V(:,1)]);
    min_h2o_V = min([AqV_V(:,1);HV_V(:,1)]);

    clearvars -except AqV_func HV_func AqH_func  AqHV_func mean_n2 max_n2 min_n2 mean_h2o_Aq max_h2o_Aq mean_h2o_V min_h2o_V AqV_check  HV_check AqH_check
    s = v2struct;
    clearvars -except s
end
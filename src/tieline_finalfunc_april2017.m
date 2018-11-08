function [z_P1,z_P2] = tieline_finalfunc_april2017(P1,P2,z,dist_func)

% Make sure z adds to one and there are no negative values.
z(z<0)=0;
z = z./sum(z);
alpha = zeros(2,1);

P_out = uniquetol([P1 P2],1e-7,'byrows',true);
P1 = P_out(:,1:4);
P2 = P_out(:,5:8);
masked_nonzero = all(P1(:,z==0)==0,2);

P1 = P1(masked_nonzero,:);
P2 = P2(masked_nonzero,:);
Ks = P1./P2;

distances1 = max(1e-20,sqrt(sum((bsxfun(@minus,z,P1) - ...
    bsxfun(@times,sum(bsxfun(@minus,z,P1).* ...
    (P2 - P1),2)./sum((P2 - P1).^2,2),P2-P1)).^2,2)));

distances2 = max(1e-20,sqrt(sum((bsxfun(@minus,z(2:end),P1(:,2:end)) - ...
    bsxfun(@times,sum(bsxfun(@minus,z(:,1),P1(:,1)).* ...
    (P2(:,1) - P1(:,1)),2)./sum((P2(:,1) - P1(:,1)).^2,2),P2(:,2:end)-P1(:,2:end))).^2,2)));

z_resides = false;
count = 0;
errornow = 1e6;
while ~z_resides || errornow > min(z(z>0)) || errornow > 1e-8 || count < 1
    
    if count == 0
        distances = distances1;
        ntake_param = 10;
        p = 3;
    elseif count == 1 
        distances = distances2;
        errorold = errornow;
        z_P1tmp = z_P1;
        z_P2tmp = z_P2;
        K_1over2tmp = K_1over2;
        ntake_param = 5;
        p = 4;
    else
        K_1over2 = (K_1over2/errornow + K_1over2tmp/errorold)./(1/errornow + 1/errorold);
        betaerror = 1e6;
        Res = 1e6;
        beta = alpha_tmp(1);
        while betaerror > 1e-8 && abs(Res) > 1e-8 
            Res = sum(z.*(K_1over2 - 1)./(1 + beta.*(K_1over2 - 1)));
            JRes_beta = -sum(z.*(K_1over2 - 1).^2./(1 + beta.*(K_1over2 -1 )).^2);
            betanew = beta - Res/JRes_beta;
            betaerror = abs(betanew - beta)/beta;
            beta = betanew;
        end
        
        z_P2 = z./(1 +(K_1over2 - 1)*beta);
        z_P1 = K_1over2.*z_P2;
        
        z_resides = all((z - z_P1).*(z - z_P2) < 0 | z==0);    
        if ~z_resides 
            z_P1 = z_P1tmp;
            z_P2 = z_P2tmp;
        end
        break
    end
    [dist_B,I] = sort(distances);
    P1_mod = P1(I,:);
    P2_mod = P2(I,:);
   
    if length(I)>1
        ntake = min(length(I),ntake_param);
        z_P1 = sum(bsxfun(@times,P1_mod(1:ntake,:),1./dist_B(1:ntake).^p))/(sum(1./dist_B(1:ntake).^p));
        z_P2 = sum(bsxfun(@times,P2_mod(1:ntake,:),1./dist_B(1:ntake).^p))/(sum(1./dist_B(1:ntake).^p));
    elseif ~isempty(I)
        z_P1 = P1_mod(1,:);
        z_P2 = P2_mod(1,:);
    else
        z_P1 = [-1,-1,-1,-1]/4;
        z_P2 = [1,1,1,1]/4;
    end
    
    K_1over2 = z_P1./z_P2;
    K_1over2(isnan(K_1over2)) = 0;
    alpha_tmp = ([[z_P1;z_P2]';[1 1]]\[z';1]);
    alpha = alpha_tmp./sum(alpha_tmp);
    if count <=1
        z_resides = all((z - z_P1).*(z - z_P2) < 0 | z==0);
        errornow = sqrt(sum((z' - ([z_P1;z_P2]')*alpha).^2));
    end
    if count == 1 && errorold < errornow
            z_P1 = z_P1tmp;
            z_P2 = z_P2tmp;
            K_1over2 = K_1over2tmp;
    end
    count = count + 1;
end

end
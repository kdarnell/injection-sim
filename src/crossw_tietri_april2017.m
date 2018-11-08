function [dot_prod] = crossw_tietri_april2017(z,alpha,alpha_func)

% Retrieve interpolated tie-triangles using a parameter that represents the
% normalized distance "alpha" from face to face of quaternary diagram.
% At present (04/04/17), alpha_func returns length(alpha) x 4 x 3 matrix
% where the third dimension represents Aq, V, H compositions.
% I'm assuming right now that alpha is either size(1,1), a scalar, or,
% alpha is size(length(z(:,1)),1);
% alpha might also be a totally different shape where the output should be
% length(z(:,1)) x length(alpha(1,:))
alpha_sz = length(alpha);
z_sz = size(z);

pts = alpha_func(alpha(:));
cols = 1:3; %We are neglecting nitrogen (perhaps I should be converting this to xyz space)...
A_vec = pts(:,cols,2)-pts(:,cols,1);
B_vec = pts(:,cols,3)-pts(:,cols,1);

% Returns a matrix length(alpha(:)) x 3
alpha_planes = cross(A_vec, B_vec, 2);

if (z_sz(1) == alpha_sz(1)) || (alpha_sz(1) == 1)
    % Returns a matrix length(z(:,1)) x 3
    z_vec = bsxfun(@minus,z(:,cols),pts(:,cols,1));
    dot_prod = sum(bsxfun(@times,z_vec,alpha_planes), 2);
    dot_prod((alpha > 1)) = dot_prod((alpha > 1)) + 1e3*(1 - alpha(alpha > 1)).^2;
    dot_prod((alpha < 0 )) = dot_prod((alpha < 0)) + 1e3*(alpha(alpha < 0)).^2;

else
    % Returns a matrix length(z(:,1)) x 3 x length(alpha)
    z_vec = bsxfun(@minus,z(:,cols),permute(pts(:,cols,1),[3 2 1]));
    dot_prod = squeeze(sum(bsxfun(@times,z_vec,permute(alpha_planes,[3 2 1])),2));
    dot_prod(:,(alpha > 1)) = dot_prod(:,(alpha > 1)) + 1e3*(1 - alpha(alpha > 1)).^2;
    dot_prod(:,(alpha < 0 )) = dot_prod(:,(alpha < 0)) + 1e3*(alpha(alpha > 1)).^2;
end

% dot_prod will be length(z(:,1)) x length(alpha) or length(z(:,1)) x if z
% and alpha are of the same length.

end
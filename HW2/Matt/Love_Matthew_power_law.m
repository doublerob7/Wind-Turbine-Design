% power_law provides the power law estimation for u:
% u/u_ref = (z/z_ref)^alpha
    % wherein
        % u = estimated windspeed at desired height
        % u_ref = windspeed at which data was taken
        % z = height at which windspeed will be estimated
        % z_ref = height at which data was taken
        % alpha is the roughness coefficient (the power in the power law)
function [u] = Love_Matthew_power_law(u_ref,z_ref,z,alpha)

% Power Law:
u = u_ref*(z/z_ref)^alpha;
function results = compute_control_min_energy(chi_hist,c)

    % Computes control history given time history of state+co-state
    
%     alpha_hist = NaN(3,size(chi_hist,2));
%     u_hist = NaN(1,size(chi_hist,2));
%     S_hist = NaN(1,size(chi_hist,2));
    
    alpha_hist = [];
    u_hist = [];
    S_hist = [];
    
    for i = 1:size(chi_hist,2)
        lambdaV = chi_hist(11:13,i);
        if norm(lambdaV)==0
            lambdaVhat = zeros(3,1);
        else
            lambdaVhat = lambdaV/norm(lambdaV);
        end
        alpha = -lambdaVhat;
        S = switch_func_eps1(norm(lambdaV), chi_hist(14,i), chi_hist(7,i), c);
        if S < -1
            u = 1;
        elseif abs(S) <= 1
            u = -S/2;
        elseif S > 1
            u = 0;
        end
        % Save
        alpha_hist = [alpha_hist,alpha];
        u_hist = [u_hist, u];
        S_hist = [S_hist, S];
    end
     
     %% Save results
     results = struct('u_hist',u_hist,'S_hist',S_hist,'alpha_hist',alpha_hist);
end
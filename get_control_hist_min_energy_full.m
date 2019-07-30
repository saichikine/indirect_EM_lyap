function results = get_control_hist_min_energy_full(chi_hist,c)

    % Computes control history given time history of state+co-state
    
    length_results = size(chi_hist,1);
    
    alpha_hist = NaN(3,length_results);
    u_hist = NaN(1,length_results);
    u_vec_hist = NaN(3,length_results);
    S_hist = NaN(1,length_results);
    
    for i = 1:length_results
        m = chi_hist(i,7);
        lambdam = chi_hist(i,14);
        lambdaV = chi_hist(i,11:13);
        lambdaV_mag = norm(lambdaV);
        if lambdaV_mag==0
            lambdaVhat = zeros(3,1);
        else
            lambdaVhat = lambdaV/norm(lambdaV);
        end
        alpha = -lambdaVhat;
        S = -lambdaV_mag*c/m - lambdam + 1;
        if S < -1
        u = 1;
    %     elseif abs(S) <= 1
    %         u = -S/2;
        elseif S >= -1 && S <= 1
            u = (1-S)/2;
        elseif S > 1
            u = 0;
        else
            error("S is not real, S=%d\n",S)
        end
        % Save
        alpha_hist(:,i) = alpha;
        u_hist(i) = u;
        u_vec_hist(:,i) = u.*alpha;
        S_hist(i) = S;
    end
     
     %% Save results
     results = struct('alpha_hist',alpha_hist,'u_hist',u_hist,'u_vec_hist',u_vec_hist,'S_hist',S_hist);
end
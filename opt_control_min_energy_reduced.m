function u = opt_control_min_energy_reduced(chi,m,params)

    % chi is combined state and costate vector
    % c is exhaust velocity
    
    c = params.c/params.vel_norm;
    lambdaV = chi(10:12); % velocity costate
    lambdaV_mag = norm(lambdaV);

    S = -lambdaV_mag*c/(m/params.m_norm) + 1; % switching function for min energy case
    %fprintf("S = %d\n",S);

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

end
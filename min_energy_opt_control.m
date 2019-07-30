function u = min_energy_opt_control(chi,c)

    % chi is combined state and costate vector
    % c is exhaust velocity

    m = chi(7); % mass
    lambdaV = chi(11:13); % velocity costate
    lambdaV_mag = norm(lambdaV);
    lambdam = chi(14); % mass costate

    S = -lambdaV_mag*c/m - lambdam; % switching function for min energy case

    if S < -1
        u = 1;
    elseif abs(S) <= 1
        u = -S/2;
    elseif S > 1
        u = 0;
    end

end
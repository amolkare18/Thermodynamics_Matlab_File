function mainfunction()
    % Define constants
    tc = 190.564;          % Critical temperature (K)
    pc = 45.992 * 10^5;    % Critical pressure (Pa)
    w = 0.01142;           % Acentric factor
    R = 8.314;             % Gas constant (Pa*m^3/(mol*K))
    a = 0.2283;            % Value of a
    b = 4.28 * 10^(-5);     % Value of b 
    

     % Constants for ideal gas specific heat capacity as a function of T
    a1 = 19.875;
    b1 = 5.021e-2;
    c1 = 1.268e-5;
    d1 = -11.004e-9;

    % Calculate k using the acentric factor w`
    k = 0.37464 + 1.54226 * w - 0.26992 * w^2;

    % Define a range of temperatures and pressures
    T_range = 150:10:350;      % Temperature range in Kelvin
    p_range = 1e5:1e5:10e5;   % Pressure range in Pa 

    % Initialize arrays to store results
    enthalpy_values = zeros(length(T_range), length(p_range));
    entropy_values = zeros(length(T_range), length(p_range));
    molar_volume_values = zeros(length(T_range), length(p_range));

    % Loop through each temperature and pressure, and calculate properties
    for i = 1:length(T_range)
        T = T_range(i);
        alpha = alphat(k, T, tc);

        for j = 1:length(p_range)
            p = p_range(j);
            Z = calculate_Z(p, T, a, b, R);
            enthalpy_values(i, j) = calculate_enthalpy(T, p, b, Z);
            entropy_values(i, j) = calculate_entropy(T, p, Z, b);
            molar_volume_values(i, j) = Z * R * T / p;  % Molar volume calculation

            % Display the results for each temperature and pressure
            fprintf('Temperature: %.2f K, Pressure: %.2f bar, Enthalpy: %.4f J/mol, Entropy: %.4f J/(mol*K), Molar Volume: %.7f m^3/mol\n', T, p / 1e5, enthalpy_values(i, j), entropy_values(i, j), molar_volume_values(i, j));
        end
    end

    
    % Plot p vs enthalpy, t vs entropy, and p vs vbar
    figure;

    % P vs Enthalpy
    subplot(3, 1, 1);
    plot( enthalpy_values/1e3,p_range / 1e5);  % Convert Pa to bar for x-axis
    xlabel('Enthalpy (J/mol)');
    ylabel('Pressure (bar)');
    title('P vs Enthalpy');
    legend(string(T_range), 'Location', 'best');

    % T vs Entropy
    subplot(3, 1, 2);
    plot( entropy_values,T_range); 
    xlabel('Entropy (J/(mol*K))');
    ylabel('Temperature (K)');
    title('T vs Entropy');
    legend(string(p_range / 1e5), 'Location', 'best');  % Use bar for pressure in legend

    % P vs Molar Volume
    subplot(3, 1, 3);
    plot( molar_volume_values,p_range / 1e5);
    xlabel('Molar Volume (m^3/mol)');
    ylabel('Pressure (bar)');
    title('P vs Molar Volume');
    legend(string(T_range), 'Location', 'best');

    hold off;

    % Nested functions
    function alpha = alphat(k, T, tc)
        alpha = (1 + k * (1 - sqrt(T / tc)))^2;
    end

    function Z = calculate_Z(p, T, a, b, R)
        B = b * p / (R * T);
        A = a * p / (R^2 * T^2);

        % Coefficients of the cubic equation for Z
        coeffs = [1, (-1 + B), (A - 3 * B^2 - 2 * B), -(A * B - B^2 - B^3)];

        % Solve the cubic equation numerically
        z_roots = roots(coeffs);

        % Filter out real, positive roots
        real_positive_roots = z_roots(imag(z_roots) == 0 & real(z_roots) > 0);

        
        if ~isempty(real_positive_roots)
            Z = real_positive_roots(1); 
        else
            Z = NaN;  % If no real positive solutions are found
            warning('No valid Z found for T = %.2f K and P = %.2f Pa.', T, p);
        end
    end

    function d_alpha = dalpha_dT(k, T, tc)
        d_alpha = (-k/sqrt(T*tc))*(k*(1-sqrt(T/tc))+1);
        %(-k / (2 * sqrt(T * tc))) * (1 + k * (1 - sqrt(T / tc)));
    end

    function integral_cp = cp_dt(T_start, T_end)
        integral_cp = a1 * (T_end - T_start) + (b1 / 2) * (T_end^2 - T_start^2) + (c1 / 3) * (T_end^3 - T_start^3) + (d1 / 4) * (T_end^4 - T_start^4);
                      
    end

    
       
        
      function enthalpy = calculate_enthalpy(T, p, b, Z)
        B = b * p / (R * T);

        % Reference temperature and pressure
        T_ref = 298.15;
        p_ref=10^5;
          
        % Integral of cp over the temperature range
        integral_cp = cp_dt(T_ref, T);
        
        % Calculate dAlpha/dT and delta_h
        dAlpha_dT = dalpha_dT(k, T, tc);
        delta_h = R * T * (Z - 1) + ((T * dAlpha_dT * a - a) / (2.828 * b)) * log((Z + (1 + sqrt(2)) * B) / (Z + (1 - sqrt(2)) * B));
     
        enthalpy = integral_cp + delta_h;
    end

    function delta_s_ideal_gas = integral_cp_by_t(T_start, T_end, p)
        delta_s_ideal_gas = a1 * log(T_end / T_start) + b1 * (T_end - T_start) + c1 * (T_end^2 - T_start^2) / 2 + d1 * (T_end^3 - T_start^3) / 3 - R * log(p / 1e5);
                            
    end

    function entropy = calculate_entropy(T, p, Z, b)
        B = b * p / (R * T);

        % Reference temperature and pressure
        T_ref = 298.15;
         p_ref=10^5;
       % Integral of cp/T over the temperature range
        delta_s_ideal_gas = integral_cp_by_t(T_ref, T, p);
       % Calculate dAlpha/dT and delta_s
        dAlpha_dT = dalpha_dT(k, T, tc);
        delta_s = R * (Z - B) + (dAlpha_dT * a) / (2.828 * b) * log((Z + (1 + sqrt(2)) * B) / (Z + (1 - sqrt(2)) * B));
        
        entropy = delta_s_ideal_gas + delta_s;
    end
end

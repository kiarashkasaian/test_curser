function focused_beam_stratified_medium()
    % FOCUSED_BEAM_STRATIFIED_MEDIUM
    % Calculates electric field of focused beam in stratified medium using Richards-Wolf integral
    % Includes three layers (oil, glass, water), aberrations, and Gaussian beam profile
    
    clear; close all; clc;
    
    %% Parameters
    % Objective parameters
    NA = 1.4;                    % Numerical aperture
    f = 3e-3;                    % Focal length (3mm)
    
    % Wavelength
    lambda = 532e-9;             % Wavelength in vacuum (532 nm)
    k0 = 2*pi/lambda;            % Wave number in vacuum
    
    % Gaussian beam parameter
    sigma = 0.8;                 % Sigma in units of NA
    
    % Refractive indices (design values)
    n_oil_design = 1.518;        % Oil design RI
    n_glass_design = 1.518;      % Glass design RI  
    n_water_design = 1.333;      % Water design RI
    
    % Actual refractive indices (with aberrations)
    n_oil_actual = 1.515;        % Actual oil RI
    n_glass_actual = 1.520;      % Actual glass RI
    n_water_actual = 1.333;      % Actual water RI
    
    % Layer thicknesses
    t_oil = 100e-6;              % Oil thickness (100 μm)
    t_glass = 170e-6;            % Glass thickness (170 μm)
    
    % Coordinate system
    x_range = linspace(-2*lambda, 2*lambda, 201);
    z_range = linspace(-4*lambda, 4*lambda, 401);
    [X, Z] = meshgrid(x_range, z_range);
    
    % Initialize field components
    Ex = zeros(size(X));
    Ey = zeros(size(X));
    Ez = zeros(size(X));
    
    %% Richards-Wolf Integration
    fprintf('Calculating electric field using Richards-Wolf integral...\n');
    
    % Integration parameters
    theta_max = asin(NA/n_water_actual);  % Maximum collection angle
    N_theta = 200;                        % Number of theta points
    N_phi = 400;                          % Number of phi points
    
    theta = linspace(0, theta_max, N_theta);
    phi = linspace(0, 2*pi, N_phi);
    
    % Remove theta=0 to avoid singularity
    theta = theta(2:end);
    
    [THETA, PHI] = meshgrid(theta, phi);
    
    % Calculate field for each point in x-z plane (y=0)
    for ix = 1:length(x_range)
        if mod(ix, 20) == 0
            fprintf('Progress: %d/%d\n', ix, length(x_range));
        end
        
        for iz = 1:length(z_range)
            x_pos = x_range(ix);
            z_pos = z_range(iz);
            
            % Calculate field components at this point
            [ex, ey, ez] = calculate_field_point(x_pos, 0, z_pos, THETA, PHI, ...
                n_oil_design, n_glass_design, n_water_design, ...
                n_oil_actual, n_glass_actual, n_water_actual, ...
                t_oil, t_glass, k0, f, sigma, NA);
            
            Ex(iz, ix) = ex;
            Ey(iz, ix) = ey;
            Ez(iz, ix) = ez;
        end
    end
    
    %% Plotting
    fprintf('Generating plots...\n');
    
    % Convert to micrometers for plotting
    x_um = x_range * 1e6;
    z_um = z_range * 1e6;
    
    % Create figure with subplots
    figure('Position', [100, 100, 1400, 1000]);
    
    % Ex component - amplitude
    subplot(3, 2, 1);
    imagesc(x_um, z_um, abs(Ex));
    axis image; colorbar;
    title('|E_x| Amplitude');
    xlabel('x (μm)'); ylabel('z (μm)');
    
    % Ex component - phase
    subplot(3, 2, 2);
    imagesc(x_um, z_um, angle(Ex));
    axis image; colorbar;
    title('E_x Phase');
    xlabel('x (μm)'); ylabel('z (μm)');
    
    % Ey component - amplitude
    subplot(3, 2, 3);
    imagesc(x_um, z_um, abs(Ey));
    axis image; colorbar;
    title('|E_y| Amplitude');
    xlabel('x (μm)'); ylabel('z (μm)');
    
    % Ey component - phase
    subplot(3, 2, 4);
    imagesc(x_um, z_um, angle(Ey));
    axis image; colorbar;
    title('E_y Phase');
    xlabel('x (μm)'); ylabel('z (μm)');
    
    % Ez component - amplitude
    subplot(3, 2, 5);
    imagesc(x_um, z_um, abs(Ez));
    axis image; colorbar;
    title('|E_z| Amplitude');
    xlabel('x (μm)'); ylabel('z (μm)');
    
    % Ez component - phase
    subplot(3, 2, 6);
    imagesc(x_um, z_um, angle(Ez));
    axis image; colorbar;
    title('E_z Phase');
    xlabel('x (μm)'); ylabel('z (μm)');
    
    sgtitle('Electric Field Components in x-z Cross Section');
    
    % Total intensity plot
    figure('Position', [200, 200, 800, 600]);
    I_total = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
    imagesc(x_um, z_um, I_total);
    axis image; colorbar;
    title('Total Intensity |E|^2');
    xlabel('x (μm)'); ylabel('z (μm)');
    
    fprintf('Calculation complete!\n');
end

function [ex, ey, ez] = calculate_field_point(x, y, z, THETA, PHI, ...
    n_oil_design, n_glass_design, n_water_design, ...
    n_oil_actual, n_glass_actual, n_water_actual, ...
    t_oil, t_glass, k0, f, sigma, NA)
    
    % Calculate electric field at a single point using Richards-Wolf integral
    
    % Gaussian apodization function
    sin_theta = sin(THETA);
    apodization = exp(-0.5 * (sin_theta / (sigma * NA / n_water_actual)).^2);
    
    % Transmission coefficients and phase aberrations
    [t_coeff, phase_aberr] = calculate_transmission_and_aberration(...
        THETA, n_oil_design, n_glass_design, n_water_design, ...
        n_oil_actual, n_glass_actual, n_water_actual, ...
        t_oil, t_glass, k0);
    
    % Polarization factors for x-polarized input
    cos_phi = cos(PHI);
    sin_phi = sin(PHI);
    cos_theta = cos(THETA);
    
    % Amplitude factors
    A0 = sqrt(cos_theta) .* apodization .* t_coeff .* exp(1i * phase_aberr);
    
    % Polarization components in spherical coordinates
    A_s = A0 .* cos_phi;  % s-polarization (tangential to theta direction)
    A_p = A0 .* sin_phi .* cos_theta;  % p-polarization (in theta-phi plane)
    
    % Phase factors
    k_water = k0 * n_water_actual;
    phase = k_water * (x * sin_theta .* cos_phi + y * sin_theta .* sin_phi + z * cos_theta);
    
    % Integration weights
    dtheta = THETA(2,1) - THETA(1,1);
    dphi = PHI(1,2) - PHI(1,1);
    sin_theta_weight = sin_theta * dtheta * dphi;
    
    % Integrate to get field components
    integrand_x = A_s .* cos_phi .* cos_theta - A_p .* sin_phi;
    integrand_y = A_s .* sin_phi .* cos_theta + A_p .* cos_phi;
    integrand_z = -A_s .* sin_theta;
    
    ex = sum(sum(integrand_x .* exp(1i * phase) .* sin_theta_weight));
    ey = sum(sum(integrand_y .* exp(1i * phase) .* sin_theta_weight));
    ez = sum(sum(integrand_z .* exp(1i * phase) .* sin_theta_weight));
    
    % Normalization factor
    norm_factor = -1i * k_water * f / (2*pi);
    ex = ex * norm_factor;
    ey = ey * norm_factor;
    ez = ez * norm_factor;
end

function [t_coeff, phase_aberr] = calculate_transmission_and_aberration(...
    THETA, n_oil_design, n_glass_design, n_water_design, ...
    n_oil_actual, n_glass_actual, n_water_actual, ...
    t_oil, t_glass, k0)
    
    % Calculate transmission coefficients and phase aberrations
    
    sin_theta = sin(THETA);
    cos_theta = cos(THETA);
    
    % Snell's law for actual media
    sin_theta_oil = n_water_actual * sin_theta / n_oil_actual;
    sin_theta_glass = n_water_actual * sin_theta / n_glass_actual;
    
    cos_theta_oil = sqrt(1 - sin_theta_oil.^2);
    cos_theta_glass = sqrt(1 - sin_theta_glass.^2);
    
    % Fresnel transmission coefficients (assuming s-polarization dominance)
    % Water-oil interface
    t_wo = 2 * n_water_actual * cos_theta ./ ...
           (n_water_actual * cos_theta + n_oil_actual * cos_theta_oil);
    
    % Oil-glass interface  
    t_og = 2 * n_oil_actual * cos_theta_oil ./ ...
           (n_oil_actual * cos_theta_oil + n_glass_actual * cos_theta_glass);
    
    % Total transmission coefficient
    t_coeff = t_wo .* t_og;
    
    % Phase aberrations due to RI mismatch
    % Design path length
    path_design = t_oil * n_oil_design + t_glass * n_glass_design;
    
    % Actual path length  
    path_actual = t_oil * n_oil_actual ./ cos_theta_oil + ...
                  t_glass * n_glass_actual ./ cos_theta_glass;
    
    % Phase aberration
    phase_aberr = k0 * (path_actual - path_design);
    
    % Handle total internal reflection (set transmission to zero)
    tir_mask = (sin_theta_oil > 1) | (sin_theta_glass > 1);
    t_coeff(tir_mask) = 0;
    phase_aberr(tir_mask) = 0;
end
function focused_beam_stratified_medium()
    % Calculate electric field of focused beam in stratified medium
    % using Richards-Wolf integral with aberrations
    
    close all; clear; clc;
    
    %% Parameters
    % Optical parameters
    NA = 1.4;                    % Numerical aperture
    lambda = 0.532e-6;           % Wavelength in vacuum (532 nm)
    sigma = 0.7;                 % Gaussian beam sigma in units of NA
    
    % Refractive indices (design values)
    n_oil_design = 1.518;        % Oil design RI
    n_glass_design = 1.518;      % Glass design RI  
    n_water = 1.33;              % Water RI
    
    % Actual refractive indices (with aberrations)
    n_oil_actual = 1.515;        % Actual oil RI
    n_glass_actual = 1.520;      % Actual glass RI
    
    % Layer thicknesses
    t_oil = 100e-6;              % Oil layer thickness (100 μm)
    t_glass = 170e-6;            % Glass coverslip thickness (170 μm)
    
    % Calculation grid
    x_range = 2e-6;              % ±2 μm in x
    z_range = 4e-6;              % ±4 μm in z (around focus)
    Nx = 201;                    % Number of x points
    Nz = 201;                    % Number of z points
    
    % Angular integration parameters
    N_theta = 200;               % Number of theta points
    N_phi = 100;                 % Number of phi points
    
    %% Setup coordinate system
    x = linspace(-x_range, x_range, Nx);
    z = linspace(-z_range, z_range, Nz);
    [X, Z] = meshgrid(x, z);
    
    % Angular coordinates
    theta_max = asin(NA/n_water);  % Maximum collection angle in water
    theta = linspace(0, theta_max, N_theta);
    phi = linspace(0, 2*pi, N_phi);
    [THETA, PHI] = meshgrid(theta, phi);
    
    %% Initialize field arrays
    Ex = zeros(Nz, Nx);
    Ey = zeros(Nz, Nx);
    Ez = zeros(Nz, Nx);
    
    %% Calculate wave vectors
    k0 = 2*pi/lambda;            % Wave vector in vacuum
    k_water = k0 * n_water;      % Wave vector in water
    
    %% Richards-Wolf integration
    fprintf('Calculating electric field using Richards-Wolf integral...\n');
    
    for iz = 1:Nz
        if mod(iz, 20) == 0
            fprintf('Progress: %d/%d\n', iz, Nz);
        end
        
        for ix = 1:Nx
            x_pos = x(ix);
            z_pos = z(iz);
            
            % Initialize field components at this point
            Ex_point = 0;
            Ey_point = 0;
            Ez_point = 0;
            
            % Integrate over angular spectrum
            for itheta = 1:N_theta
                for iphi = 1:N_phi
                    theta_val = theta(itheta);
                    phi_val = phi(iphi);
                    
                    % Skip if beyond NA
                    if sin(theta_val) > NA/n_water
                        continue;
                    end
                    
                    % Calculate aberration phase
                    aberration_phase = calculate_aberration_phase(theta_val, ...
                        n_oil_design, n_oil_actual, t_oil, ...
                        n_glass_design, n_glass_actual, t_glass, ...
                        n_water, k0);
                    
                    % Gaussian apodization at back focal plane
                    NA_local = n_water * sin(theta_val);
                    gauss_factor = exp(-(NA_local/(sigma*NA))^2);
                    
                    % Fresnel transmission coefficients
                    [t_s, t_p] = fresnel_coefficients(theta_val, ...
                        n_oil_actual, n_glass_actual, n_water);
                    
                    % Direction cosines in water
                    l = sin(theta_val) * cos(phi_val);
                    m = sin(theta_val) * sin(phi_val);
                    n = cos(theta_val);
                    
                    % Amplitude factors
                    A = sqrt(cos(theta_val)) * sin(theta_val) * gauss_factor;
                    
                    % Phase factor
                    phase = k_water * (l*x_pos + n*z_pos) + aberration_phase;
                    exp_phase = exp(1i * phase);
                    
                    % Polarization components (input is x-polarized)
                    % Transform to local s,p coordinates and back
                    cos_phi = cos(phi_val);
                    sin_phi = sin(phi_val);
                    
                    % s-polarization unit vector
                    s_x = -sin_phi;
                    s_y = cos_phi;
                    s_z = 0;
                    
                    % p-polarization unit vector  
                    p_x = cos(theta_val) * cos_phi;
                    p_y = cos(theta_val) * sin_phi;
                    p_z = -sin(theta_val);
                    
                    % Input field components in s,p basis
                    E_s = cos_phi;  % x-pol projected to s
                    E_p = -sin_phi * sin(theta_val);  % x-pol projected to p
                    
                    % Apply transmission coefficients
                    E_s_trans = E_s * t_s;
                    E_p_trans = E_p * t_p;
                    
                    % Transform back to Cartesian coordinates
                    Ex_contrib = A * (E_s_trans * s_x + E_p_trans * p_x) * exp_phase;
                    Ey_contrib = A * (E_s_trans * s_y + E_p_trans * p_y) * exp_phase;
                    Ez_contrib = A * (E_s_trans * s_z + E_p_trans * p_z) * exp_phase;
                    
                    % Add contributions
                    Ex_point = Ex_point + Ex_contrib;
                    Ey_point = Ey_point + Ey_contrib;
                    Ez_point = Ez_point + Ez_contrib;
                end
            end
            
            % Store results
            Ex(iz, ix) = Ex_point;
            Ey(iz, ix) = Ey_point;
            Ez(iz, ix) = Ez_point;
        end
    end
    
    % Normalize by integration area
    dtheta = theta_max / N_theta;
    dphi = 2*pi / N_phi;
    Ex = Ex * dtheta * dphi;
    Ey = Ey * dtheta * dphi;
    Ez = Ez * dtheta * dphi;
    
    %% Plot results
    plot_results(x*1e6, z*1e6, Ex, Ey, Ez);
    
    fprintf('Calculation complete!\n');
end

function aberration_phase = calculate_aberration_phase(theta, ...
    n_oil_design, n_oil_actual, t_oil, ...
    n_glass_design, n_glass_actual, t_glass, ...
    n_water, k0)
    % Calculate aberration phase due to RI mismatch
    
    % Snell's law to find angles in each medium
    sin_theta_water = sin(theta);
    sin_theta_glass_design = n_water * sin_theta_water / n_glass_design;
    sin_theta_glass_actual = n_water * sin_theta_water / n_glass_actual;
    sin_theta_oil_design = n_water * sin_theta_water / n_oil_design;
    sin_theta_oil_actual = n_water * sin_theta_water / n_oil_actual;
    
    % Check for total internal reflection
    if sin_theta_glass_actual > 1 || sin_theta_oil_actual > 1
        aberration_phase = 0;
        return;
    end
    
    theta_glass_design = asin(sin_theta_glass_design);
    theta_glass_actual = asin(sin_theta_glass_actual);
    theta_oil_design = asin(sin_theta_oil_design);
    theta_oil_actual = asin(sin_theta_oil_actual);
    
    % Optical path difference
    opd_glass = t_glass * (n_glass_actual * cos(theta_glass_actual) - ...
                          n_glass_design * cos(theta_glass_design));
    opd_oil = t_oil * (n_oil_actual * cos(theta_oil_actual) - ...
                      n_oil_design * cos(theta_oil_design));
    
    % Total aberration phase
    aberration_phase = k0 * (opd_glass + opd_oil);
end

function [t_s, t_p] = fresnel_coefficients(theta_water, n_oil, n_glass, n_water)
    % Calculate Fresnel transmission coefficients through the layers
    
    % Angles in each medium (Snell's law)
    sin_theta_water = sin(theta_water);
    sin_theta_glass = n_water * sin_theta_water / n_glass;
    sin_theta_oil = n_water * sin_theta_water / n_oil;
    
    % Check for total internal reflection
    if sin_theta_glass > 1 || sin_theta_oil > 1
        t_s = 0;
        t_p = 0;
        return;
    end
    
    theta_glass = asin(sin_theta_glass);
    theta_oil = asin(sin_theta_oil);
    
    cos_theta_water = cos(theta_water);
    cos_theta_glass = cos(theta_glass);
    cos_theta_oil = cos(theta_oil);
    
    % Fresnel coefficients for each interface
    % Water-glass interface
    t_s_wg = 2 * n_water * cos_theta_water / ...
             (n_water * cos_theta_water + n_glass * cos_theta_glass);
    t_p_wg = 2 * n_water * cos_theta_water / ...
             (n_glass * cos_theta_water + n_water * cos_theta_glass);
    
    % Glass-oil interface  
    t_s_go = 2 * n_glass * cos_theta_glass / ...
             (n_glass * cos_theta_glass + n_oil * cos_theta_oil);
    t_p_go = 2 * n_glass * cos_theta_glass / ...
             (n_oil * cos_theta_glass + n_glass * cos_theta_oil);
    
    % Total transmission (product of individual transmissions)
    t_s = t_s_wg * t_s_go;
    t_p = t_p_wg * t_p_go;
end

function plot_results(x_um, z_um, Ex, Ey, Ez)
    % Plot the electric field components
    
    figure('Position', [100, 100, 1400, 1000]);
    
    % Field components to plot
    fields = {Ex, Ey, Ez};
    field_names = {'E_x', 'E_y', 'E_z'};
    
    for i = 1:3
        field = fields{i};
        
        % Amplitude plot
        subplot(3, 2, 2*i-1);
        imagesc(x_um, z_um, abs(field));
        axis image;
        colorbar;
        xlabel('x (μm)');
        ylabel('z (μm)');
        title(sprintf('|%s| (Amplitude)', field_names{i}));
        colormap(gca, 'hot');
        
        % Phase plot
        subplot(3, 2, 2*i);
        imagesc(x_um, z_um, angle(field));
        axis image;
        colorbar;
        xlabel('x (μm)');
        ylabel('z (μm)');
        title(sprintf('arg(%s) (Phase)', field_names{i}));
        colormap(gca, 'hsv');
        caxis([-pi, pi]);
    end
    
    sgtitle('Electric Field Components in x-z Cross Section');
    
    % Additional plot: Total intensity
    figure('Position', [200, 200, 800, 600]);
    I_total = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
    imagesc(x_um, z_um, I_total);
    axis image;
    colorbar;
    xlabel('x (μm)');
    ylabel('z (μm)');
    title('Total Intensity |E|^2');
    colormap('hot');
    
    % Line plots through focus
    figure('Position', [300, 300, 1200, 400]);
    
    % Find center indices
    center_x = round(length(x_um)/2);
    center_z = round(length(z_um)/2);
    
    % Plot along z-axis (x=0)
    subplot(1, 2, 1);
    plot(z_um, abs(Ex(center_x, :)), 'r-', 'LineWidth', 2); hold on;
    plot(z_um, abs(Ey(center_x, :)), 'g-', 'LineWidth', 2);
    plot(z_um, abs(Ez(center_x, :)), 'b-', 'LineWidth', 2);
    xlabel('z (μm)');
    ylabel('|E| (a.u.)');
    title('Field Amplitude along z-axis (x=0)');
    legend('|E_x|', '|E_y|', '|E_z|');
    grid on;
    
    % Plot along x-axis (z=0)
    subplot(1, 2, 2);
    plot(x_um, abs(Ex(:, center_z)), 'r-', 'LineWidth', 2); hold on;
    plot(x_um, abs(Ey(:, center_z)), 'g-', 'LineWidth', 2);
    plot(x_um, abs(Ez(:, center_z)), 'b-', 'LineWidth', 2);
    xlabel('x (μm)');
    ylabel('|E| (a.u.)');
    title('Field Amplitude along x-axis (z=0)');
    legend('|E_x|', '|E_y|', '|E_z|');
    grid on;
end

% Run the calculation
focused_beam_stratified_medium();
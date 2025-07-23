% FOCUSED_BEAM_STRATIFIED.M
% ------------------------------------------------------------
% Calculates the focal electric field distribution of a high-NA
% objective when the focus lies in water but the immersion/oil and
% cover-glass indices deviate from their design values.  The Richards–Wolf
% vectorial diffraction integral is evaluated in a stratified medium
% (oil → glass → water).  A Gaussian pupil illumination with user-selectable
% width (given in NA units) and linear x-polarisation is assumed.
% The script visualises x–z cross-sections through the focus showing
% amplitude and phase of the E-field components.
% ------------------------------------------------------------

clear; clc;

%% 1. INPUT PARAMETERS ---------------------------------------------------

% Wavelength in vacuum [m]
lambda0      = 532e-9;          % 532 nm (green)
% Numerical aperture of the objective
NA           = 1.40;
% Focal length of the objective [m] (needed only for scaling constants).
% The exact value is irrelevant for the normalised focal distribution, but
% choose something realistic (e.g. 3 mm) so that the physical prefactor has
% the correct magnitude if desired.
f_obj        = 3e-3;            % 3 mm objective focal length

% Refractive indices (actual)
n_oil        = 1.515;           % immersion oil (actual)
n_glass      = 1.520;           % cover glass (actual)
n_water      = 1.333;           % sample medium (actual)  – focus is here

% Refractive indices that the objective was designed for
n_oil_des    = 1.518;           % design value oil
n_glass_des  = 1.525;           % design value glass

% Layer thicknesses [m] (along the optical axis)
d_oil        = 100e-6;          % axial thickness of the immersion oil gap
% For a high-NA objective the nominal cover-glass thickness is 170 μm.
d_glass      = 170e-6;          % thickness of the glass coverslip

% Width (1/e) of the Gaussian amplitude at the pupil, given in units of NA
sigmaNA      = 0.50;            % 0.5 × NA (≈ moderate over-filling)

% Sampling of the x–z plane that will be plotted
x_range      = 1.5e-6;          % half-width in x [m] around focus
z_range      = 3.0e-6;          % half-width in z [m]
Nx           = 301;             % #points in x
Nz           = 401;             % #points in z

% Angular sampling of the Richards–Wolf integral
Ntheta       = 250;             % polar angle samples (0…θ_max)
Nphi         = 360;             % azimuth samples (0…2π)

%% 2. PRE-COMPUTE CONSTANTS ---------------------------------------------

k0     = 2*pi/lambda0;          % vacuum wavenumber
k_w    = k0 * n_water;          % wavenumber in water (image space)

% Aperture half-angle inside the immersion oil
theta_max_oil = asin(NA / n_oil);

% Discretise pupil angles in oil space
theta = linspace(0, theta_max_oil, Ntheta);             % 1×Ntheta
phi   = linspace(0, 2*pi,        Nphi  );               % 1×Nphi
[THETA, PHI] = meshgrid(theta, phi);                    % size Nphi×Ntheta

% Convert angles into the water layer using Snell’s law
sin_theta_oil   = sin(THETA);
% Oil → glass
sin_theta_glass = (n_oil / n_glass) * sin_theta_oil;
% Glass → water
sin_theta_water = (n_oil / n_water) * sin_theta_oil;

% Stop rays hitting total internal reflection at either interface
mask_valid = (sin_theta_glass<=1) & (sin_theta_water<=1);

% Replace invalid rays by NaNs so that they contribute zero
sin_theta_glass(~mask_valid)  = NaN;
sin_theta_water(~mask_valid)  = NaN;
THETA(~mask_valid)            = NaN;  % needed later for cos()

cos_theta_oil   = sqrt(1 - sin_theta_oil .^2);
cos_theta_glass = sqrt(1 - sin_theta_glass.^2);
cos_theta_water = sqrt(1 - sin_theta_water.^2);

%% 3. GAUSSIAN AMPLITUDE & PHASE ABERRATION ------------------------------

% Normalised pupil radius ρ = sinθ / sinθ_max
theta_max_norm = sin(theta_max_oil);
rho = sin_theta_oil ./ theta_max_norm;                 % Nphi×Ntheta

% Gaussian illumination amplitude in the pupil
A_pupil = exp( -(rho.^2) / (2*sigmaNA^2) );

% Phase aberration from refractive-index mismatch (Gibson & Lanni-style)
% Optical path inside each layer relative to design condition
OPD = k0 * ( ...
    n_oil   * d_oil   ./ cos_theta_oil   + ...
    n_glass * d_glass ./ cos_theta_glass - ...
    n_oil_des   * d_oil   ./ cos_theta_oil   - ...
    n_glass_des * d_glass ./ cos_theta_glass );

phase_aberr = exp(1i * OPD);

% Total pupil function (scalar for amplitude + phase)
Pupil = A_pupil .* phase_aberr .* mask_valid;          % invalid rays → 0

%% 4. POLARISATION VECTORS AT THE PUPIL ----------------------------------
% Incident linear x-polarisation expressed in cylindrical basis at pupil
% (Richards & Wolf 1959, Eqs. (59))
ex_p = cos(PHI);       % projections of x-polarised vector
ey_p = sin(PHI);
ez_p = zeros(size(PHI));

%% 5. PREPARE OBJECT SPACE GRID ------------------------------------------

x = linspace(-x_range, x_range, Nx);          % 1×Nx
z = linspace(-z_range, z_range, Nz);          % 1×Nz
[X,Z] = meshgrid(x, z);                       % Nz×Nx
Y = zeros(size(X));                           % x–z plane ⇒ y = 0

%% 6. EVALUATE THE RICHARDS–WOLF INTEGRAL --------------------------------

fprintf('Integrating %d×%d angular samples onto %d×%d spatial grid…\n', ...
    Nphi, Ntheta, Nz, Nx);

% Differential elements
dtheta = theta(2) - theta(1);
dphi   = phi(2)   - phi(1);

% Pre-allocate complex fields
Ex = complex(zeros(size(X)));
Ey = complex(zeros(size(X)));
Ez = complex(zeros(size(X)));

% Loop over all angular samples.  Vectorised over the entire x–z grid for
% each (θ,φ) point because the spatial grid is usually smaller.
for ii = 1:numel(THETA)
    if ~mask_valid(ii), continue; end   % skip rays with TIR

    th  = THETA(ii);
    ph  = PHI(ii);
    ct  = cos_theta_water(ii);
    st  = sin_theta_water(ii);

    % Phase factor (Fresnel propagation to each (x,z) position)
    % k·r = k (x sinθ cosφ + y sinθ sinφ + z cosθ) ; here y = 0
    phase_term = exp(1i * k_w * ( X * st * cos(ph) + Z * ct ));

    % Vectorial Richards–Wolf weighting coefficients (Richards & Wolf 1959)
    % c_th = cosθ in image space (water) ; st = sinθ
    common  = Pupil(ii) * st * sqrt(ct) * dtheta * dphi;    % sqrt(ct) → apodisation
    % Components of the electric field integrand
    Ex = Ex + common * ( (1 + ct) * cos(ph)^2 + (1 - ct) * sin(ph)^2 ) .* phase_term;
    Ey = Ey + common * ( (1 + ct) * cos(ph) * sin(ph) - (1 - ct) * cos(ph) * sin(ph) ) .* phase_term;
    Ez = Ez - 2 * common * ct * cos(ph) .* phase_term;
end

% Global prefactor (Richards & Wolf)
Ex = -1i * k_w * f_obj / (2*pi) * Ex;
Ey = -1i * k_w * f_obj / (2*pi) * Ey;
Ez = -1i * k_w * f_obj / (2*pi) * Ez;

%% 7. NORMALISATION -------------------------------------------------------
% Normalise so that max(|E|) = 1 for easier visual comparison
E_mag = sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2);
Ex = Ex ./ max(E_mag(:));
Ey = Ey ./ max(E_mag(:));
Ez = Ez ./ max(E_mag(:));

%% 8. VISUALISATION -------------------------------------------------------

figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.8]);

subplot(2,3,1);
imagesc(x*1e6, z*1e6, abs(Ex));
axis image; title('|E_x|'); xlabel('x [µm]'); ylabel('z [µm]'); colorbar;

subplot(2,3,2);
imagesc(x*1e6, z*1e6, angle(Ex));
axis image; title('arg(E_x)'); xlabel('x [µm]'); ylabel('z [µm]'); colorbar;

subplot(2,3,3);
imagesc(x*1e6, z*1e6, abs(Ey));
axis image; title('|E_y|'); xlabel('x [µm]'); ylabel('z [µm]'); colorbar;

subplot(2,3,4);
imagesc(x*1e6, z*1e6, angle(Ey));
axis image; title('arg(E_y)'); xlabel('x [µm]'); ylabel('z [µm]'); colorbar;

subplot(2,3,5);
imagesc(x*1e6, z*1e6, abs(Ez));
axis image; title('|E_z|'); xlabel('x [µm]'); ylabel('z [µm]'); colorbar;

subplot(2,3,6);
imagesc(x*1e6, z*1e6, angle(Ez));
axis image; title('arg(E_z)'); xlabel('x [µm]'); ylabel('z [µm]'); colorbar;

sgtitle('Focused field inside water (oil–glass–water system)');

fprintf('Done.\n');
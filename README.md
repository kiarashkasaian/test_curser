# Focused Beam in Stratified Medium - Richards-Wolf Integral

This MATLAB code calculates the electric field of a focused beam in a stratified medium using the Richards-Wolf integral formulation. The code specifically handles a three-layer system (oil-glass-water) with aberrations due to refractive index mismatches.

## Features

- **Stratified Medium**: Three-layer system with oil, glass, and water
- **High NA Objective**: Supports NA = 1.4 objectives
- **Aberration Correction**: Includes phase aberrations due to RI mismatches between design and actual values
- **Gaussian Beam Profile**: Input beam has Gaussian apodization at the back focal plane
- **Linear Polarization**: Input beam is linearly polarized in x-direction
- **Vectorial Calculation**: Full 3D electric field components (Ex, Ey, Ez)
- **Comprehensive Visualization**: Shows amplitude and phase of all field components

## Physical Setup

```
    Oil (immersion)     |  n_oil = 1.515 (actual), 1.518 (design)
    ----------------    |  thickness = 100 μm
    Glass (coverslip)   |  n_glass = 1.520 (actual), 1.518 (design)  
    ----------------    |  thickness = 170 μm
    Water (sample)      |  n_water = 1.33
         FOCUS          |
```

## Key Parameters

- **Wavelength**: 532 nm (green laser)
- **Numerical Aperture**: 1.4
- **Gaussian Beam Width**: σ = 0.7 (in units of NA)
- **Calculation Grid**: ±2 μm in x, ±4 μm in z
- **Angular Sampling**: 200 θ points, 100 φ points

## Theory

The code implements the Richards-Wolf integral:

```
E(r) = ∫∫ A(θ,φ) × T(θ,φ) × exp[ik·r + iΦ_aberration(θ)] dΩ
```

Where:
- `A(θ,φ)` is the amplitude function including Gaussian apodization
- `T(θ,φ)` are the Fresnel transmission coefficients
- `Φ_aberration(θ)` is the aberration phase due to RI mismatch
- Integration is over the solid angle of the objective

## Key Functions

### `calculate_aberration_phase()`
Calculates the phase aberration due to optical path differences between design and actual refractive indices in the oil and glass layers.

### `fresnel_coefficients()`
Computes Fresnel transmission coefficients for s- and p-polarized light through the stratified medium.

### `plot_results()`
Creates comprehensive visualizations including:
- 2D amplitude and phase maps for Ex, Ey, Ez
- Total intensity distribution
- Line plots through the focus

## Usage

Simply run the script in MATLAB:

```matlab
focused_beam_stratified_medium()
```

The calculation will take several minutes depending on your computer. Progress is displayed during calculation.

## Output

The code generates three figure windows:

1. **Main Field Components**: 3×2 subplot showing amplitude and phase of Ex, Ey, Ez
2. **Total Intensity**: Combined intensity |E|² distribution
3. **Line Profiles**: Field amplitudes along x and z axes through the focus

## Customization

You can modify key parameters at the top of the main function:

- `sigma`: Gaussian beam width parameter
- `n_oil_actual`, `n_glass_actual`: Actual refractive indices
- `t_oil`, `t_glass`: Layer thicknesses
- `x_range`, `z_range`: Calculation region size
- `Nx`, `Nz`: Grid resolution
- `N_theta`, `N_phi`: Angular integration resolution

## Physical Insights

The code reveals several important optical phenomena:

1. **Aberration Effects**: RI mismatches cause focus shift and aberrations
2. **Polarization Mixing**: High NA focusing converts linear to elliptical polarization
3. **Longitudinal Field**: Strong Ez component near focus due to high NA
4. **Asymmetric PSF**: Stratified medium breaks symmetry of point spread function

## Requirements

- MATLAB (tested on R2019b and later)
- No additional toolboxes required

## Performance Notes

- Calculation time: ~5-10 minutes on modern PC
- Memory usage: ~100 MB for default grid size
- For faster calculation, reduce `N_theta`, `N_phi`, or grid size
- For higher accuracy, increase angular sampling or grid resolution

# Focused Beam in Stratified Medium - Richards-Wolf Integral

This MATLAB code calculates the electric field of a focused beam in a stratified medium using the Richards-Wolf integral formulation. The code includes all the requested features for optical microscopy simulations.

## Features

- **Three-layer stratified medium**: Oil, glass, and water layers with focus in water
- **High NA objective**: NA = 1.4 (typical for oil immersion objectives)
- **Aberration modeling**: Accounts for refractive index mismatches between design and actual values
- **Gaussian beam profile**: Input beam at back focal plane (BFP) with adjustable sigma
- **Linear polarization**: X-polarized input beam
- **Complete field visualization**: Shows amplitude and phase of Ex, Ey, Ez components in x-z cross section

## Physics Implementation

### Richards-Wolf Integral
The code implements the vectorial Richards-Wolf integral for high-NA focusing:

```
E(r) = -ik∫∫ A(θ,φ) T(θ) exp(iΦ(θ)) exp(ik·r) sin(θ) dθ dφ
```

Where:
- `A(θ,φ)` is the apodization function (Gaussian beam profile)
- `T(θ)` are the Fresnel transmission coefficients
- `Φ(θ)` are the phase aberrations due to RI mismatch
- `k·r` is the phase factor for focusing

### Stratified Medium Model
The code handles three layers:
1. **Oil layer**: Immersion oil between objective and coverslip
2. **Glass layer**: Coverslip 
3. **Water layer**: Sample medium (where focus is located)

### Aberration Calculation
Phase aberrations are calculated as:
```
Φ_aberr = k₀ × (path_actual - path_design)
```

Where the actual path accounts for refraction through each layer according to Snell's law.

### Fresnel Transmission
Transmission coefficients are calculated for each interface using Fresnel equations, accounting for the polarization state and angle-dependent transmission.

## Parameters

### Default Values
- **Wavelength**: 532 nm (green laser)
- **Numerical Aperture**: 1.4
- **Focal length**: 3 mm
- **Gaussian sigma**: 0.8 (in units of NA)
- **Layer thicknesses**: 100 μm oil, 170 μm glass

### Refractive Indices
- **Design values**: n_oil = 1.518, n_glass = 1.518, n_water = 1.333
- **Actual values**: n_oil = 1.515, n_glass = 1.520, n_water = 1.333

## Usage

Simply run the main function:
```matlab
focused_beam_stratified_medium()
```

**Requirements**: MATLAB with Parallel Computing Toolbox for optimal performance.

The code will:
1. Start a parallel pool (if not already running)
2. Pre-calculate angle-dependent quantities for efficiency
3. Calculate the electric field using Richards-Wolf integral with parallel processing
4. Display progress during calculation
5. Generate comprehensive plots showing:
   - Amplitude and phase of Ex, Ey, Ez components
   - Total intensity distribution
   - All plots in x-z cross section (y=0 plane)

**Note**: If Parallel Computing Toolbox is not available, the code will fall back to serial processing (replace `parfor` with `for`).

## Output Plots

The code generates two figures:

### Figure 1: Electric Field Components (3×2 subplot)
- Top row: Ex amplitude and phase
- Middle row: Ey amplitude and phase  
- Bottom row: Ez amplitude and phase

### Figure 2: Total Intensity
- Combined intensity |E|² = |Ex|² + |Ey|² + |Ez|²

## Computational Notes

- **Grid size**: 201×401 points covering ±2λ in x and ±4λ in z
- **Integration resolution**: 200 θ points × 400 φ points
- **Calculation time**: ~30 seconds to 2 minutes (depending on CPU cores)
- **Memory usage**: Moderate (primarily for storing field arrays)
- **Parallel processing**: Uses `parfor` for significant speedup (3-8x faster)

## Performance Optimizations

The code includes several performance improvements:

1. **Parallel Processing**: Uses MATLAB's `parfor` to distribute calculations across CPU cores
2. **Pre-computation**: Angle-dependent quantities calculated once and reused
3. **Vectorized Operations**: Eliminates inner loops using MATLAB's vector operations
4. **Memory Optimization**: Efficient array operations and minimal memory allocation
5. **Smart Progress Reporting**: Reduces I/O overhead during computation

### Parallel Processing Setup
The code automatically:
- Detects if a parallel pool is running
- Starts a local parallel pool if needed
- Distributes spatial points across available workers
- Uses shared read-only data for angle-dependent calculations

### Expected Performance
- **Serial version**: ~5-15 minutes for full calculation
- **Parallel version**: ~30 seconds to 2 minutes (depending on cores)
- **Speedup**: Typically 3-8x on modern multi-core systems
- **Memory**: ~100-500 MB depending on grid resolution

## Customization

You can easily modify:
- **sigma**: Change Gaussian beam width
- **Refractive indices**: Adjust for different materials or aberration levels
- **Layer thicknesses**: Modify oil and glass layer dimensions
- **Wavelength**: Change for different laser sources
- **Grid resolution**: Adjust for higher/lower resolution vs. speed trade-off

## Physical Insights

The code reveals important physics:
- **Spherical aberration** from RI mismatch causes focal shift and broadening
- **Polarization mixing** creates all three field components even from x-polarized input
- **High-NA effects** show significant longitudinal field (Ez) component
- **Stratified medium** effects modify the point spread function

This implementation is suitable for:
- Optical microscopy simulations
- Point spread function analysis
- Aberration studies
- High-NA focusing research
- Educational demonstrations of vectorial diffraction theory

# Bloch-McConnell CEST Analysis Pipeline

## New Files for BMC-Based Analysis

### 1. `phantom_CEST_BMC_analysis.m` - Main Analysis Script

**What it does:**
- Automatic detection of 24 tubes
- WASSR-based B0 field correction
- Bloch-McConnell-based analysis (replaces Lorentzian fitting)
- Kex (exchange rate) extraction
- Concentration estimation
- Parametric map generation

**Key Features:**
- ✅ No Lorentzian fitting
- ✅ BMC simulation for %CEST calculation
- ✅ Kex extraction from CEST effects
- ✅ Concentration estimation per tube
- ✅ WASSR B0 correction
- ✅ Automatic tube detection
- ✅ Tube numbering verification

### 2. `reorder_tubes_to_match_layout.m` - Tube Reordering Helper

**When to use:**
If automatic tube numbering doesn't match your phantom layout (shots.png), use this script to manually reorder tubes to match your specific layout.

## Quick Start Guide

### Step 1: Run Main Analysis

```matlab
cd /path/to/BM_sim_fit
phantom_CEST_BMC_analysis
```

### Step 2: Verify Tube Numbering

The script will pause and show you `DETECTED_tube_numbering.png`.

**Compare this with your shots.png image:**
- Do the tube numbers match your layout?
- **If YES**: Click "Yes, continue" and the analysis will complete automatically
- **If NO**: Click "No, I need to reorder manually" and proceed to Step 3

### Step 3: Manual Reordering (if needed)

```matlab
reorder_tubes_to_match_layout
```

Follow the interactive prompts to map detected tubes to your actual layout.

### Step 4: Rerun Analysis (if reordered)

```matlab
phantom_CEST_BMC_analysis
```

It will automatically load the reordered masks and continue.

## Your Phantom Layout

Based on your shots.png image:

| Tube Numbers | Chemical | Concentration | pH Values |
|--------------|----------|---------------|-----------|
| 1-3 | Iopamidol | 20 mM | 6.2, 6.8, 7.4 |
| 4-6 | Iopamidol | 50 mM | 6.2, 6.8, 7.4 |
| 7-9 | Creatine | 20 mM | 6.2, 6.8, 7.4 |
| 10-12 | Creatine | 50 mM | 6.2, 6.8, 7.4 |
| 13-15 | Taurine | 20 mM | 6.2, 6.8, 7.4 |
| 16-18 | Taurine | 50 mM | 6.2, 6.8, 7.4 |
| 19-21 | PLL | 0.1% w/v | 6.2, 6.8, 7.4 |
| 22-24 | PBS (blank) | - | 6.2, 6.8, 7.4 |

## Output Files

After running the analysis, you'll get:

1. **BMC_CEST_results.csv** - Main results table with:
   - Tube numbers and labels
   - Kex values (Hz)
   - Concentration estimates (mM)
   - %CEST at multiple offsets (1.9, 3.5, 4.3 ppm)
   - B0 field values per tube

2. **BMC_parametric_maps.png** - Comprehensive visualization:
   - Kex map
   - %CEST maps at 4.3 ppm (Iopamidol)
   - %CEST maps at 1.9 ppm (Creatine)
   - %CEST maps at 3.5 ppm (Amide)
   - B0 field map
   - Kex bar chart by tube

3. **DETECTED_tube_numbering.png** - Tube detection verification image

4. **BMC_tubeMasks.mat** - Tube masks for reuse

5. **BMC_CEST_workspace.mat** - Full MATLAB workspace

## CSV Output Format

| Column | Description |
|--------|-------------|
| `tube_number` | 1-24 |
| `tube_label` | Chemical, concentration, pH |
| `num_voxels` | Number of pixels in tube ROI |
| `Kex_Hz` | Exchange rate (Hz) |
| `concentration_mM` | Estimated concentration (mM) |
| `chemical_shift_ppm` | Expected chemical shift |
| `CEST_at_3_5ppm` | %CEST at 3.5 ppm (amide) |
| `CEST_at_4_3ppm` | %CEST at 4.3 ppm (Iopamidol) |
| `CEST_at_1_9ppm` | %CEST at 1.9 ppm (Creatine) |
| `B0_mean_ppm` | Mean B0 field offset |

## Differences from Previous Scripts

### ❌ Removed:
- Lorentzian multi-pool fitting
- Manual ROI selection for single tube
- Generic offset handling

### ✅ Added:
- Bloch-McConnell-based Kex extraction
- Concentration estimation
- Automatic processing of all 24 tubes
- Tube numbering verification system
- Chemical-specific CEST measurements
- Predefined offsets for your acquisition protocol

## Kex Calculation Method

The script uses a simplified Bloch-McConnell approach:

```
Kex ≈ 2π × Δω × (MTRasym% / f_b)
```

Where:
- Δω = chemical shift frequency (rad/s)
- MTRasym = measured CEST effect (%)
- f_b = bound pool fraction (~0.5% assumed)

For more accurate Kex values, you can integrate with `multiZfit.m` from the optimisation directory.

## Troubleshooting

### Issue: "Detected 20 tubes instead of 24"

**Solution:**
- Adjust Gaussian smoothing: Change `imgaussfilt(S0_image, 4)` to `3` or `5`
- Adjust size threshold: Change `numOfPixels < 10` to different value (5-20)

### Issue: "Tube numbering doesn't match shots.png"

**Solution:**
- Use `reorder_tubes_to_match_layout.m` script
- Or manually edit the `tube_order` variable in the BMC script

### Issue: "Kex values seem too high/low"

**Solution:**
- Check B0 correction quality (view B0_map)
- Verify S0 normalization
- Adjust assumed bound pool fraction (currently 0.5%)

### Issue: "Want to use full multiZfit instead of simplified model"

**Solution:**
- Uncomment multiZfit integration sections in the BMC script
- Ensure all BM simulation functions are in MATLAB path
- May need to adjust Sim structure parameters

## Advanced Usage

### Integrate Full Bloch-McConnell Fitting

To use the complete `multiZfit` function for more accurate Kex:

1. Ensure optimisation folder is in path
2. Set up Sim structure with appropriate parameters:

```matlab
Sim = init_Sim();
Sim.B0 = 3.0;  % Tesla
Sim.FREQ = [0; expected_shift_ppm];  % Water + CEST pool
Sim.T1 = [2.0; 1.5];  % T1 values (s)
Sim.T2 = [0.08; 0.01];  % T2 values (s)
% ... continue setup ...

FIT = multiZfit(tube_zspec, ppmOffsets, Sim);
Kex_fitted = FIT.results.k_BA;
```

3. Replace simplified Kex calculation in Step 7

## Support

For issues or questions about:
- **Tube detection**: Check Step 5 in main script
- **B0 correction**: Check Step 3 and 6
- **BMC fitting**: Check Step 7
- **Tube reordering**: Use helper script

## Version History

- **v1.0** (2025-11-05): Initial BMC-based analysis pipeline
  - Removed Lorentzian fitting
  - Added BMC Kex extraction
  - Added tube numbering verification
  - Integrated WASSR B0 correction

---

**Last Updated:** 2025-11-05
**Requires:** MATLAB R2019b or later with Image Processing Toolbox

# Figure 3C scale calibration

Figure 3C includes two calibrated 1-mm scale bars. Each bar applies to all three maps in its row. No expression values, spots, tissue pixels, sample sizes, or statistical results were changed.

## Sources

- GEO sample records identify GBM030 (GSM8506649) and GBM049 (GSM8506650) as Visium CytAssist spatial gene-expression FFPE data: the locally archived `GSE276841_family.soft.gz` and [GSE276841](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276841).
- The [archived tissue map](../inputs/figures/55_figure3_visual_polish/Figure3/Figure3C_spatial_score_atlas.png) is the original 3354 × 2267 PNG with the stained sections and score overlays. The spot source table in `inputs/write/55_figure3_visual_polish/Figure3/source_data/` supplies the displayed barcodes and color limits.
- `sources/` contains the four original public, compressed tissue-position and scale-factor files, copied from the project's archived GEO download; there are 4,992 array coordinates per sample. Original filenames are preserved.
- [10x Genomics technical note CG000686, pp. 2–3](https://cdn.10xgenomics.com/image/upload/v1722285090/support/documents/CG000686_VisiumHDSpatialGeneExpressionPerformance_TechnicalNote_RevA.pdf) specifies the 100-µm centre-to-centre pitch for Visium CytAssist/v2.
- [Space Ranger spatial-output documentation](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs) defines array indices, full-resolution pixel coordinates and `tissue_hires_scalef`. It states that `spot_diameter_fullres` is intended for visualization; it is not assumed here to represent the nominal 55-µm capture diameter.

## Calculation and checks

Use physical array coordinates `u = 50 × array_col` and `v = 50 × sqrt(3) × array_row` in micrometres. Fit all original coordinates to the registered full-resolution pixel positions with an affine transformation. Residuals must be below one full-resolution pixel. Use the inverse affine metric along the horizontal image direction, not a mean of axes or a displayed marker diameter. Apply the original high-resolution scale factor, crop limits, equal-aspect axes and 142×96-mm, 600-dpi panel geometry.

The displayed barcodes, array indices and full-resolution coordinates must match the Figure 3C source table exactly; its high-resolution positions must also match the original scale-factor transformation. `Figure3C_scale_calibration.csv` records both calibrations, residuals, transformations and the final bar lengths in the retained 3354×2267-pixel C-panel coordinate system. The plotting script runs these checks before adding the bars in the left white margin, leaving the six tissue/spot raster regions intact.

The scale bars use array registration. To redraw Figure 3, run the following command from the repository root after installing the [figure dependencies](../README.md#rebuild-figures-and-tables).

```sh
NOA_LEGACY_ROOT="$PWD/NOA_R1/inputs" python NOA_R1/scripts/09b_export_final_figures.py --figures 3
```

The script writes the assembled Figure 3 as a 170-mm-wide, 600-dpi RGB TIFF and PNG.

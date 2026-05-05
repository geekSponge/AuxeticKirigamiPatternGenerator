# Kirigami Bistable Auxetic Pattern Generator

This project implements parametric laser-cutting patterns for bistable auxetic kirigami structures inspired by:
Rafsanjani, A., & Pasini, D. (2016). *Bistable auxetic mechanical metamaterials inspired by ancient geometric motifs*. Extreme Mechanics Letters, 9, 291-296. https://doi.org/10.1016/j.eml.2016.09.001

Six patterns are available:
- Square + tilted
- Square + circular
- Square + parallel
- Triangular + tilted
- Triangular + circular
- Triangular + parallel

The tilted-pattern tiling logic was also compared against the public FabAcademy/ULB workshop generator:
Fablab ULB. (2018). *FAB14 Workshop: Digital Mechanical Material by Design - Bistable Auxetic Kirigami*. https://fabacademy.org/2018/labs/fablabulb/FAB14_workshop_bistable.html

Code was developed with assistance from OpenAI Codex/ChatGPT.




## Files

- `generate_kirigami_pattern.py`  
  Main command-line generator for all six pattern categories.

- `kirigami_geometry_io.py`  
  Shared geometry, clipping, SVG, DXF, and PNG export utilities.

- `generate_verification_figures.py`  
  Optional script that makes labeled verification figures and clean 50 x 50 mm examples for all six categories.

- `requirements.txt`  
  Python dependency for PNG previews and verification figures.


## Main Command Pattern

Automatic naming mode:

```powershell
python .\generate_kirigami_pattern.py `
  --family square `
  --motif tilted `
  --width 50 `
  --height 50 `
  --cell-size 20 `
  --t 1 `
  --a-over-l 0.5
```

This writes SVG, DXF, and PNG preview files to `.\output` with a meaningful filename containing the settings.

You can also choose a different automatic output folder:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif parallel --width 60 --height 60 --cell-size 20 --t 1 --a-over-l 0.5 --output-dir .\my_patterns
```

Manual file naming mode:

```powershell
python .\generate_kirigami_pattern.py `
  --family square `
  --motif tilted `
  --width 50 `
  --height 50 `
  --cell-size 20 `
  --t 1 `
  --a-over-l 0.5 `
  --svg .\output\square_tilted.svg `
  --dxf .\output\square_tilted.dxf `
  --preview-png .\output\square_tilted_preview.png
```

For laser cutting, use the `.svg` or `.dxf`. The preview PNG is only for checking by eye.

## Common Parameters

These work for all six categories:

| Parameter | Meaning |
|---|---|
| `--family square` or `--family triangular` | Base grid and rotating-unit family |
| `--motif tilted`, `circular`, or `parallel` | Cut motif category |
| `--width` | Sheet width in mm |
| `--height` | Sheet height in mm |
| `--cell-size` | Building-block side length `l` in mm |
| `--t` | Finite hinge / ligament gap parameter in mm |
| `--a-over-l` | Ideal rotating-unit side-length ratio `a/l` |
| `--offset-x`, `--offset-y` | Optional pattern phase shift in mm |
| `--cut-width` | SVG stroke width in mm, default `0.2` |
| `--output-dir` | Automatic output folder, default `.\output` |
| `--svg` | Output SVG path |
| `--dxf` | Output DXF path |
| `--preview-png` | Optional preview image path |
| `--handle` | Add centered rectangular handle tab(s) to the outer contour |
| `--handle-side` | `top`, `bottom`, `left`, `right`, `top-bottom`, `left-right`, `both`, or `all`; default `both` |
| `--handle-width` | Handle tab width in mm; default is 25% of the shorter sheet side |
| `--handle-height` | Handle tab protrusion in mm; default is 15% of the shorter sheet side |
| `--no-outline` | Do not include rectangular sheet outline |

If `--svg`, `--dxf`, and `--preview-png` are all omitted, the script automatically writes all three output types.

## Optional Handle Tabs

By default the outer contour is the rectangle around the patterned sheet. Add `--handle` when you want the contour to include tab handles for gripping or mounting. The internal kirigami cuts are unchanged; only the outer `OUTLINE` contour changes.

Centered top and bottom tabs using automatic tab size:

```powershell
python .\generate_kirigami_pattern.py --family square --motif tilted --width 60 --height 60 --cell-size 20 --t 1 --a-over-l 0.5 --handle
```

Centered top and bottom tabs with explicit dimensions:

```powershell
python .\generate_kirigami_pattern.py --family square --motif tilted --width 60 --height 60 --cell-size 20 --t 1 --a-over-l 0.5 --handle --handle-width 12 --handle-height 8
```

Top tab only:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif parallel --width 60 --height 60 --cell-size 20 --t 1 --a-over-l 0.5 --handle --handle-side top
```

Left and right tabs:

```powershell
python .\generate_kirigami_pattern.py --family square --motif tilted --width 60 --height 60 --cell-size 20 --t 1 --a-over-l 0.5 --handle --handle-side left-right --handle-width 12 --handle-height 8
```

Tabs on all four sides:

```powershell
python .\generate_kirigami_pattern.py --family square --motif tilted --width 60 --height 60 --cell-size 20 --t 1 --a-over-l 0.5 --handle --handle-side all --handle-width 12 --handle-height 8
```

## Motif-Specific Parameters

### Tilted Motifs

Optional parameter:

```text
--theta-deg
```

If `--theta-deg` is omitted, the code derives `theta` from `a/l`.

Definitions:

| Type | Relation |
|---|---|
| Square + tilted | `a/l = cos(theta) - sin(theta)` |
| Triangular + tilted | `a/l = cos(theta) - sqrt(3) sin(theta)` |

Typical paper/manufacturing values:

```text
Square tilted: theta about 20 deg to 24 deg
Triangular tilted: theta about 12 deg to 15.5 deg
```

### Circular Motifs

Optional parameter:

```text
--radius-over-l
```

If `--radius-over-l` is omitted, the code derives `R/l` from `a/l`.

Definitions:

| Type | Relation |
|---|---|
| Square + circular | `a/l = (-1 + sqrt(4(R/l)^2 - 1)) / sqrt(2)` |
| Triangular + circular | `a/l = (-1 + sqrt(12(R/l)^2 - 3)) / 2` |

Useful defaults for `a/l = 0.5`:

```text
Square circular: R/l approximately 0.989
Triangular circular: R/l approximately 0.764
```

### Parallel Motifs

Optional parameter:

```text
--width-over-l
```

If `--width-over-l` is omitted, the code derives `w/l` from `a/l`.

Definitions:

| Type | Relation |
|---|---|
| Square + parallel | `a/l = (1 - 2w/l) / sqrt(2)` |
| Triangular + parallel | `a/l = 1 - 2sqrt(3)(w/l)` |

Useful defaults for `a/l = 0.5`:

```text
Square parallel: w/l approximately 0.146
Triangular parallel: w/l approximately 0.144
```

## Example Commands For All Six Types

Create an output folder:

```powershell
mkdir .\output
```

### 1. Square + Tilted

Using `a/l`:

```powershell
python .\generate_kirigami_pattern.py --family square --motif tilted --width 50 --height 50 --cell-size 20 --t 1 --a-over-l 0.5 --svg .\output\square_tilted.svg --dxf .\output\square_tilted.dxf --preview-png .\output\square_tilted.png
```

Using direct angle:

```powershell
python .\generate_kirigami_pattern.py --family square --motif tilted --width 50 --height 50 --cell-size 20 --t 1 --theta-deg 20 --svg .\output\square_tilted_theta20.svg --dxf .\output\square_tilted_theta20.dxf --preview-png .\output\square_tilted_theta20.png
```

### 2. Square + Circular

Using `a/l`:

```powershell
python .\generate_kirigami_pattern.py --family square --motif circular --width 50 --height 50 --cell-size 20 --t 1 --a-over-l 0.5 --svg .\output\square_circular.svg --dxf .\output\square_circular.dxf --preview-png .\output\square_circular.png
```

Using direct `R/l`:

```powershell
python .\generate_kirigami_pattern.py --family square --motif circular --width 50 --height 50 --cell-size 20 --t 1 --radius-over-l 0.989 --svg .\output\square_circular_R.svg --dxf .\output\square_circular_R.dxf --preview-png .\output\square_circular_R.png
```

### 3. Square + Parallel

Using `a/l`:

```powershell
python .\generate_kirigami_pattern.py --family square --motif parallel --width 50 --height 50 --cell-size 20 --t 1 --a-over-l 0.5 --svg .\output\square_parallel.svg --dxf .\output\square_parallel.dxf --preview-png .\output\square_parallel.png
```

Using direct `w/l`:

```powershell
python .\generate_kirigami_pattern.py --family square --motif parallel --width 50 --height 50 --cell-size 20 --t 1 --width-over-l 0.146 --svg .\output\square_parallel_w.svg --dxf .\output\square_parallel_w.dxf --preview-png .\output\square_parallel_w.png
```

### 4. Triangular + Tilted

Using `a/l`:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif tilted --width 50 --height 50 --cell-size 20 --t 1 --a-over-l 0.5 --svg .\output\triangular_tilted.svg --dxf .\output\triangular_tilted.dxf --preview-png .\output\triangular_tilted.png
```

Using direct angle:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif tilted --width 50 --height 50 --cell-size 20 --t 1 --theta-deg 12 --svg .\output\triangular_tilted_theta12.svg --dxf .\output\triangular_tilted_theta12.dxf --preview-png .\output\triangular_tilted_theta12.png
```

### 5. Triangular + Circular

Using `a/l`:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif circular --width 50 --height 50 --cell-size 20 --t 1 --a-over-l 0.5 --svg .\output\triangular_circular.svg --dxf .\output\triangular_circular.dxf --preview-png .\output\triangular_circular.png
```

Using direct `R/l`:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif circular --width 50 --height 50 --cell-size 20 --t 1 --radius-over-l 0.764 --svg .\output\triangular_circular_R.svg --dxf .\output\triangular_circular_R.dxf --preview-png .\output\triangular_circular_R.png
```

### 6. Triangular + Parallel

Using `a/l`:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif parallel --width 50 --height 50 --cell-size 20 --t 1 --a-over-l 0.5 --svg .\output\triangular_parallel.svg --dxf .\output\triangular_parallel.dxf --preview-png .\output\triangular_parallel.png
```

Using direct `w/l`:

```powershell
python .\generate_kirigami_pattern.py --family triangular --motif parallel --width 50 --height 50 --cell-size 20 --t 1 --width-over-l 0.144 --svg .\output\triangular_parallel_w.svg --dxf .\output\triangular_parallel_w.dxf --preview-png .\output\triangular_parallel_w.png
```

## Generate Verification Figures

This command generates labeled figures and clean 50 x 50 mm examples for all six types:

```powershell
python .\generate_verification_figures.py
```

The outputs will be written to:

```text
.\verification_results
```

## Notes

- `a` means the ideal rotating-unit side length.
- `t` is the finite uncut ligament / cut-gap parameter.
- For tilted motifs, finite `t` changes the visible cut endpoints, so the endpoint polygon is not always exactly equal to the ideal `a`.
- Keep the labeled PNGs for checking only. Use the clean `.svg` or `.dxf` files for laser cutting.

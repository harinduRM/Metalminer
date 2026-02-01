# MetalMiner CLI Options and Tags (Detailed)

This file documents all command-line flags and configuration options used by
MetalMiner, with detailed explanations of what each option does and how it
affects the pipeline. It is intended as a wiki-style reference for users.

## CLI entry point

```bash
python -m metalminer.cli [OPTIONS]
```

## Notes on list inputs

- Most list flags accept space-separated values (e.g., `--target-metals Pu Th Ce`)
- Comma-separated values are also accepted (e.g., `--target-metals Pu,Th,Ce`)
- Mixed spacing and commas are supported (e.g., `Pu, Th Ce`)

## Boolean flags

Many flags are defined as `BooleanOptionalAction`, so both forms exist:

- `--flag` to enable
- `--no-flag` to disable

## Primary inputs

### `--target-metals` (list[str])

**What it does:**
Defines which element symbols the CSD search targets. Each metal in the
list becomes a separate search and output set.

**When to use:**
Use this for broad searches (e.g., all structures containing Pu or Th).

**Interactions:**
If `--target-csd-ids` is provided, the explicit IDs are used instead of
search, but the metal list still sets output file naming.

**Default:** `["Pu"]`

### `--target-csd-ids` (list[str] | null)

**What it does:**
Bypasses the substructure search and processes only the given CSD
refcodes (entries) in the list.

**When to use:**
Use this to re-run on a known set of CSD IDs or for curated lists.

**Interactions:**
When provided, search is skipped and only those entries are loaded.

**Default:** `None`

## Extraction + filtering

### `--extraction-method` (str)

**What it does:**
Chooses the extraction strategy used to build the coordination
environment around the target metal.

**Typical values:**
"Topological" (default) or "Geometric".

**Interactions:**
Geometric mode uses `--geometric-radius` and cutoffs more heavily.

**Default:** `"Topological"`

### `--r-factor-limit` (float)

**What it does:**
Filters out structures with R-factor greater than or equal to this value.
Higher R-factor implies lower structural quality.

**When to use:**
Lower this for stricter quality control (e.g., 5).

**Default:** `10`

### `--process-limit` (int)

**What it does:**
Caps how many unique metal sites are processed in a run.

**When to use:**
Use a smaller number to test quickly or debug.

**Default:** `10000`

### `--filter-polymeric` / `--no-filter-polymeric` (bool)

**What it does:**
If enabled, polymeric structures are skipped.

**When to use:**
Enable for discrete complexes only; disable for extended networks.

**Default:** `True`

### `--filter-powder` / `--no-filter-powder` (bool)

**What it does:**
If enabled, powder diffraction structures are skipped.

**When to use:**
Disable if you want powder structures included.

**Default:** `True`

### `--filter-all-disorder` / `--no-filter-all-disorder` (bool)

**What it does:**
If enabled, structures with global disorder are removed.

**When to use:**
Enable to avoid complex disorder cases.

**Default:** `False`

### `--filter-primary-disorder` / `--no-filter-primary-disorder` (bool)

**What it does:**
If enabled, structures with primary-sphere disorder are removed.

**When to use:**
Enable to avoid complicated primary coordination environments.

**Default:** `False`

### `--correct-primary-disorder` / `--no-correct-primary-disorder` (bool)

**What it does:**
Attempts to resolve primary-sphere disorder before analysis.

**When to use:**
Enable for automatic cleanup; disable if you want raw disordered cases.

**Default:** `True`

## Extraction geometry controls

### `--geometric-radius` (float)

**What it does:**
Sets the radius (in Angstrom) for selecting neighboring atoms in
Geometric extraction.

**When to use:**
Increase if coordination shells are larger; decrease for stricter shells.

**Default:** `3.8`

### `--limit-nm-nm` (float)

**What it does:**
Cutoff for non-metal to non-metal distances used during pruning.

**Default:** `2.8`

### `--limit-m-nm` (float)

**What it does:**
Cutoff for metal to non-metal distances used during pruning.

**Default:** `3.5`

### `--limit-h-x` (float)

**What it does:**
Cutoff for hydrogen to X distances used during cleanup.

**Default:** `1.3`

### `--extraction-cutoff-json` (str | null)

**What it does:**
Overrides the three cutoff values in a single JSON string.

**Example:**
```
{"LIMIT_NM_NM":2.8,"LIMIT_M_NM":3.5,"LIMIT_H_X":1.3}
```

**When to use:**
Use this if you want to pass a full set of cutoffs at once.

**Default:** `None`

### `--num-metal-layers` (int)

**What it does:**
Controls how many metal layers deep are included in the coordination
environment analysis.

**Typical values:**
1 = primary sphere only, 2+ = include further metal shells.

**Default:** `1`

### `--hydrogen-addition-method` (str)

**What it does:**
Chooses how missing hydrogens are added during processing.

**Typical values:**
"Geometric".

**Default:** `"Geometric"`

### `--disorder-resolve-method` (str)

**What it does:**
Selects the strategy used to resolve disorder in the structure.

**Typical values:**
"Hybrid".

**Default:** `"Hybrid"`

### `--metalloligands` (list[str])

**What it does:**
List of element symbols that should be treated as metalloligands during
ligand analysis and oxidation-state evaluation.

**When to use:**
Modify this list if your chemistry includes unusual metalloligands.

**Default:** `["Cr", "V", "Mo", "W", "S", "Co"]`

## Oxidation states

### `--oxidation-states-filter` (list[str])

**What it does:**
Filters results to only these oxidation states. Use "all" to keep every
oxidation state.

**Examples:**
```
--oxidation-states-filter 3 4 5
--oxidation-states-filter all
```

**Default:** `["all"]`

### `--get-abstract` / `--no-get-abstract` (bool)

**What it does:**
Enables extraction and parsing of abstract text (when available) to
support oxidation state estimation.

**Implication:**
When disabled, only non-text methods are used.

**Default:** `True`

## Manual refinement

### `--edit-manual` / `--no-edit-manual` (bool)

**What it does:**
Runs the manual refinement step after processing (interactive or
user-assisted corrections, where available).

**When to use:**
Disable for fully automated runs.

**Default:** `True`

## Visualization

### `--visualize` (bool)

**What it does:**
Displays 3D visualizations during processing.

**When to use:**
Enable for debugging or inspection; disable for batch runs.

**Default:** `False`

### `--visualization-limit` (int)

**What it does:**
Limits the number of structures that will be visualized in a run.

**Default:** `15`

## Timeouts and retries

### `--site-timeout-seconds` (int | null)

**What it does:**
Sets a per-site timeout in seconds. If a site exceeds this, it is marked
as a timeout and can be retried later.

**Notes:**
Set to 0 or omit to disable timeouts.

**Default:** `None`

### `--retry-timeouts` / `--no-retry-timeouts` (bool)

**What it does:**
When resuming from a spool, timed-out sites are retried if enabled.

**Default:** `True`

## Typical usage example

```bash
python -m metalminer.cli \
  --target-metals Pa Am Cm Bk Cf \
  --oxidation-states-filter all \
  --extraction-method Topological \
  --r-factor-limit 5 \
  --process-limit 100000 \
  --no-visualize \
  --no-filter-polymeric \
  --filter-powder \
  --no-filter-all-disorder \
  --no-filter-primary-disorder \
  --correct-primary-disorder \
  --hydrogen-addition-method Geometric \
  --disorder-resolve-method Hybrid \
  --num-metal-layers 1 \
  --get-abstract \
  --limit-nm-nm 2.8 \
  --limit-m-nm 3.5 \
  --limit-h-x 1.3 \
  --metalloligands Cr V Mo W Co \
  --no-edit-manual \
  --site-timeout-seconds 240
```

## Python entry point equivalent

```python
from metalminer.main import main as run_metalminer

if __name__ == "__main__":
    run_metalminer(
        TARGET_METAL_list=["Pa", "Am", "Cm", "Bk", "Cf"],
        target_csd_ids=None,
        OXIDATION_STATES_FILTER=["all"],
        EXTRACTION_METHOD="Topological",
        R_FACTOR_LIMIT=5,
        PROCESS_LIMIT=100000,
        VISUALIZE=False,
        FILTER_POLYMERIC=False,
        FILTER_POWDER=True,
        FILTER_ALL_DISORDER=False,
        FILTER_PRIMARY_DISORDER=False,
        CORRECT_PRIMARY_DISORDER=True,
        Hydrogen_Addition_method="Geometric",
        DISORDER_RESOLVE_METHOD="Hybrid",
        num_metal_layers=1,
        GET_ABSTRACT=True,
        Extraction_Cut_off_distances={
            "LIMIT_NM_NM": 2.8,
            "LIMIT_M_NM": 3.5,
            "LIMIT_H_X": 1.3,
        },
        metalloligands=["Cr", "V", "Mo", "W", "Co"],
        Edit_manual=False,
        SITE_TIMEOUT_SECONDS=240,
    )
```

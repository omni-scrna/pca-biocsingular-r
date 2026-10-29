# pca-biocsingular-r

`BiocSingular`-backed PCA module for omnibenchmark scRNA pipelines. `BiocSingular::runSVD`
dispatches on a `BSPARAM`, so the SVD solver is a benchmark parameter:

| `--solver` | `BSPARAM` | algorithm |
|---|---|---|
| `exact` | `ExactParam()` | LAPACK full SVD, truncated to `k` |
| `random` | `RandomParam()` | Halko randomized SVD via `rsvd` |
| `irlba` | `IrlbaParam()` | implicitly-restarted Lanczos via `irlba` |

## Setup

```sh
pixi install
pixi run check        # prints OK if all runtime libraries load
pixi run export-env   # regenerates envs/pca-biocsingular-r.yml from pixi.toml
```

## Usage

```sh
pixi run Rscript pca.R \
  --output_dir <dir> \
  --name <name> \
  --normalized_selected_h5 <normalized_selected.h5> \
  --solver {exact|random|irlba} \
  --n_components <int> \
  --random_seed <int> \
  [--dense {true|false}] \
  [--deferred {true|false}] \
  [--blas_threads <int>] \
  [--oversampling <int>] \
  [--power_iters <int>]
```

Outputs in `<output_dir>`:

- `<name>_embedding.tsv`: `cell_id`, `PC1..PC{n_components}`, one row per cell
- `<name>_loadings.tsv`: `gene`, `PC1..`, one row per gene

## Parameters

- `--dense`: materialise the matrix before PCA. Only this path does real BLAS-3 work, so it is
  the one that interacts with `--blas_threads`.
- `--blas_threads`: set via `RhpcBLASctl`, overriding Snakemake's `OMP_NUM_THREADS=1`. `0` leaves
  it alone. Results differ at ~1e-15 across thread counts; compare with `all.equal`.
- `--oversampling` (`p`, default 10), `--power_iters` (`q`, default 2): `random` only, passed to
  `rsvd::rsvd`. Ignored with a warning for other solvers.
- `--deferred`: center as an operator instead of materialising the centered matrix. Numerically
  equivalent; only lowers peak memory for `random` (~660 → 450 MB on be1-fixture).
- `--random_seed`: used by `random` only. `exact` is deterministic; `irlba` converges to `tol`.

### `q` matters

rsvd defaults to `q=2`; scikit-learn's `randomized_svd` uses `n_iter=7` at these shapes.
Mean per-PC `|cor|` vs `exact` on be1-fixture (2000 × 1715, `k=50`):

| `random` config | PC1-20 | PC21-30 | PC31-40 | PC41-50 |
|---|---|---|---|---|
| `q=2, p=10` (default) | 0.98 | 0.619 | 0.275 | 0.115 |
| `q=7, p=10` | 1.00 | 0.994 | 0.866 | 0.222 |
| `q=7, p=20` | 1.00 | 0.998 | 0.954 | 0.365 |
| scanpy randomized | 1.00 | 0.996 | 0.710 | 0.302 |

The error at low `q` is bias, so averaging seeds does not fix it. Past PC30 the spectrum is flat
and no setting recovers individual components.

## Citation

See `CITATION.cff`. Please also cite `BiocSingular` (Lun) and the solver used: Baglama & Reichel
(2005) for `irlba`, Halko, Martinsson & Tropp (2011) for `random`.

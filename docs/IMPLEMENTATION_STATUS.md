# Implementation status and review items

## What is supplied

- Two existing R analysis scripts, copied unchanged from the retrieved files.
- SOP, README, exact input guide and metadata templates.
- R project file, dependency setup and syntax/package check script.
- Git ignore rules that keep source datasets and generated results out of routine commits.
- Source-file SHA-256 checksums for provenance.

## Items to resolve before claiming an end-to-end reproduction

| Item | Evidence in included scripts | Required action |
| --- | --- | --- |
| Species mismatch | SDM has nine species without gambiae; overlap has five including gambiae | Identify the source model for each overlap map and decide the intended analysis scope |
| Local paths | Original scripts use D:/ and E:/ paths | Edit settings on the local machine using the README |
| Predictor names | Several raster stems differ | Document explicit mapping and verify source equivalence |
| Aspect processing | SDM transforms before alignment; overlap fallback after alignment | Reconcile with the fitted model's actual preprocessing |
| Different overlap methods | SDM has PCA-env; separate script has SVM/centroid metrics | Keep methods, outputs and interpretations distinct |
| Classification settings | SVM 0.80 is exploratory; other cutoffs unset | Report continuous metrics or justify and record all categorical rules |
| Future scenarios | FUTURE is empty | Supply real scenario inputs before describing projected results |
| Data provenance | Input files and source metadata not supplied | Complete product/version/period/units/access records |
| R environment | No R binary or input data available during preparation | Run the included check locally and create a validated lockfile |
| Model execution | No models run during repository preparation | Perform preparation audit and actual local modelling; retain logs and review |
| Licensing and citation | Not assigned | Confirm with project contributors before release |

No universal threshold, exact spatial resolution, scenario period, final species scope or publication result has been inferred to fill these gaps. Retain the original scripts until the implementation is reconciled with the actual study version.

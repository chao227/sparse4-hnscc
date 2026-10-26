
# Known limitations

1. sparse4 is a post hoc parsimonious model. Candidate procedures and external
   cohort results had been inspected before the four-gene rule was adopted.
2. Internal validation is conditional on the historical full-cohort Yellow
   module membership and the locked full-cohort VST/top-5,000 realization; it
   is not end-to-end validation of module discovery.
3. External cohorts are post hoc evaluations, not pristine outcome-unseen
   validation cohorts. GSE159067 contains only two of four genes and is a
   partial exploratory evaluation.
4. The M1-M4 and calibration results are same-cohort apparent analyses and do
   not correct all model-development optimism.
5. Spatial results concern complete Yellow-module activity and TLS-associated
   expression scores. They are not a spatial validation of sparse4, direct
   pathology confirmation, or causal evidence.
6. Historical B-state labels remain provisional where reference agreement and
   marker detection were weak. Cross-sectional ordering is not lineage proof.
7. Public raw data and participant-level derivatives are not redistributed.
   A complete clean-machine full-data run therefore requires users to retrieve
   the public deposits and reconstruct the documented inputs.
8. Some exact historical package versions were not recoverable. The environment
   files distinguish currently observed versions from historical evidence.
9. Final Fig. 7A and Fig. 7D are author-edited schematics. Their final pixel
   assets and a claimed executable schematic generator are not distributed;
   the underlying numeric analyses remain mapped individually.
10. Data-intensive execution requires users to retrieve the public datasets,
    reconstruct the documented local inputs, and provide the computational
    resources required by the single-cell and spatial analyses.
11. Historical and secondary scripts are retained under `archive_exploratory/`
    for transparency but are not part of the ordered full workflow. Their
    presence documents the analysis history and does not authorize or imply a
    new execution of those procedures.
12. Several downstream scripts use workflow-specific reconstructed objects,
    including regulon, missingness/annotation, table-assembly, and figure-source
    inputs. Expected paths and formats are documented in the relevant scripts.

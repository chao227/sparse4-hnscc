
# Exploratory and historical code

These scripts were used during exploration, procedure comparison, or secondary
reconstruction and are retained for transparency. They do not define the
current locked four-gene score, do not produce the repository's final-results
workflow, and are not called by `run_all.sh`. Their execution order is not
guaranteed. The current score definition is controlled only by
`config/model_lock.json`, `config/coefficients.tsv`, and
`docs/sparse4_scoring_spec.md`.

The separate lambda-min procedure and its completed permutation evidence apply
only to that procedure. They are not inferential evidence for sparse4.

The archived donor-paired script records the historical computation underlying
an already reported uncertainty summary. It is retained for auditability, is
not part of the current executable workflow, and is never called by
`run_all.sh`.

# Reproduce Figure S16

This directory runs on its own. Install Python with `numpy`, `scipy` and
`pandas`, and R with `jsonlite`, `Matrix`, `matrixStats` and Cairo support.
Use Arial for the published figure's typography. From this directory, run

```sh
python reproduce.py
```

The script downloads both upstream versions from GitHub and checks their
SHA-256 digests. It uses the three ADSP LD matrices included here. The PDF
and the alpha and PIP results go to `output/`. Use `--out` to write them elsewhere. The downloaded code and intermediate files live in a temporary
directory.

## Comparison

Table S2 gives two ways to add the same $z$-score update to the 2019 RSS-lambda
code. In `rss_fits.R`, `rss2019_in_new_bottle()` calls the original function one
sweep at a time.
`rss2019_internal.R` adds the update inside a copy of the fitting loop. Both
implementations use the same starting values, priors and stopping rule as the
RSparsePro fits in this comparison.

We fix lambda at 0.1, 0.5, 1 and 5 in 36 cases with ADSP sketch, singular or
identity LD. We compare single effect inclusion probabilities (alpha) and
summarize each variant's PIP as the maximum across effects, rounded to four
decimal places as in RSparsePro.

## Original sources

- susieR **0.7.1.0487**, March 21, 2019, commit
  [`9837e100`](https://github.com/stephenslab/susieR/tree/9837e1007808afb4380ecdc4270f57bd847824c3).
  This is the first R implementation of the RSS-lambda model.
- RSparsePro commit
  [`5990896`](https://github.com/zhwm/RSparsePro_LD/tree/5990896d58053a78d4260f011aa83949ae770c0a).
  The Python source in The ADSP input contains three LD matrices from the Release 5 European
ancestry chromosome 22 sketch. Each uses all 10,000 stored projections.

Gao Wang wrote the benchmark with assistance from GPT-6.1 Sol in Codex for macOS.

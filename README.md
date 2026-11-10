# ychannel-adjoint-inverse

MATLAB code and data for the simultaneous regularized recovery of spatially varying diffusion
coefficients D(x, y) and source terms S(x, y) from NIR concentration images of an acid–base
reaction in a microfluidic Y-channel.

Code for the article *Simultaneous Regularized Field Recovery of Diffusion Coefficients and
Source Terms in Microfluidic Acid–Base Reactions* (R. Abdurrahman, K. Matsunaga, Y. Abe,
N. Kakuta), submitted to Chemical Engineering Science.

## Usage

Set `condition` (`'T28'`, `'T20'` or `'T37'`) in each script and run in MATLAB:

1. `step1_initial_estimates.m`: initial estimates of D and S (a few minutes)
2. `step2_cg_estimation.m`: adjoint-based CG recovery with FR, PR, PB and DY (several hours)
3. `step3_article_figures.m`: figures and numbers of the article

The default settings (300 CG iterations) reproduce the article for all three data sets.
Results are written to `results/<condition>/`.

## Contents

- `data/`: measured concentration maps (T20, T28, T37)
- `core/`: solver, adjoint and CG functions
- `tests/`: synthetic tests

## Requirements

MATLAB R2020a or newer. Step 1 also needs the Image Processing and Curve Fitting Toolboxes.

## License

[MIT](LICENSE)

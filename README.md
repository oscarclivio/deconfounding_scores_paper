1) Delete all files in `results/`; this folder contains individual runs' results and runs are not re-run if their corresponding results are already saved in the corresponding folder.
2) In each `bash run_[...].sh` file, in the line `cat "$COMBO_FILE" | parallel --line-buffer -j X --colsep ' ' '`, `X` is the number of parallel threads. Please adapt `X` to your computational resources. WARNING: HC-MNIST can be computationally intensive, thus the default lower value of `X` of 1 in `bash run_hcmnist.sh`.
3) Run each `bash run_[..].sh` file.
4) For simulated datasets, retrieve results through `Rscript ProcessResultsSimulated.R`. It prints a tabular given a table of results as in Table 1, and saves a figure as in Figure 2 in `combined_plot.pdf`. For semi-synthetic datasets, retrieve results through `Rscript ProcessResults.R`. It prints a tabular given a table of results as in Table 2. The tabulars can be visualized in LaTeX by defining `\newcommand{\red}[1]{\textcolor{red}{#1}}` and `\newcommand{\green}[1]{\textcolor{teal}{#1}}`.

Sources: 
- https://github.com/afranks86/reducer (original code)
- https://github.com/anndvision/quince (HCMNIST)

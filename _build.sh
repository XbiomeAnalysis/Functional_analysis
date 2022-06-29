#!/bin/sh

set -ev

/home/zhuzhengnong/software/HPC/miniconda3/envs/r4-base/bin/Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
/home/zhuzhengnong/software/HPC/miniconda3/envs/r4-base/bin/Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
/home/zhuzhengnong/software/HPC/miniconda3/envs/r4-base/bin/Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::epub_book')"


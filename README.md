# Denmark Forest Quality Projections Code Repository

**Summary reports, data visualisations and downlads can be found on:** https://jakobjassmann.github.io/dk_forest_lidar_2021/

Last update: March 2022.<br>

Jakob J. Assmann, Pil B. M. Pedersen, Jesper E. Moeslund, Cornelius Senf, Urs A. Treier, Derek Corcoran, Zsófia Koma, Signe Normand 

Code repository accompanying Assmann et al. (in prep) "Airborne laser scanning data fusion and machine learning identify forests of high conservation value in Denmark".

## Content

This repository contains all scripts and workflows required to regenerate the analysis presented in Assmann et al (in prep). 

- [1. Folder structure overview](#1-folder-structure-overview)
- [2. Required source data](#2-required-source-data)
- [3. Analysis workflow](#3-analysis-workflow)
- [4. Acknowledgements](#4-acknowledgements)
- [5. Citation](#5-citation)
- [6. License](#6-license) 

## 1. Folder structure overview

This code repostiory has the following folder sturcture:

```
data/                                   - Files and placeholder folders containing the data required for the analysis.
  |- models/                            - Compressed R data files containt the fitted models
  |- predictor_data/                    - Placeholder folder for the predictor data
  |- projections/                       - Placeholder folder for the nationwide projections
  |- response_data/                     - Annotated forest polygons on which the training dataset is based
  |- stratification/                    - Vector/raster source files for the stratifications tested (BIOWIDE and Derek)
  |- training_data/                     - Pixel centre coordinated and extracted predictor values for the model training
  |- variograms/                        - Variogram samples for predictor variables
docs/                                   - Figure outputs, R Markdown reports and website documents
  |- figures/                           - Figrues for the manuscript
  |- gbm_models_performance_files/      - Auxiliary files for the gbm_models_performance.html report
  |- js/                                - Folder containing JavaScript source files required for the webapp
  |- ranger_models_performance_files/   - Auxiliary files for the ranger_models_performance.html report
  |_ training_annotations_files/        - Auxiliary files for the training_annotation.html report
  |- variograms/                        - Plots for the predictor variograms
scripts/                                - Scripts for data prep, model fitting, projections and analysis
  |- figures/                           - Scripts to generate the manuscript figures
  |- models/                            - Scripts for model fitting
  |- prep/                              - Scripts for data prep
  |- projections/                       - Scripts for projecting the models across Denmark
```
[\[back to content\]](#content)

## 2. Required source data

[\[back to content\]](#content)

## 3. Analysis workflow

A conceptual overview of the workflow can be found [here](https://jakobjassmann.github.io/dk_forest_lidar_2021/workflow.html).<br><br>
The following steps describe the scripts contained in this reporsiory and the order in wich they need to be run.

*Note: All random numbers are generated as pseudo-random numbers using a seed to ensure reproducibility of the analysis.*

**Data preparation**

After download (see [Section 2.](#2-required-source-data) above), the polygon annotations, predictor data and forests masks need to be perpared as follows. 
<br><br>
All scripts are found in `scripts/prep`.

1. Data wrangling (Stage I): Run scripts `1a_*.R - 1d_*.R` to prepare: the forest mask (1a), the BIOWIDE stratification (1b), the foliage height diversity predictor (1c) and the tree type predictor (1d).
2. Data wrangling (Stage II): Run scripts `2a_*.R - 2c_*.R` to: reproject the predictors (2a), generate the variograms (2b - optional!), generate the focal predictors (2c). 
3. Pixel training data prep: Run scripts `3a_*.R and 3b_*.R `to: generate the pixel training sample (3a) and split the training sample into a 80/20 split using the geographic stratification (3b). 
4. Data wrangling (Stage III): Run scripts `4a_*.R and 4b_*.R` to: prepare the forest polygons for the leaflet we app (4a) and the forest change masks for the projections and web app (4b). 

**Model fitting**

The forest quality models are fitted and fine tuned using the scripts contained in `scripts/models`. The fitted models are saved in `data/models`. <br><br> All models presented in the final publication are contained as compressed R data files in this repository.

- Fit gradient boosting models (scripts: `/scripts/models/pixel_gbm_\*.R".
- Fit random forest models (scripts: `/scripts/models/pixel_ranger_\*.R".

**Projection across Denmark**

The models can then be used to project the forest quality predictions across Denmark. For the publication, we only used the models based on the "BIOWIDE" stratification, as these models performed better. The projection script also prepares the outputs as cloud optimised geotiffs for webhosting and visualisation using the web app. 

- Projection script: `scripts/projections/project_pixels.R`

**Web app and visualisation**

The JS and HTML source code for the HTML file that generates the web app can be found in [docs/data_viz.html](https://github.com/jakobjassmann/dk_forest_lidar_2021/blob/main/docs/data_vis.html). We have documented the basic principles of how it works in [this guide](https://jakobjassmann.github.io/dk_forest_lidar_2021/cog_guide.html). 

[\[back to content\]](#content)

## 4. Acknowledgements

We would like to acknowledge our funders. We woul also like to thank Thomas Nord Larsen, Bjørn Ole Ejlersen, Bjarne Aabrandt Jensen and Ane Brunbjerg for sharing tabular, vector and raster data with us that were key to faciliating the model predictions.

[\[back to content\]](#content)

## 5. Citation

Jakob J. Assmann, Pil B. M. Pedersen, Jesper E. Moeslund, Cornelius Senf, Urs A. Treier, Derek Corcoran, Zsófia Koma, Signe Normand. (in prep). "Airborne laser scanning data fusion and machine learning identify forests of high conservation value in Denmark". https://github.com/jakobjassmann/dk_forest_lidar_2021

[\[back to content\]](#content)

## 6. License

All code (.R files) are licensed under a modified MIT license (see license text [ADD]()).<br>
All other content of this repository, including: figures, raster, tabular and vector data are licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a><a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a>

When using code or other content from this repostiory in scientific work, please use the citation provided above to acknwoledge the authors. 

[\[back to content\]](#content)

# SnowTrendsAnalysis

In the following repository are included all the data, plots, codes and operations used for the paper "Diverging snowfall trends across months and elevation in the northeastern Italian Alps" by Bertoldi et al. 2023. Please note that when the data are used, both the dataset and the paper should be cited. See the section **DOI** to access to the dataset and the full article. 

The repository contains:
- **PreliminaryData_1**, **PreliminaryData_2**: *.zip* folders with all the preliminary operations on the raw data, as the data quality check and the gap filling procedure. Plots of the results coming from the data analysis are also included
- **FinalData**: folder with the final version data of fresh snow (HN), snow depth (HS), precipitation (P), mean air temperature (TMEAN), north atlantic oscillation index (NAO) and arctic oscillation index (AO). *_S* stays for seasonal data. In addition it is also included the **TopographicData** folder, with *DEM* and *.shp* files of the Trentino - South Tyrol region
- **R**: folder with all the R codes used for the analysis. Each step corresponds to a different R script. In every script can be found some comments that helps to understand better the code. To make the statistical analyses, many functions were used: *pettitt.test*, *sens.slope*, *mk.test* (included in *trend* package) and *cor.test*, *anova* (included in *stats* package)
- **Plots**: folder with the results coming from the analysis of the trends and the attribution validation. Each folder presents the result of each R script in the **R** folder. It also included in the folder **QGIS** the two maps of the HN climatology and study area made using the QGIS software

### Data Analysis
The data analysis was suddived in two parts. The first part was entirely dedicated to the collection and manipulation of fresh snow data in Trentino - South Tyrol. In particular, a monthly aggregation was performed starting from daily data. The dataset was then updated to 2019/2020 season and subjected to a gap filling procedure. At this point, the dataset was analysed and a quality check procedure was conducted.

Instead, the second part was devoted to the real analysis of the data. More in detail, first of all a trend analysis was conducted using statistical tests such as the Mann-Kendall trend test (Mann 1945; Kendall 1975) and the Sen's slope estimator (Sen 1968). Then, the trend attribution was evaluated, first using a correlation analysis and then an explained variance analysis.

### R code
The code was organized in different scripts. Each script made a particular operation. *00_Data.R* refers to all the preliminary data operations and analyses contained in **PreliminaryData_1** and **PreliminaryData_2**. Scripts from *01_TrendHN.R* to *05_Attribution.R* contain all the trend computations and plotting operations. *06_CorAnalysis.R* refers to the correlation analysis conducted using the Pearson correlation coefficient (Freedman et al. 2007). *07_GLMAnalysis.R*, *08_R2Analysis.R* are used as test scripts to have a comparison between a multiple linear regression model (v1) and a general regression model (family gamma) (v2) applied to HN using the drivers as dependent variables. *09_EVAnalysis.R* computes the one-way analysis of variance (ANOVA) for a multiple linear regression model (Miller 1997). *99_AvgHN.R* is used to calculate the HN climatology.

### Plots
The outputs coming from the computations contained in the R code are organized in the same structure of the code, so that every folders contained the correspond result of the script. The only exception regarts the data analysis results that are stored in the preliminary *.zip* folders.



_____
#### References
Freedman, D., Pisani, R., and Purves, R. (2007). “Statistics (international student edition)”. In: Pisani, R. Purves, 4th edn.WWNorton & Company, New York.

Kendall, M.G. (1975). "Rank correlation methods". Griffin, London, UK.

Mann, H.B. (1945). “Nonparametric tests against trend”. In: Econometrica: Journal of the econometric society, pp. 245–259.

Miller, R.G. (1997). "Beyond ANOVA: basics of applied statistics". CRC press.

Sen, P. K. (1968). “Estimates of the regression coefficient based on Kendall’s tau”. In: Journal of the American statistical association 63.324, pp. 1379–1389.



_____
#### DOI
Dataset: [![DOI](https://zenodo.org/badge/447151538.svg)](https://zenodo.org/badge/latestdoi/447151538)

Full article: https://doi.org/10.1002/joc.8002



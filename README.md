# SnowTrendsAnalysis

In the following repository are included all the data, plots, codes and operations used for the snow trends analysis. More in detail, the repository contains:
- **PreliminaryData_1**, **PreliminaryData_2**: *.zip* folders with all the preliminary operations on the raw data, as the data quality check and the gap filling procedure
- **FinalData**: folder with the final version data of fresh snow (HN), snow depth (HS), precipitation (P), mean air temperature (TMEAN), north atlantic oscillation index (NAO) and arctic oscillation index (AO). *_S* stays for seasonal data. In addition it is also included the **TopographicData** folder, with *DEM* and *.shp* files of Trentino - South Tyrol region
- **R**: folder with all the R codes used for the analysis. Each step corresponds to a different R script. *00_Data.R* refers to all the preliminary data operations and analyses contained in **PreliminaryData_1** and **PreliminaryData_2**
- **Plots**: folder with the results coming from the analysis of the trends and the attribution validation. Each folder presents the result of each R script in the **R** folder. It also included in the folder **QGIS** the two maps of the HN climatology and study area made using the QGIS software

### Data Analysis
The data analysis was suddived in two parts. The first part was entirely dedicated to the collection and manipulation of fresh snow data in Trentino - South Tyrol. In particular, a monthly aggregation was performed starting from daily data. The dataset was then updated to 2019/2020 season and subjected to a gap filling procedure. At this point, the dataset was analysed and a quality check procedure was conducted.

Instead, the second part was devoted to the real analysis of the data. More in detail, first of all a trend analysis was conducted using statistical tests such as the MannKendall trend test (cit.) and the Sen's slope estimator (cit.). Then, the trend attribution was evaluated, first using a correlation analysis and then an explained variance analysis.

# DRY_NAM
Drought is the major stressor affecting and defining arid and semi-arid regions worldwide. Characterizing drought onset timing, duration, and intensity is important for understanding how drought events impact environmental or agricultural systems. The arid and semi-arid region of the Southwestern United States (SWUS) and Northwestern Mexico (NWMX) is an area that is affected by the North American Monsoon (NAM), which splits the hot-dry season with a period of rainfall that typically begins in mid-summer and lasts through early fall. The bimodal yearly soil water deficit, coupled with multiyear droughts when winter rainfall remains insufficient to refill soil water content capacities under the NAM climate, challenges the generation of standard drought event characterization, yet of critical important agro-environmental impacts and climate change assessment. Here  To address this complexity, we present here a dataset that provides a comprehensive characterization of annual drought facets across the NAM region from 1960 to 2024, based on a new methodological analytical framework applied toof daily soil water content time series. The dataset provides 45 annual drought features, including onset, offset, duration, peak, drying and wetting rates, and severity, for both the pre- and post-monsoon drought periods at 0.1° spatial resolution. The dataset is derived from the daily Soil Water Content, simulated by a generic solar-radiation driven water budget model, which was evaluated against a regional soil water probes network in the region and the standard SPEI-1. Data are stored in 40,825 lat-long pixels (0.1° resolution) CSV files (0.1° resolution), each containingwith 65 rows (years) and 45 columns (representing drought features), enabling detailed spatiotemporal analyses and agro-environmental drought impact assessments. 
ref publication
<!-- badges: start -->

<!-- badges: end -->

## Project files
>  ### "Images"
>Contains images on the code process and drought features. **Figure_3:** "SOLar radiation-driven Drought Index (SOLDI) time series and LOESS smoothing as a function of day of the year for the year 2012 of the >pixel (-96.25,31.12). Black lines show daily SOLDI values; the red line is the corresponding LOESS smoothing curve, and asterisks indicate the detected peaks by the code." ; **FIgure_4:** "Fitted curves over >SOLar radiation-driven Drought Index (SOLDI) and drought facets for the pixel (-96.25,31.12). Solid black line represent the daily SOLDI values;  red curve indicates the fitting of the 1st drought drying phase; >dark blue curve indicates the fitting curve of the 1st drought wetting phase; orange curve indicates the fitting curve of the 2nd drought drying phase; light blue curve indicates the fitting curve of the 2nd >drought wetting phase; Triangles indicate drought onsets; crosses indicate drought offsets; black and grey horizontal dotted lines indicate drought thresholds; dark red, orange and yellow bold solid lines >indicate respectively extreme, moderate and low drought durations; red and orange asterisks indicate 1st and 2nd drought peaks respectively. (a): Wetting.Phase.Duration: (b): Monsoon.Peak; (c): >Dry.Monsoon.Duration; (d): Mod.Monsoon.Duration; (e): Wet.Monsoon.Duration; (f): Monsoon.Peak.DOY."; **Figure_5:** Characterization of a four-year multiyear drought over the 2010-2014 period of the pixel >-104.25,33.62. The black curve represents daily SOLDI values; red curves indicate onset phases; blue curves indicate offset phases; horizontal yellow, orange and red lines represent respectively the low, >moderate and extreme drought durations. Red asterisks are drought peaks; triangles represent onsets; black crosses represent offsets; horizontal black and grey dashed lines are the drought thresholds. ; >**Figure_6:** "Flowchart describing the processing of SOLDI data up to the final table."
>  
>  ### "SWC_example_files"
>  Contains the 3 raw SWC daily values used in the reference article and applicable to the code as examples. **SWC_-101.55,18.12.csv:** mainly unimodal drought dynamics (see Figire_2 (b)) ; >**SWC_-104.25,33.62.csv:** example of a multiyear drought from 2011 to 2014 (see Figure_5) ; **SWC_-96.25,31.12.csv :** mainly bimodal drought dynamics (see Figire_2 (a)).
>  
>  ### "DRY_NAM_2026_03_05.R"
>  The code used to extract the 44 drought features of the DRY_NAM dataset. It contains these main parts :
>  - `1. Initializing the data frame`
>  The final data frame is created
>  - `2. Load the functions`
>  The 2 functions necessary to execute the process. first one 'peakfunction' based on the RamanMP package (not available anymore) to detect the peaks of the drought & and the second one 'process_file' to >extract the drought facets
>  - `3.a. Execute the functions with 'lapply'`
>  To execute the functions in a simple way
>  - `3.b. Execute the functions with 'parallel'`
>  To execute the functions in a parallel way

## DRY_NAM_2026_03_05.R code description
The code extracts 44 yearly drought features (see Table below) from a model-based solar radiation-derived daily drought index (called SOLDI) whether the drought dynamics is unimodal or bimodal. This code has been applied on the North American Monsoon region over 1960-2024. This dataset is called DRY_NAM and can be found here: **dataset link**
More information in the reference data descriptor : **data paper link**

![Texte alternatif](https://github.com/LudwigBoussaroque/DRY_NAM/blob/main/Images/Figure_6_cropped.pdf)

*Figure 1 : Flowchart of the DRY_NAM code process*


## Features names, unit, and description
| Feature                   | Units         | Description                                                                                                 |
| ------------------------- | ------------- | ----------------------------------------------------------------------------------------------------------- |
| First.Drying.Rate         | mm/day        | 1st drought maximum daily soil water loss                                                                   |
| First.Wetting.Rate        | mm/day        | 1st drought maximum daily soil water recovery                                                               |
| First.Peak.Value          | mm            | 1st drought maximum SOLDI value over the year                                                               |
| First.Peak.DOY            | DOY           | 1st drought day of the year when SOLDI reaches its peak value                                               |
| First.Low.D.Onset         | DOY           | When 1st drought period exceeds the SOLDI 25% maximum value                                                 |
| First.Mod.D.Onset         | DOY           | When 1st drought period exceeds the SOLDI 50% maximum value                                                 |
| First.Extreme.D.Onset     | DOY           | When 1st drought period exceeds the SOLDI 75% maximum value                                                 |
| First.Low.D.Offset        | DOY           | When 1st drought period drops below the SOLDI 25% maximum value                                             |
| First.Mod.D.Offset        | DOY           | When 1st drought period drops below the SOLDI 50% maximum value                                             |
| First.Extreme.D.Offset    | DOY           | When 1st drought period drops below the SOLDI 75% maximum value                                             |
| First.Low.D.Duration      | Days          | 1st drought period spanning between Low Drought Onset & Offset                                              |
| First.Mod.D.Duration      | Days          | 1st drought period spanning between Moderate Drought Onset & Offset                                         |
| First.Extreme.D.Duration  | Days          | 1st drought period spanning between Extreme Drought Onset & Offset                                          |
| First.Low.D.S             | mm            | 1st drought SOLDI sum between Low Drought Onset & Offset                                                    |
| First.Mod.D.S             | mm            | 1st drought SOLDI sum between Moderate Drought Onset & Offset                                               |
| First.Extreme.D.S         | mm            | 1st drought SOLDI sum between Extreme Drought Onset & Offset                                                |
| Second.Drying.Rate        | mm/day        | 2nd drought maximum daily soil water loss                                                                   |
| Second.Wetting.Rate       | mm/day        | 2nd drought maximum daily soil water recovery                                                               |
| Second.Peak.Value         | mm            | 2nd drought maximum SOLDI value observed over the year                                                      |
| Second.Peak.DOY           | DOY           | 2nd drought day of the year when SOLDI reaches its peak value                                               |
| Second.Low.D.Onset        | DOY           | When 2nd drought period exceeds the SOLDI 25% maximum value                                                 |
| Second.Mod.D.Onset        | DOY           | When 2nd drought period exceeds the SOLDI 50% maximum value                                                 |
| Second.Extreme.D.Onset    | DOY           | When 2nd drought period exceeds the SOLDI 75% maximum value                                                 |
| Second.Low.D.Offset       | DOY           | When 2nd drought period drops below the SOLDI 25% maximum value                                             |
| Second.Mod.D.Offset       | DOY           | When 2nd drought period drops below the SOLDI 50% maximum value                                             |
| Second.Extreme.D.Offset   | DOY           | When 2nd drought period drops below the SOLDI 75% maximum value                                             |
| Second.Low.D.Duration     | Days          | 2nd drought period spanning between Low Drought Onset & Offset                                              |
| Second.Mod.D.Duration     | Days          | 2nd drought period spanning between Moderate Drought Onset & Offset                                         |
| Second.Extreme.D.Duration | Days          | 2nd drought period spanning between Extreme Drought Onset & Offset                                          |
| Second.Low.D.S            | mm            | 2nd drought SOLDI sum between Low Drought Onset & Offset                                                    |
| Second.Mod.D.S            | mm            | 2nd drought SOLDI sum between Moderate Drought Onset & Offset                                               |
| Second.Extreme.D.S        | mm            | 2nd drought SOLDI sum between Extreme Drought Onset & Offset                                                |
| Wetting.Phase.Duration    | Days          | Duration between First.Peak.DOY and Peak.Monsoon.DOY                                                        |
| Dry.Monsoon.Duration      | Days          | Duration between First.Extreme.D.Offset (or First.Peak.DOY) and Second.Extreme.D.Onset (or Second.Peak.DOY) |
| Mod.Monsoon.Duration      | Days          | Duration between First.Mod.D.Offset (or First.Peak.DOY) and Second.Mod.D.Onset (or Second.Peak.DOY)         |
| Wet.Monsoon.Duration      | Days          | Duration between First.Low.D.Offset (or First.Peak.DOY) and Second.Low.D.Onset (or Second.Peak.DOY)         |
| Peak.Monsoon              | mm            | Difference between First.Peak.Value and the SOLDI value at Peak.Monsoon.DOY                                 |
| Peak.Monsoon.DOY          | DOY           | Lowest SOLDI value between the two droughts                                                                 |
| Cumulative.Year           | mm            | Total sum of the SOLDI values over the year                                                                 |
| Unimodal                  | TRUE or FALSE | Whether or not the drought of the year is unimodal                                                          |


## Pre requisites
The code has been created with:
- **R** (version 4.4.1 (2024-06-14))
- **R Packages** :
  - `terra` (1.8-54)
  - `lubridate` (1.9.4)
  - `tidyr` (1.3.1)
  - `phenofit` (0.3.10)
  - `parallel`
  - 
To install packages:
```r
install.packages(c("terra", "lubridate", "tidyr", "phenofit", "parallel"))
```

## License
This code is under the CC0 1.0 Universal license. See [LICENSE](LICENSE) for more details.


## Example
- Define the working directory, where the SWC files are located
- Select the SWC file(s) from which to extract the features
- Define the output directory, where the final output DF_DRY_NAM files will will be located
```r
#setwd("C:/Users/UM/DRY_NAM/SWC")
setwd("C:/Users/bouss/Desktop/PhD/Data_Monsoon")

# Define the output directory
output_dir <- "C:/Users/bouss/Desktop/Data_Monsoon/"
  
# Load the files to process
#files_to_process <-  c("SWC_-96.25,31.12.csv","SWC_-104.25,33.62.csv","SWC_-101.55,18.12.csv")
```
  
- Execute the `1. Initializing the data frame` command lines to create the table.
- Execute the `2. Load the functions` commmand lines to load the functions.
- Then decide whether to execute `3.a. Execute the functions with 'lapply'` or `3.b. Execute the functions with 'parallel'`.


## Reference data descriptor
**Data paper link**

## Contact
ludwig.boussaroque@umontpellier.fr



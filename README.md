# DRY_NAM
Drought is the major stressor affecting and defining arid and semi-arid regions worldwide. Characterizing drought onset timing, duration, and intensity is important for understanding how drought events impact environmental or agricultural systems. The arid and semi-arid region of the Southwestern United States (SWUS) and Northwestern Mexico (NWMX) is an area that is affected by the North American Monsoon (NAM), which splits the hot-dry season with a period of rainfall that typically begins in mid-summer and lasts through early fall. The bimodal yearly soil water deficit, coupled with multiyear droughts when winter rainfall remains insufficient to refill soil water content capacities under the NAM climate, challenges the generation of standard drought event characterization, yet of critical important agro-environmental impacts and climate change assessment. Here  To address this complexity, we present here a dataset that provides a comprehensive characterization of annual drought facets across the NAM region from 1960 to 2024, based on a new methodological analytical framework applied toof daily soil water content time series. The dataset provides 45 annual drought features, including onset, offset, duration, peak, drying and wetting rates, and severity, for both the pre- and post-monsoon drought periods at 0.1° spatial resolution. The dataset is derived from the daily Soil Water Content, simulated by a generic solar-radiation driven water budget model, which was evaluated against a regional soil water probes network in the region and the standard SPEI-1. Data are stored in 40,825 lat-long pixels (0.1° resolution) CSV files (0.1° resolution), each containingwith 65 rows (years) and 45 columns (representing drought features), enabling detailed spatiotemporal analyses and agro-environmental drought impact assessments. 
ref publication
<!-- badges: start -->

<!-- badges: end -->

## Project structure
### Data
Contains the 3 raw SWC daily values used in the reference article and applicable to the code as examples

### DRY_NAM.R
The code contains these main parts :
- `1. Initializing the data frame`
The final data frame is created
- `2. Load the functions`
The 2 functions necessary to execute the process. first one 'peakfunction' based on the RamanMP package (not available anymore) to detect the peaks of the drought & and the second one 'process_file' to extract the drought facets
- `3.a. Execute the functions with 'lapply'`
To execute the functions in a simple way
- `3.b. Execute the functions with 'parallel'`
To execute the functions in a parallel way

## DRY_NAM.R code description
The code attached (file name.R) extracts 45 yearly drought facets (see Table below) from a model-based solar radiation-derived daily drought index (called SOLDI) whether the drought isuni modal or bimodal. This code has been applied on the North American Monsoon region over 1960-2024 where the data can be find here : Zenodo link 
Reference publication : data paper link

![Texte alternatif](https://github.com/LudwigBoussaroque/DRY_NAM/blob/main/Images/Figure_6_cropped.pdf)

*Figure 1 : Flowchart of the DRY_NAM code process*


## Facets description
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
- **Paquets R** :
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
This code is ...
### License
Ce projet est sous licence **MIT**. Voir le fichier [LICENSE](LICENSE) pour plus de détails.


## Example
- Define the working directory, where the SWC files are located
- Select the SWC file(s) from which to extract facets
- Define the output directory, where the final output DF_DRY_NAM files will will be located
```r
#setwd("C:/Users/UM/DRY_NAM/SWC")
setwd("C:/Users/bouss/Desktop/PhD/Data_Monsoon")

# Define the output directory
output_dir <- "C:/Users/bouss/Desktop/"
  
# Load the files to process
#files_to_process <-  c("")
```
  
- Execute the `1. Initializing the data frame`
- Execute the `2. Load the functions`
- Then execute wether `3.a. Execute the functions with 'lapply'` or `3.b. Execute the functions with 'parallel'` wether 


## Citation

## Contact


![R Version](https://img.shields.io/badge/R-4.2-blue)
![Licence](https://img.shields.io/badge/Licence-MIT-green)


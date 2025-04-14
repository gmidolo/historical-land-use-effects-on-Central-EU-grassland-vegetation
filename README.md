This repository contains data and code associated with the manuscript:

> Nineteenth-century land use shapes the current occurrence of some plant species, but weakly affects the richness and total composition of Central European grasslands

published in *Landscape Ecology* (https://doi.org/10.1007/s10980-024-02016-6)

## Authors

* Midolo, Gabriele<sup>1</sup>\*
* Skokanová; Hana<sup>2</sup>
* Clark, Adam Thomas<sup>3</sup>
* Vymazalová, Marie<sup>2</sup>
* Chytrý, Milan<sup>4</sup>
* Dullinger, Stefan<sup>5</sup>
* Essl, Franz<sup>6</sup>
* Šibík, Jozef<sup>7</sup>
* Keil, Petr<sup>1</sup>

## Affiliations

1: Department of Spatial Sciences, Faculty of Environmental Sciences, Czech University of Life Sciences Prague, Praha-Suchdol, Czech Republic

2: Silva Tarouca Research Institute for Landscape and Ornamental Gardening, Department of Landscape Ecology, Brno, Czech Republic

3: Department of Biology, University of Graz, Graz, Austria

4: Department of Botany and Zoology, Faculty of Science, Masaryk University, Brno, Czech Republic

5: Division of Biodiversity Dynamics and Conservation, Department of Botany, University Vienna, Vienna, Austria

6: Division of BioInvasions, Global Change & Macroecology, Department of Botany and Biodiversity Research, University of Vienna, Vienna, Austria

7: Plant Science and Biodiversity Center, Slovak Academy of Sciences, Bratislava, Slovakia

\* **Correspondence:** Gabriele Midolo; e-mail: midolo@fzp.czu.cz; ORCID: https://orcid.org/0000-0003-1316-2546

## Repository Contents

The repository contains the following material:

* **IndVal.all.habitats.csv:** The results of the IndVal statistics (De Cáceres et al. 2010) for 1,498 species for the historical land use categories calculated across the entire dataset.

* **IndVal.separate.habitats.csv:** The results of the IndVal statistics for 1,498 species for the historical land use categories calculated for each habitat type (dry grasslands, mesic grasslands, wet grasslands) separately.

* **ecological.and.disturbance.values.csv:** The original Ellenberg-type and disturbance indicator values (data obtained from Tichý et al. 2023 and Midolo et al. 2023; data can be accessed at the [FloraVeg.eu](http://www.floraveg.eu/) website), and the varimax-rotated components (‘RC’) used in the analysis for 1,461 species.

* **Rcode folder:**

    * **Rdata.RDS:** An R data file (.RDS) containing the data to reproduce the analyses in R.  It contains the following objects:

        * `indic.val`: A data.frame with 1,461 observations and 15 variables; it contains the same data as 'ecological.and.disturbance.values.csv'.

        * `plot.data`: A list of 3 data.frames, each containing vegetation plot data from the European Vegetation Archive (EVA; Chytrý et al. 2016) for dry (R1), mesic (R2) and wet (R3) grasslands. Data includes plot IDs, historical land use, plot size (m2), bioclimatic variables (‘bio’; Karger et al. 2017), soil pH (Hengl et al. 2017), and country (Austria or Czechia/Slovakia). Geographic coordinates are not shared.

        * `species.data`: A list of 3 matrices. Each is a site x species matrix containing the relative abundance of species (columns) and plot sites (rows) for dry (R1), mesic (R2) and wet (R3) grasslands elaborated from EVA.

    * **Scripts (\*.R):** R scripts for reproducing main analyses on species richness, species composition, and species indicator analyses. R scripts are also rendered in .html with R Markdown.

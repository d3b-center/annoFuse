## About `shinyFuse`

`ShinyFuse` is an application for exploration of fusion output derived from [annoFuse](https://github.com/d3b-center/annoFuse).

This app has two main sections, each with specific functionality, described here in brief.

#### About FusionExplorer

FusionExplorer is an interactive way to explore the putative oncogenic fusion results from annoFuse. 

**Required input:** `PutativeOncogenicFusion.tsv`

**Features:** 

_Filtering_: Select a column and enter one or more attributes on which to filter the rows. <br>
_Data export_: Export the filtered table using the download button above.<br>
_Fusion visualization_: To view the gene fusion, exons, protein domains, and breakpoints:

 1. Load the `exons` and `pfam` domains by clicking the buttons on the upper right side. 
 1. Select any row in the data table to visualize that specific fusion. Plots for the 5' gene, 3' gene, and both genes will be generated to the right of the table.

_Plot export_: To save the fusion plot, select the green download arrow below the plot.

#### About FusionSummary

**Features:** 

_Recurrent fusions_: To analyze recurrent fusions within a cohort:

 1. To the left of the table, select a `Grouping Column`. This will be the column used for analysis of recurrent fusions. For example, if your dataset has multiple histologies and you are interested in recurrent fusions by histology, this value would be histology.
 1. Select the number of top fusions to display in the plot.
 1. Select a `Counting column`. This is usually a patient-level identifier.

_Fusion visualization_: To view the recurrent fusions results, click on the `FusionSummary` tab, where plots for recurrent fusions and recurrently-fused genes will be generated.<br>
_Plot export_: To save the fusion plot, select the green download arrow below the plot.

For a more interactive way to learn how to use `shinyFuse`, please take the tour, which can be activated by clicking on the question mark in the upper right corner, and test out the demo dataset by clicking on `Load demo data` in the upper left corner. 
Thanks for visiting!

## News

[Date of Launch]. Welcome to ShinyFuse, an application for exploration of fusion output derived from [annoFuse](https://github.com/d3b-center/annoFuse)! 

## Contact

Please file application issues, suggestions, or bug reports [here](https://github.com/d3b-center/annoFuse/issues).

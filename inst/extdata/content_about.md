## About `shinyFuse`

`ShinyFuse` is an application for exploration of fusion output derived from [annoFuse](https://github.com/d3b-center/annoFuse).

This app has two main sections, each with specific functionality, described here in brief.

#### About FusionExplorer

FusionExplorer is an interactive way to explore the putative oncogenic fusion results from annoFuse. 

**Required input:** `PutativeDriverAnnoFuse.tsv`

**Features:** 

_Filtering_: On the left, select columns to display:

1. If filtering a non-numeric column, select a level from the dropdown on which to filter rows.
2. If filtering a numeric column, use the sliding scale or enter numeric logic according to the following (Eg: to filter for a fusion `SpanningFragCount` greater than 10, write `10 ... max` (number, space, ellipsis, space, number) in the filter box, where `max` is the highest number on the sliding scale). <br>

_Data export_: Export the filtered table using the download button above.<br>
_Fusion visualization_: To view the gene fusion, exons, protein domains, and breakpoints:

 1. Load the `exons` and `pfam` domains by clicking the buttons on the upper right side. 
 2. Select any row in the data table to visualize a specific fusion. Plots for the 5' gene, 3' gene, and both genes will be generated to the right of the table. Plots can be visualized by breakpoint (row), sample, or the entire cohort of samples.

_Plot export_: To save the fusion plot, select the green download arrow below the plot.

#### About FusionSummary

**Features:** 

_Recurrent fusions_: To analyze recurrent fusions within a cohort:

 1. To the left of the table, select a `Grouping Column`. This will be the column used for analysis of recurrent fusions. For example, if your dataset has multiple histologies and you are interested in recurrent fusions by histology, this value would be histology.
 2. Select the number of top fusions to display in the plot.
 3. Select a `Counting column`. This is usually a patient-level identifier.
 4. Select additional plotting filters below:
 
 | Filter                  | Description                                                                                                                                                                                                                     | Options                                                                                                                                           |
|-------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| Fusion type             | Consequence on reading frame of fusion                                                                                                                                                                                          | in-frame, frameshift, other                                                                                                                       |
| Caller                  | Algorithm used to call fusions                                                                                                                                                                                                  | Arriba, STARFusion                                                                                                                                |
| Confidence              | Confidence for fusion call (specific to Arriba)                                                                                                                                                                                 | high, medium, low, NA                                                                                                                             |
| Breakpoint Location     | Location of fusion breakpoint                                                                                                                                                                                                   | Genic = within gene body of both genes, Intragenic = within one gene, Intergenic = one breakpoint is within a gene and one is between two genes |
| Spanning Fragment Count | [number of RNA-Seq fragments that encompass the fusion junction such that   one read of the pair aligns to a different gene than the other paired-end read of that fragment](https://github.com/STAR-Fusion/STAR-Fusion/wiki) | numeric                                                                                                                                           |
| Junction Read Count     | [number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction](https://github.com/STAR-Fusion/STAR-Fusion/wiki)                                           | numeric                                                                                                                                           |
| Caller Count            | Number of algorithms from which fusion was called                                                                                                                                                                               | max N callers                                                                                                                                     |


_Fusion visualization_: To view the recurrent fusions results, click on the `FusionSummary` tab, where plots for recurrent fusions and recurrently-fused genes will be generated.<br>
_Plot export_: To save the fusion plot, input the desired dimensions on the left and/or select the green download arrow below the plot.

For a more interactive way to learn how to use `shinyFuse`, please take the tour, which can be activated by clicking on the question mark in the upper right corner, and test out the demo dataset by clicking on `Load demo data` in the upper left corner. 
Thanks for visiting!

## News

[Date of Launch]. Welcome to ShinyFuse, an application for exploration of fusion output derived from [annoFuse](https://github.com/d3b-center/annoFuse)! 

## Contact

Please file application issues, suggestions, or bug reports [here](https://github.com/d3b-center/annoFuse/issues).

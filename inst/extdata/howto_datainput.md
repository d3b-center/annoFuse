- `shinyFuse` accepts **tab-separated text files**, as commonly produced by the workflow of `annoFuse`
- The file should contain a **header row**, with minimal column names: `Sample`, `LeftBreakpoint`, `RightBreakpoint`, `FusionName`, `Gene1A`, and `Gene1B`, with rows formatted as shown in the example below.
- If your data in stored in Excel, export the file as tab-separated format for compatibility with `shinyFuse`

Below is a screenshot of for a minimal file which will be correctly loaded by `shinyFuse`.
To take full advantage of all features, we recommend running fusion analysis, annotation, followed by `annoFuse` as noted [here](https://github.com/d3b-center/annoFuse#install-package).
For a fully annotated input file example, see the [`PutativeDriverAnnoFuse.tsv`](https://github.com/d3b-center/annoFuse/blob/master/inst/extdata/PutativeDriverAnnoFuse.tsv) provided bundled with the package.

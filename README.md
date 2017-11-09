Extract barcode sequences and correct sequencing errors
````
bxcheck trim -l known-barcodes.txt input.fq.gz -o output
````

Plot some stats
````
plot-bxcheck.py output.txt -d fq-plots/
````

Run stats on mapped 10x BAM, requires BX tags
````
bxcheck stats -l known-barcodes.txt mapped.bam > stats.txt
````

Plot the stats
````
plot-bxcheck.py stats.txt -d bam-plots/
# or
plot-bxcheck.py output.txt stats.txt -d fq-and-bam-plots/
````

Plot summary stats highlighting differences between multiple samples
````
cat summary.conf 
  fSpaAur1    /path/to/fSpaAur1/trim.txt
  fSpaAur1    /path/to/fSpaAur1/bxcheck.txt
  fCotGob3    /path/to/fCotGob3/trim.txt
  fCotGob3    /path/to/fCotGob3/bxcheck.txt
  mOncTor1    /path/to/mOncTor1/trim.txt
  mOncTor1    /path/to/mOncTor1/bxcheck.txt
  fErpCal1    /path/to/fErpCal1/trim.txt
  fErpCal1    /path/to/fErpCal1/bxcheck.txt
  fAnaTes1    /path/to/fAnaTes1/trim.txt
  fAnaTes1    /path/to/fAnaTes1/bxcheck.txt
  fMasArm1    /path/to/fMasArm1/trim.txt
  fMasArm1    /path/to/fMasArm1/bxcheck.txt
plot-summary.py -i summary.conf summary.html
````

An example of [plot-bxcheck output](https://pd3.github.io/bxcheck/example-plot-bxcheck.html)
and [plot-summary output](https://pd3.github.io/bxcheck/example-plot-summary.html).


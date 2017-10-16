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

An example of the output is shown [here](https://pd3.github.io/bxcheck/example-plot-bxcheck.html).

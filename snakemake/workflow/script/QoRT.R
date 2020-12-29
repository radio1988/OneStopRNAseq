library("DESeq2")
library("QoRTs")
#Read in the QC data:
print(getwd())
res <- read.qc.results.data("bam_qc/QoRTs/",
                            decoder.files = "meta/decoder.txt",
                            calc.DESeq2 = TRUE, calc.edgeR = FALSE);
out <- "bam_qc/QoRTs_MultiPlot/"

if (!dir.exists(out))
{
  dir.create(out, recursive = TRUE)
}

makeMultiPlot.all(res,
                  outfile.dir = out,
                  plot.device.name = "pdf");

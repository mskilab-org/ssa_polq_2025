library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gUtils)
library(GxG)

seq.similarity <- function(x){
  return(1 - (x / ((pad * 2) + stride)))
}

#### LARGE AMPLICON ###
#seq1 <- "TTACCTCTCTAGTCTGTGCTAGCTCttccagccccctgtcatggcatcttccaggggtccgagagctcagctagtcttcttcctccaacccgggcccctatgtccacttcaggacagcatgtttgctgcctccagggatcctgtgtccccgagctgggaccaccttatattcccagggccggttaatgtggctctggttctgggtacttttatctgtcccctccaccccacagtggggccactagggacaggattggtgacagaaaagccccatccttaggcctcctccttcctagtctcctgatattgggtctaacccccacctcctgttaggcagattccttatctggtgacacacccccatttcctggagccatctctctccttgccagaacctctaaggtttgcttacgatggagccagagaggatcctgggagggagagcttggcagggggtgggagggaagggggggatgcgtgacctGCCCGGTTCTCAGTGGCCA"

seq1_coords <- GRanges(seqnames = "chr19",
        ranges = IRanges(start = 55626871, end = 55627373),
        strand = "-")

real_coords1 <- GRanges(seqnames = "AMPLICON",
        ranges = IRanges(start = -251, end = 251),
        strand = "*")
         
gr1 <- GRanges(seq1_coords) %>% gr.chr()
strand(gr1) <- "-"
real_coords1$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, seq1_coords)

hm1 <- homeology(real_coords1, 
pad = pad, 
stride = stride, 
verbose = T)

hm1 <- hm1$transform(seq.similarity)

plot(hm1$gtrack(), windows = GRanges("AMPLICON:-100-100")) %>% skitools::ppdf("LARGE.pdf")

custom_colors <- rev(scico::scico(n = 19, palette = "batlowK"))
hm1$gtrack(col = custom_colors) %>% plot(windows = GRanges("AMPLICON:-100-100")) %>% skitools::ppdf("BATLOW_LARGE.pdf")

### SMALL AMPLICON ###

seq2_coords <- GRanges(seqnames = "chr19",
        ranges = IRanges(start = 55626947, end = 55627192),
        strand = "-")

real_coords2 <- GRanges(seqnames = "AMPLICON",
        ranges = IRanges(start = -69, end = 176),
        strand = "*")

seqinfo(seq2_coords) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)["chr19"]
gr2 <- GRanges(seq2_coords) %>% gr.chr()
strand(gr2) <- "-"
real_coords2$seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, seq2_coords)


hm2 <- homeology(real_coords2, 
pad = pad, 
stride = stride, 
verbose = T)

hm2 <- hm2$transform(seq.similarity)

plot(hm2$gtrack(), windows = GRanges("AMPLICON:-69-69")) %>% skitools::ppdf("SMALL.pdf")

custom_colors <- rev(scico::scico(n = 19, palette = "batlowK"))
hm2$gtrack(col = custom_colors) %>% plot(windows = GRanges("AMPLICON:-69-69")) %>% skitools::ppdf("BATLOW_SMALL.pdf")


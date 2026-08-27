# SuSiE 95% credible sets for the 8 GEUVADIS genes:
# is the published lead (or a perfect proxy) inside a CS, and how large is it?
# Run from the project root.
suppressPackageStartupMessages(library(susieR))
D <- "data/geuvadis"; RES <- "results/bench_full/26_geuvadis"
GENES <- data.frame(
  tag  = c("ENSG00000164308","ENSG00000197728","ENSG00000166750","ENSG00000203875",
           "ENSG00000198468","ENSG00000124587","ENSG00000230658","ENSG00000174652"),
  full = c("ENSG00000164308.12","ENSG00000197728.5","ENSG00000166750.4","ENSG00000203875.4",
           "ENSG00000198468.2","ENSG00000124587.9","ENSG00000230658.1","ENSG00000174652.12"),
  gene = c("ERAP2","RPS26","SLFN5","SNHG5","FLVCR1-AS1","PEX6-region","TRA2A-AS","ZNF266"),
  chr  = c("5","12","17","6","1","6","7","19"),
  lead_pos = c(96252589L,56401085L,33571546L,86387888L,213049214L,42944850L,23143113L,9544276L),
  stringsAsFactors = FALSE)
hdr <- strsplit(readLines(gzfile(file.path(D,"GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")), n=1L),"\t")[[1]]
gt2d <- function(v){d<-integer(length(v)); d[v%in%c("0|1","1|0")]<-1L; d[v%in%c("1|1")]<-2L; d}
rows <- list()
for (i in seq_len(nrow(GENES))) {
  g <- GENES[i,]
  ord <- system(sprintf("bcftools query -l -S %s --force-samples %s 2>/dev/null",
    shQuote(file.path(D,"geuvadis_eur_overlap.txt")),
    shQuote(sprintf("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", g$chr))), intern=TRUE)
  sp <- strsplit(readLines(file.path(D,"loci",paste0(g$tag,".gt"))),"\t",fixed=TRUE)
  pos <- as.integer(vapply(sp,`[`,"",1L))
  G <- vapply(sp, function(r) gt2d(r[-1L]), integer(length(ord)))
  maf <- colMeans(G)/2; maf <- pmin(maf,1-maf)
  keep <- maf >= 0.05 & apply(G,2,sd) > 0
  G <- G[,keep,drop=FALSE]; pos <- pos[keep]
  con <- gzfile(file.path(D,"GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")); open(con)
  y_row <- NULL
  while (length(l <- readLines(con, n=4000L))) {
    hit <- grep(paste0("^",g$full,"\t"), l)
    if (length(hit)) { y_row <- strsplit(l[hit[1]],"\t")[[1]]; break } }
  close(con)
  expr <- as.numeric(y_row[-(1:4)]); names(expr) <- hdr[-(1:4)]
  y <- as.numeric(scale(as.numeric(expr[ord])))
  X <- scale(G); X[!is.finite(X)] <- 0
  fit <- susie(X, y, L = 10)
  cs <- susie_get_cs(fit, X = X, coverage = 0.95)$cs
  j_lead <- match(g$lead_pos, pos)
  g_lead <- G[, j_lead]
  in_cs <- FALSE; cs_size <- NA_integer_
  for (c1 in cs) {
    hit <- j_lead %in% c1 ||
           any(vapply(c1, function(s) cor(G[,s], g_lead)^2 >= 0.8, logical(1)))
    if (hit) { in_cs <- TRUE; cs_size <- length(c1); break }
  }
  rows[[i]] <- data.frame(gene = g$gene, n_cs = length(cs),
                          lead_in_cs = in_cs, cs_size = cs_size)
  cat(sprintf("%-12s n_CS=%d  lead(or r2>=0.8 proxy) in CS: %s  size=%s\n",
              g$gene, length(cs), in_cs, cs_size))
}
out <- do.call(rbind, rows)
write.csv(out, file.path(RES, "susie_cs_summary.csv"), row.names = FALSE)
cat("saved susie_cs_summary.csv\n")

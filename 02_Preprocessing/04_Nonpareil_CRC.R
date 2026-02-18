#############################################
## NONPAREIL ANALYSIS FOR YOUR COHORTS
## USING .npo FILES INSIDE PROJECT FOLDERS
#############################################

set.seed(45)
suppressPackageStartupMessages({
  library(Nonpareil)
  library(dplyr)
  library(readr)
  library(tibble)
})

## ------------------------------------------------------------
## 1. MAP PROJECT FOLDER TO COHORT LABEL
## ------------------------------------------------------------

proj_to_cohort <- c(
  "PRJDB4176"    = "Cohort_1",
  "PRJEB27928"   = "Cohort_2",
  "PRJEB6070"    = "Cohort_3",
  "PRJEB72525"   = "Cohort_4",
  "PRJEB72526"   = "Cohort_5",
  "PRJEB7774"    = "Cohort_6",
  "PRJNA1167935" = "Cohort_7",
  "PRJNA389927"  = "Cohort_8"
)

## The folders you have:
projects <- names(proj_to_cohort)

## Base directory containing the project folders
base_dir <- "."

## ------------------------------------------------------------
## 2. COLORS FOR COHORTS (distinct & consistent)
## ------------------------------------------------------------
cohort_colors <- c(
  "Cohort_1" = "#1f77b4",
  "Cohort_2" = "#ff7f0e",
  "Cohort_3" = "#2ca02c",
  "Cohort_4" = "#d62728",
  "Cohort_5" = "#9467bd",
  "Cohort_6" = "#8c564b",
  "Cohort_7" = "#e377c2",
  "Cohort_8" = "#7f7f7f"
)

## ------------------------------------------------------------
## 3. BUILD SAMPLE LIST
## ------------------------------------------------------------
samples <- bind_rows(lapply(projects, function(pid) {
  
  folder <- file.path(base_dir, paste0("Nonpareil_", pid))
  if (!dir.exists(folder)) return(NULL)
  
  files <- list.files(folder, pattern = "\\.npo$", full.names = TRUE)
  if (!length(files)) return(NULL)
  
  tibble(
    File    = files,
    Name    = basename(files),
    Cohort  = proj_to_cohort[[pid]],
    Color   = cohort_colors[[proj_to_cohort[[pid]]]]
  )
}))

if (nrow(samples) == 0) {
  stop("No .npo files found in Nonpareil_* project folders.")
}

write.table(samples, "nonpareil_samples.tsv", sep="\t",
            row.names=FALSE, quote=FALSE)


## ------------------------------------------------------------
## 4. GLOBAL NONPAREIL PLOT (all cohorts)
## ------------------------------------------------------------
png("Nonpareil_All_Cohorts.png", width=4, height=4, units="in", res=600)

nps <- Nonpareil.set(
  files     = samples$File,
  col       = samples$Color,
  labels    = samples$Cohort,
  plot      = TRUE,
  plot.opts = list(
    plot.observed   = TRUE,
    plot.model      = TRUE,
    plot.dispersion = FALSE,
    legend          = FALSE,
    lwd             = 1.5,
    main            = ""
  )
)

axis(1, lwd=2.5, lwd.ticks=1.2)
axis(2, lwd=2.5, lwd.ticks=1.2)
box(lwd=2.5)

legend("bottomright",
       legend = names(cohort_colors),
       col    = cohort_colors,
       lwd    = 2,
       cex    = 0.5,
       bg     = "white", bty="n")

dev.off()

## ------------------------------------------------------------
## 5. SUMMARY TABLE FOR ALL SAMPLES
## ------------------------------------------------------------
sum_all <- as.data.frame(summary(nps)) %>%
  mutate(
    Name   = samples$Name,
    File   = samples$File,
    Cohort = samples$Cohort
  ) %>%
  select(Cohort, Name, File, kappa, C, LR, LRstar, modelR, diversity)

write_csv(sum_all, "nonpareil_summary_all_cohorts.csv")

## ------------------------------------------------------------
## 6. PER-COHORT PLOTS AND SUMMARY
## ------------------------------------------------------------

dir.create("nonpareil_per_cohort", showWarnings = FALSE)

for (ct in names(cohort_colors)) {
  
  sub <- samples %>% filter(Cohort == ct)
  if (nrow(sub) == 0) next
  
  out_png <- file.path("nonpareil_per_cohort", paste0(ct, ".png"))
  
  png(out_png, width=4, height=4, units="in", res=600)
  
  Nonpareil.set(
    files     = sub$File,
    col       = sub$Color,
    labels    = sub$Name,
    plot      = TRUE,
    plot.opts = list(
      plot.observed   = TRUE,
      plot.model      = TRUE,
      plot.dispersion = FALSE,
      legend          = FALSE,
      lwd             = 1.4,
      main            = ct
    )
  )
  
  axis(1, lwd=2.5, lwd.ticks=1.2)
  axis(2, lwd=2.5, lwd.ticks=1.2)
  box(lwd=2.5)
  
  dev.off()
  
  ## Summary for this cohort
  s_proj <- as.data.frame(summary(Nonpareil.set(sub$File, plot=FALSE))) %>%
    mutate(
      Name   = sub$Name,
      File   = sub$File,
      Cohort = ct
    ) %>%
    select(Cohort, Name, File, kappa, C, LR, LRstar, modelR, diversity)
  
  write_csv(s_proj, file.path("nonpareil_per_cohort", paste0("summary_", ct, ".csv")))
}


###############################################################################
## 7. AUTOMATIC SUMMARY STATISTICS (GLOBAL + COHORT) — CLEAN VERSION
##     mean ± sd and median ONLY
##############################################################################

library(dplyr)
library(readr)

message("Computing summary statistics…")

# Reload the main summary file
df <- read_csv("nonpareil_summary_all_cohorts.csv")

numeric_cols <- c("kappa","C","LR","LRstar","modelR","diversity")

# =======================================================
# GLOBAL SUMMARY: mean ± sd and median
# =======================================================
global_stats <- df %>%
  summarise(
    n_samples = n(),
    across(all_of(numeric_cols),
           list(
             mean  = ~mean(.x, na.rm=TRUE),
             sd    = ~sd(.x, na.rm=TRUE),
             median = ~median(.x, na.rm=TRUE)
           ),
           .names="{.col}_{.fn}")
  )

write_csv(global_stats, "nonpareil_global_summary_stats.csv")


# =======================================================
# PER-COHORT SUMMARY: mean ± sd and median
# =======================================================
cohort_stats <- df %>%
  group_by(Cohort) %>%
  summarise(
    n_samples = n(),
    across(all_of(numeric_cols),
           list(
             mean  = ~mean(.x, na.rm=TRUE),
             sd    = ~sd(.x, na.rm=TRUE),
             median = ~median(.x, na.rm=TRUE)
           ),
           .names="{.col}_{.fn}")
  ) %>%
  ungroup()

write_csv(cohort_stats, "nonpareil_per_cohort_summary_stats.csv")

message("Summary statistics generated:")
message(" - nonpareil_global_summary_stats.csv")
message(" - nonpareil_per_cohort_summary_stats.csv")

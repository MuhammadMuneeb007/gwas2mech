# Step09_Config.py - parameters of Step09_Colocalize_GTEx_SuSiE.py
# Edit here; the command-line flags --qtl-types / --study-regex / --candidate-p / --strong-pp4 /
# --suggestive-pp4 override single values for one run. Plan again (and use --force to redo finished loci)
# after changing anything here.

# ---- which QTL contexts (eQTL Catalogue datasets) ---------------------------------------------------
QTL_TYPES = ["eQTL", "sQTL", "pQTL", "isoQTL", "exonQTL"]            # eQTL, sQTL, pQTL, isoQTL, exonQTL, apaQTL
STUDY_REGEX = r"(?i)gtex|^sun_2018$"
# Download the same contexts with Step00 (the LBF files are sorted + indexed there):
#   envs/pipeline/bin/python Step00_Download_Resources.py --groups eqtl_catalogue qtl_catalogue qtl_catalogue_susie
#   (other studies: add --catalogue-study-regex '<the same regex>')
TISSUE_REGEX = ""                  # optional, e.g. r"brain|substantia" (tissue_label / sample_group)
DATASET_IDS = []                   # optional explicit list, e.g. ["QTD000116", "QTD000136"]; overrides STUDY_REGEX

# ---- discovery in the dense QTL summary statistics ---------------------------------------------------
CANDIDATE_P = 1e-5                 # molecular traits with P <= this in the GWAS locus window are tested
MAX_CANDIDATES_PER_CONTEXT_LOCUS = 100

# ---- coloc.bf_bf priors and thresholds ------------------------------------------------------------------
P1 = 1e-4                          # prior: variant associated with the GWAS trait
P2 = 1e-4                          # prior: variant associated with the molecular trait
P12 = 5e-6                         # prior: variant shared (default classification)
P12_LOW = 1e-6                     # prior sensitivity (PRIOR_ROBUST_STRONG = still strong at this p12)
P12_HIGH = 1e-5
STRONG_PP4 = 0.80                  # PP.H4 >= this: STRONG_SHARED_SIGNAL
SUGGESTIVE_PP4 = 0.50              # PP.H4 >= this: SUGGESTIVE_SHARED_SIGNAL
OVERLAP_MIN = 0.5                  # coloc.bf_bf overlap.min
TRIM_BY_POSTERIOR = True           # coloc.bf_bf trim_by_posterior
SAVE_VARIANT_PP_MIN_PP4 = 0.50     # write per-variant PP.H4 for traits at or above this PP.H4

# ---- SLURM (None = take it from the #SBATCH header of HPC_CPU.sh; command-line flags win) --------------
SLURM = {
    "cpus-per-task": None,
    "mem": None,                   # e.g. "32G"
    "time": None,                  # e.g. "24:00:00"
    "partition": None,
    "account": None,
    "qos": None,
}

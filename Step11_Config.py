# Step11_Config.py - parameters of Step11_Custom_QTL_Colocalization.py (external molecular QTL, coloc.abf)
# Edit here; --qtl-types / --candidate-p override single values for one run. Plan again (and use --force to
# redo finished loci) after changing anything here.

# ---- resources --------------------------------------------------------------------------------------
RESOURCE_ROOT = "resources/molecular_qtl"            # scanned by --scan (Step00 --groups molecular_qtl meqtl)
REGISTRY = "config/Step11_Custom_QTL_Registry.tsv"   # written by --scan; enable rows by hand after checking them
QTL_TYPES = ["sQTL", "pQTL", "haQTL", "caQTL", "mQTL"]

# ---- discovery and variant matching -----------------------------------------------------------------
CANDIDATE_P = 1e-5           # molecular traits with P <= this in the GWAS locus window are tested
MAX_TRAITS = 100             # per dataset and locus (smallest P first)
MIN_SHARED = 50              # shared GWAS / QTL variants needed for coloc.abf
MIN_OVERLAP = 0.5            # shared variants / variants of the smaller side
MAX_REGION_ROWS = 300000     # safety limit on rows read per file and locus

# ---- coloc.abf priors and thresholds (single causal variant per trait and locus) ---------------------
P1 = 1e-4
P2 = 1e-4
P12 = 5e-6                   # default classification
P12_LOW = 1e-6               # PRIOR_ROBUST_STRONG = still strong at this p12
P12_HIGH = 1e-5
STRONG_PP4 = 0.80
SUGGESTIVE_PP4 = 0.50
R_TIMEOUT = 3600             # seconds for the R coloc.abf run of one locus

# ---- SLURM (None = take it from the #SBATCH header of HPC_CPU.sh; command-line flags win) --------------
SLURM = {
    "cpus-per-task": None,
    "mem": None,
    "time": None,
    "partition": None,
    "account": None,
    "qos": None,
}

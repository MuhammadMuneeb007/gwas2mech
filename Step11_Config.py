# Step11_Config.py - parameters of Step11_Custom_QTL_Colocalization.py (external molecular QTL)
# Edit here; --qtl-types / --candidate-p override single values for one run. After a change: --scan (resources),
# plan again, and --force to redo finished loci (parameters that change results also invalidate the resume hash).

# ---- resources (found through the Step00 manifests; nothing is searched recursively) --------------------
RESOURCE_ROOT = "resources/molecular_qtl"
REGISTRY = "config/Step11_Custom_QTL_Registry.tsv"     # resource table written by --scan (one row per dataset)
FILES_TABLE = "config/Step11_Custom_QTL_Files.tsv"     # one row per file (manifest md5 / size, present, indexed)

# ADSP FunGen xQTL (NIAGADS NG00184, Step00 --groups molecular_qtl / meqtl):
#   <RESOURCE_ROOT>/<folder>/adsp_fungen_brain/ADSP_FunGen_xQTL.v1.<QTL>.all.manifest.txt -> all/<manifest path>
ADSP_MODALITIES = {"sQTL": "sqtl", "pQTL": "pqtl", "haQTL": "haqtl", "mQTL": "mqtl"}
ADSP_PROVIDER_DIR = "adsp_fungen_brain"
NIAGADS_MANIFEST = "resources/molecular_qtl/manifests/NG00184.v1.csv"   # fallback file listing

# Files in eQTL Catalogue summary-statistics format (named header; sample_size from the catalogue tables)
EQTL_CATALOGUE_FORMAT_DATASETS = [
    {"dataset_id": "QTD100018", "modality": "caQTL", "source_study": "Kumasaka_2019",
     "provider": "eQTL Catalogue re-analysis (Zenodo 13848268)",
     "path": "resources/molecular_qtl/caqtl/kumasaka_grch38/QTD100018.all.tsv.gz",
     "trait_metadata": "resources/molecular_qtl/caqtl/kumasaka_grch38/QTD100018_peak_metadata.tsv.gz"},
    {"dataset_id": "QTD000584", "modality": "pQTL", "source_study": "Sun_2018",
     "provider": "eQTL Catalogue (INTERVAL plasma pQTL)",
     "path": "resources/molecular_qtl/pqtl/interval_sun2018/QTD000584.all.tsv.gz"},
]
EQTL_CATALOGUE_TABLES = [
    "resources/coloc/eqtl_catalogue/tabix_ftp_paths.tsv",
    "resources/coloc/eqtl_catalogue/tabix_ftp_paths_imported.tsv",
    "resources/molecular_qtl/manifests/tabix_ftp_paths.tsv",
    "resources/molecular_qtl/manifests/tabix_ftp_paths_imported.tsv",
]
EQTL_CATALOGUE_SUSIE_DIR = "resources/coloc/eqtl_catalogue/susie_provider"   # provider SuSiE -> Step09 / Step10
EXTERNAL_REGISTRIES = ["resources/molecular_qtl/resource_registry.tsv", "resources/coloc/custom_qtl_registry.tsv"]
REFERENCE_FASTA = "resources/genome/GRCh38/GRCh38.primary_assembly.genome.fa"   # build / allele verification

QTL_TYPES = ["sQTL", "pQTL", "haQTL", "caQTL", "mQTL"]
DISABLED_DATASETS = []        # dataset ids (DATASET_ID of the resource table) to leave out

# ---- documented values the files / metadata do not carry (never guessed) ------------------------------
# QTL sample size: {"DATASET_ID": (N, "source, e.g. paper table / portal page")}
QTL_N = {}
# QTL trait SD (beta / SE without an allele frequency): {"DATASET_ID": (sdY, "source")}
QTL_SDY = {}
# Column roles for headerless files or non-standard names: {"DATASET_ID" or "bed3plus17_qtl": {role: column}}
# roles: chrom start end pos variant ref alt trait gene beta se p freq n an ac z   (headerless columns are c1, c2, ...)
COLUMN_OVERRIDES = {}
ALLOW_VERIFIED_AN_HALF = True          # N = AN/2 only for the eQTL Catalogue schema ('an' = total allele number)

# ---- GWAS side --------------------------------------------------------------------------------------
ALLOW_UNREVIEWED_GWAS_TYPE = False     # Unknown / NEEDS_REVIEW studies stay NOT_TESTED

# ---- overlap with Step09 / Step10 ---------------------------------------------------------------------
ALLOW_DUPLICATES_OF_STEP09_10 = False  # same dataset id already in a Step09 / Step10 snapshot -> skipped

# ---- discovery and harmonisation ----------------------------------------------------------------------
CANDIDATE_P = 1e-5           # molecular traits with P <= this in the locus window are tested
MAX_TRAITS = 100             # per dataset and locus (smallest P first)
MIN_SHARED = 50              # allele-matched shared variants needed for coloc.abf
MIN_OVERLAP = 0.5            # shared variants / variants of the smaller side
INCLUDE_INDELS = True        # indels only match with identical alleles (no left-normalisation)
PALINDROMIC_POLICY = "drop_high_maf"   # keep | drop_high_maf (QTL MAF > PALINDROMIC_MAX_MAF or unknown) | drop
PALINDROMIC_MAX_MAF = 0.42
MAX_REGION_ROWS = 3000000    # safety limit on candidate-trait rows read per dataset and locus
VALIDATION_ROWS = 2000       # rows checked per dataset by --scan
REF_CHECK_MIN_MATCH = 0.95   # REF = GRCh38 base in at least this fraction of SNPs
STATS_MAX_INCONSISTENT = 0.10   # P not reproduced by beta / SE in more rows -> wrong columns
COMPLEX_LD_REGIONS = {       # GRCh38; loci overlapping these get ABF_ASSUMPTION_WARNING
    "MHC": (6, 28510120, 33480577),
    "17q21.31_MAPT_inversion": (17, 45300000, 46700000),
    "8p23.1_inversion": (8, 8100000, 12000000),
}

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

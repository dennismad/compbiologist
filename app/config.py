from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_DATA_PATH = DATA_DIR / "raw" / "uniprot_kinase.tsv"
PROCESSED_DATA_PATH = DATA_DIR / "processed" / "processed_proteins.csv"
SUMMARY_PATH = DATA_DIR / "processed" / "summary.json"
SAMPLE_DATA_PATH = BASE_DIR / "app" / "sample_data" / "uniprot_kinase_sample.tsv"

GEO_RAW_JSON_PATH = DATA_DIR / "raw" / "geo_last_search.json"
GEO_PROCESSED_PATH = DATA_DIR / "processed" / "geo_datasets.csv"
GEO_SUMMARY_PATH = DATA_DIR / "processed" / "geo_summary.json"
GEO_LOADED_JSON_PATH = DATA_DIR / "processed" / "geo_loaded.json"
GEO_LOADED_CSV_PATH = DATA_DIR / "processed" / "geo_loaded.csv"
GEO_MATRIX_CACHE_DIR = DATA_DIR / "raw" / "geo_matrix"
GEO_CACHE_DB_PATH = DATA_DIR / "processed" / "geo_cache.json"
GEO_SAMPLE_PATH = BASE_DIR / "app" / "sample_data" / "geo_sample.json"

ANALYSIS_RESULT_PATH = DATA_DIR / "processed" / "analysis_result.json"
ANALYSIS_DGE_PATH = DATA_DIR / "processed" / "analysis_dge.csv"

UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_QUERY = "reviewed:true AND organism_id:9606 AND keyword:KW-0418"
DEFAULT_FETCH_SIZE = 250

GEO_ESEARCH_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
GEO_ESUMMARY_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GEO_DEFAULT_QUERY = "single cell breast cancer"
GEO_DEFAULT_FETCH_SIZE = 15

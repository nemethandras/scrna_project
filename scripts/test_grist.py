from dotenv import load_dotenv
load_dotenv()   # reads .env from the project root automatically

# scripts/test_grist.py
import requests, os
from grist_api import GristDocAPI

# Assign to variables first
server  = os.environ["GRIST_SERVER"]
doc_id  = os.environ["GRIST_DOC_ID"]
api_key = os.environ["GRIST_API_KEY"]

api = GristDocAPI(doc_id, server=server, api_key=api_key)

resp = requests.get(
    f"{server}/api/docs/{doc_id}/tables/Runs/columns",
    headers={"Authorization": f"Bearer {api_key}"}
)

# Print just the id and label for each column — easy to read
for col in resp.json()["columns"]:
    print(f"  id: {col['id']:<30} label: {col['fields']['label']}")

for table in ["Samples", "QC_Summary"]:
    resp = requests.get(
        f"{server}/api/docs/{doc_id}/tables/{table}/columns",
        headers={"Authorization": f"Bearer {api_key}"}
    )
    print(f"\n{table}:")
    for col in resp.json()["columns"]:
        print(f"  id: {col['id']:<35} label: {col['fields']['label']}")

""" # Try pushing one test row to Runs
test_row = api.add_records("Runs", [{
    "sample_name":             "TEST_SAMPLE",
    "run_id":                  "test_run_001",
    "pipeline_version":        "0.1.0",
    "reference_genome":        "GRCh38",
    "sequencing_mode":         "single",
    "run_date":                "2026-03-25",
    "mapping_rate":            5.48,
    "total_raw_variants":      2855,
    "total_filtered_variants": 2124,
    "vcf_file_path":           "results/SRR1258218_chr22/vcf/SRR1258218.filtered.vcf",
}])

print(f"Row added with ID: {test_row}")

# Fetch it back to confirm
rows = api.fetch_table("Runs")
print(f"Rows in Runs table: {len(rows)}")
print(rows[-1]) """
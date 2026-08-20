import pandas as pd
from bs4 import BeautifulSoup
import os
import sys

# Usage: python extract_amplicon_tables.py report.html

# Parse arguments
if len(sys.argv) < 2:
    print("Usage: python extract_amplicon_tables.py <report.html> [output_folder]")
    sys.exit(1)

html_path = sys.argv[1]
outdir = sys.argv[2] if len(sys.argv) > 2 else "extracted_tables"
os.makedirs(outdir, exist_ok=True)

# Load HTML
with open(html_path, "r", encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

# Find all tables
tables = soup.find_all("table")
print(f"Found {len(tables)} tables.")

# Extract tables
for i, table in enumerate(tables, 1):
    df = pd.read_html(str(table))[0]
    out_file = os.path.join(outdir, f"table_{i}.csv")
    df.to_csv(out_file, index=False)
    print(f"Saved: {out_file}")

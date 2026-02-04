#!/usr/bin/env bash
set -euo pipefail

mkdir tnbc_manifest_output
cd tnbc_manifest_output

###(1) DOWNLOAD CLINICAL MATRIX FROM XENA###
echo "(1/5) Downloading TCGA-BRCA clinical matrix from UCSC Xena..."

URL=(
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
)

for url in "${URL[@]}"; do
    echo "From URL: $url"

    curl -L -H "User-Agent: Mozilla/5.0 (research-download)" \
         -H "Accept: text/plain" \
         -o BRCA_clinicalMatrix.txt \
         "$url"

    # Check if download was successful
    if [[ -s BRCA_clinicalMatrix.txt ]] && ! grep -q "Access Denied" BRCA_clinicalMatrix.txt; then
        echo "Successfully downloaded from: $url"
        break
    else
        echo "Failed to download from: $url" >&2
        rm -f BRCA_clinicalMatrix.txt
    fi
done

###(2) EXTRACT TRIPLE-NEGATIVE IHC SAMPLES###
echo "(2/5) Extracting triple-negative samples (ER-, PR-, HER2-) from Xena matrix..."

# Find samples with triple-negative IHC results
awk -F'\t' '
  BEGIN {OFS="\t"}
  NR==1 {
    # From columns
    for (i=1; i<=NF; i++) {
      if ($i=="ER_Status_nature2012") er=i
      if ($i=="PR_Status_nature2012") pr=i
      if ($i=="HER2_Final_Status_nature2012") her2=i
    }
  }
  NR>1 {
    if ($er=="Negative" && $pr=="Negative" && $her2=="Negative") {
      print $1
    }
  }
' BRCA_clinicalMatrix.txt | sort -u > tnbc_xena_barcodes.txt

echo "TNBC sample count from Xena: $(wc -l < tnbc_xena_barcodes.txt)"

if [[ ! -s tnbc_xena_barcodes.txt ]]; then
  echo "ERROR: No TNBC samples detected. Check column names in Xena file." >&2
  exit 1
fi

###(3) CONVERT XENA BARCODES TO GDC UUIDS###
echo "(3/5) Converting Xena TNBC barcodes to GDC case UUIDs..."

# Convert Xena barcodes to case-level TCGA barcodes
cut -d- -f1-3 "tnbc_xena_barcodes.txt" \
  | sed 's/\r$//' \
  | sort -u > tnbc_case_barcodes.txt

# Check for case-level barcodes
if [[ ! -s tnbc_case_barcodes.txt ]]; then
  echo "ERROR: No case-level barcodes extracted" >&2
  exit 1
fi

# Build JSON array for GDC API
jq_query=$(jq -R -s -c 'split("\n")[:-1]' tnbc_case_barcodes.txt)

# Query GDC /cases endpoint
cases_json=$(curl -s -X POST "https://api.gdc.cancer.gov/cases" \
  -H "Content-Type: application/json" \
  -d "{
        \"filters\": {
          \"op\": \"in\",
          \"content\": {
            \"field\": \"submitter_id\",
            \"value\": $jq_query
          }
        },
        \"fields\": \"case_id\",
        \"size\": 5000
      }")

# Extract case UUIDs
output_file="tnbc_case_uuids.txt"
jq -r '.data.hits[].case_id' <<< "$cases_json" | sort -u > "$output_file"

if [[ ! -s "$output_file" ]]; then
  echo "ERROR: No case UUIDs found for submitted barcodes." >&2
  exit 1
fi

echo "Mapped TNBC cases from Xena to GDC UUIDs: $(wc -l < "$output_file") UUIDs"

###(4) PULL MAF FILES FOR TNBC UUIDS###
echo "(4/5) Querying MAF files for TNBC cases..."

INPUT="tnbc_case_uuids.txt"

# Clean and batch case UUIDs (25 per batch to avoid GDC limits)
sed 's/\r$//' "$INPUT" | sort -u > tnbc_case_uuids.clean.txt
split -l 25 tnbc_case_uuids.clean.txt tnbc_cases_batch_

# Define function to run query by ancestry
run_by_race () {
  RACE="$1"
  OUT_PREFIX="$2"

# Output intermediate and final files
  PAIRS="${OUT_PREFIX}_case_file_pairs.txt"
  FINAL="${OUT_PREFIX}_one_per_case.txt"

  > "$PAIRS"
  > "$FINAL"

# Query GDC in batches
  for batch in tnbc_cases_batch_*; do
    case_ids_json=$(jq -R -s -c 'split("\n")[:-1]' "$batch")

    curl -s -X POST https://api.gdc.cancer.gov/files \
      -H "Content-Type: application/json" \
      -d "{
            \"filters\": {
              \"op\": \"and\",
              \"content\": [
                {
                  \"op\": \"in\",
                  \"content\": {
                    \"field\": \"cases.case_id\",
                    \"value\": $case_ids_json
                  }
                },
                {
                  \"op\": \"=\",
                  \"content\": {
                    \"field\": \"cases.project.project_id\",
                    \"value\": \"TCGA-BRCA\"
                  }
                },
                {
                  \"op\": \"=\",
                  \"content\": {
                    \"field\": \"analysis.workflow_type\",
                    \"value\": \"Aliquot Ensemble Somatic Variant Merging and Masking\"
                  }
                },
                {
                  \"op\": \"=\",
                  \"content\": {
                    \"field\": \"access\",
                    \"value\": \"open\"
                  }
                },
                {
                  \"op\": \"=\",
                  \"content\": {
                    \"field\": \"cases.demographic.race\",
                    \"value\": \"$RACE\"
                  }
                }
              ]
            },
            \"fields\": \"file_id,cases.case_id\",
            \"size\": 1000
          }" \
    | jq -r '.data.hits[] | "\(.cases[0].case_id)\t\(.file_id)"' \
    >> "$PAIRS"
  done

# Collapse to exactly one MAF file per case
  sort -u -k1,1 "$PAIRS" | cut -f2 > "$FINAL"

# Report results
  echo "TNBC cases queried: $(wc -l < tnbc_case_uuids.clean.txt)"

  echo "$RACE TNBC MAF files (one per case): $(wc -l < "$FINAL")"
}

run_by_race "white" "tnbc_mutations_white"
run_by_race "black or african american" "tnbc_mutations_black"

###(5) GENERATE GDC MANIFEST FILE FOR TNBC SAMPLES###

echo "(5/5) Generating GDC manifest file..."

# Define race-specific manifest function
make_manifest () {
  INPUT_IDS="$1"
  OUT_MANIFEST="$2"

# Convert to JSON array
  jq -R -s -c 'split("\n")[:-1]' "$INPUT_IDS" > file_ids.json

#Query metadata
  curl -s -X POST https://api.gdc.cancer.gov/files \
    -H "Content-Type: application/json" \
    -d "{
          \"filters\": {
            \"op\": \"in\",
            \"content\": {
              \"field\": \"file_id\",
              \"value\": $(cat file_ids.json)
            }
          },
          \"fields\": \"file_id,file_name,md5sum,file_size\",
          \"size\": 500
        }" \
  | jq -r '
    .data.hits[]
    | [.file_id, .file_name, .md5sum, .file_size]
    | @tsv
  ' > manifest_body.tsv

# Add manifest header
  printf "id\tfilename\tmd5\tsize\n" > "$OUT_MANIFEST"
  cat manifest_body.tsv >> "$OUT_MANIFEST"

  if [[ $(wc -l < "$OUT_MANIFEST") -le 1 ]]; then
    echo "ERROR: Unable to generate manifest file: $OUT_MANIFEST" >&2
    exit 1
  fi

  echo "Created GDC manifest file: $OUT_MANIFEST ($(($(wc -l < "$OUT_MANIFEST") - 1)) files)"

  rm -f manifest_body.tsv file_ids.json

}

# Generate manifest files
make_manifest \
  tnbc_mutations_white_one_per_case.txt \
  gdc_manifest_tnbc_white_somatic_mutations.tsv

make_manifest \
  tnbc_mutations_black_one_per_case.txt \
  gdc_manifest_tnbc_black_somatic_mutations.tsv

##############################################################

#Need to download from manifest files in separate directories
#gdc-client download -m gdc_manifest_tnbc...

#Need to gunzip all files in Bash with gunzip:
#find . -type f -name "*.gz" -exec gunzip {} +

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
import os
from datetime import datetime

class BRCAMatrixExtractor:
    def __init__(self):
        self.clinical_data = None
        self.extracted_data = None
        self.missing_placeholder = "~"  # Unused keyboard symbol for missing values
        
    def read_clinical_matrix(self, clinical_file):
        """
        (1) Read TCGA-BRCA clinical matrix file
        """
        print(f"(1)Reading clinical matrix: {clinical_file}")
        
        if not os.path.exists(clinical_file):
            print(f"Clinical matrix file not found: {clinical_file}")
            
            """
            # Download matrix from Xena via Bash, if needed
            URL=(
                "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
            )

            for url in "${URL[@]}"; do
               echo "From URL: $url"

               curl -L -H "User-Agent: Mozilla/5.0 (research-download)" \
                 -H "Accept: text/plain" \
                 -o BRCA_clinicalMatrix.txt \
                 "$url"
             """

            return False
        
        try:
            # Read the clinical matrix (tab-separated)
            self.clinical_data = pd.read_csv(clinical_file, sep='\t', index_col=0, low_memory=False)
            print(f"Loaded clinical data: {self.clinical_data.shape[0]} samples, {self.clinical_data.shape[1]} features")
            
            # Check for required IHC columns
            required_ihc_columns = [
                'ER_Status_nature2012',
                'PR_Status_nature2012', 
                'HER2_Final_Status_nature2012'
            ]
            
            missing_columns = [col for col in required_ihc_columns if col not in self.clinical_data.columns]
            
            if missing_columns:
                print(f"Missing required IHC columns: {missing_columns}")
                return False
            else:
                print(f"All required IHC columns found.")
                
            return True
            
        except Exception as e:
            print(f"Error reading clinical matrix: {e}")
            return False
    
    def validate_user_columns(self, user_columns):
        """
        (2) Validate that user-defined columns exist in the clinical matrix
        """
        print(f"(2) Validating user-defined columns...")
        
        available_columns = self.clinical_data.columns.tolist()
        valid_columns = []
        invalid_columns = []
        
        for col in user_columns:
            if col in available_columns:
                valid_columns.append(col)
                # Count non-empty values
                non_empty_count = self.clinical_data[col].dropna().shape[0]
                non_empty_non_blank = (
                    (self.clinical_data[col].notna()) & 
                    (self.clinical_data[col] != '') & 
                    (self.clinical_data[col] != 'NA')
                ).sum()
                print(f"   {col}: {non_empty_non_blank} non-empty values")
            else:
                invalid_columns.append(col)
                print(f"   {col}: Column not found")
        
        if invalid_columns:
            print(f"\nAvailable columns containing similar terms:")
            for invalid_col in invalid_columns:
                # Try to find similar column names
                similar_cols = [col for col in available_columns 
                               if any(term.lower() in col.lower() 
                                     for term in invalid_col.replace('_', ' ').split())]
                if similar_cols:
                    print(f"   For '{invalid_col}', similar columns found:")
                    for sim_col in similar_cols[:5]:  # Show first 5
                        print(f"     - {sim_col}")
        
        return valid_columns, invalid_columns
    
    def filter_complete_ihc_samples(self):
        """
        (3) Filter samples with non-empty values in all required IHC columns
        """
        print(f"(3) Filtering samples with complete IHC data...")
        
        # Define the required columns
        ihc_columns = [
            'ER_Status_nature2012',
            'PR_Status_nature2012', 
            'HER2_Final_Status_nature2012'
        ]
        
        # Start with all samples
        total_samples = len(self.clinical_data)
        print(f"Total samples in clinical matrix: {total_samples}")
        
        # Filter for non-empty values in IHC columns
        mask = True
        for col in ihc_columns:
            col_mask = (
                self.clinical_data[col].notna() & 
                (self.clinical_data[col] != '') & 
                (self.clinical_data[col] != 'NA')
            )
            mask = mask & col_mask
            
            non_empty_count = col_mask.sum()
            print(f"   {col}: {non_empty_count} non-empty values")
        
        # Apply the filter
        self.filtered_clinical = self.clinical_data[mask].copy()
        print(f"Samples with complete IHC data: {len(self.filtered_clinical)}")
        
        return True
    
    def standardize_receptor_status(self, value):
        """
        (3.5) Standardize receptor status values to +/- format
        """
        if pd.isna(value) or value == '' or value == 'NA':
            return 'Unknown'
        
        value = str(value).upper().strip()
        
        # Common positive indicators
        if any(pos in value for pos in ['POSITIVE', 'POS', '+', '1', 'YES', 'TRUE']):
            return '+'
        # Common negative indicators  
        elif any(neg in value for neg in ['NEGATIVE', 'NEG', '-', '0', 'NO', 'FALSE']):
            return '-'
        else:
            return value  # Keep original if can't standardize
    
    def create_receptor_pattern(self, er, pr, her2):
        """
        (3.6) Create standardized receptor status pattern (e.g., +/+/-)
        """
        er_std = self.standardize_receptor_status(er)
        pr_std = self.standardize_receptor_status(pr)
        her2_std = self.standardize_receptor_status(her2)
        
        return f"{er_std}/{pr_std}/{her2_std}"
    
    def extract_data(self, user_columns):
        """
        (4) Extract the required data into the output format
        """
        print(f"(4) Extracting data for output...")
        
        results = []
        
        for sample_id in self.filtered_clinical.index:
            # Get IHC data
            er_status = self.filtered_clinical.loc[sample_id, 'ER_Status_nature2012']
            pr_status = self.filtered_clinical.loc[sample_id, 'PR_Status_nature2012']
            her2_status = self.filtered_clinical.loc[sample_id, 'HER2_Final_Status_nature2012']
            
            # Standardize receptor status
            er_std = self.standardize_receptor_status(er_status)
            pr_std = self.standardize_receptor_status(pr_status)
            her2_std = self.standardize_receptor_status(her2_status)
            
            # Create receptor pattern
            receptor_pattern = self.create_receptor_pattern(er_status, pr_status, her2_status)
            
            # Start building the row
            row_data = {
                'SampleID': sample_id,
                'ER_Status': er_std,
                'PR_Status': pr_std,
                'HER2_Status': her2_std,
                'Receptor_Pattern': receptor_pattern
            }
            
            # Add user-defined columns
            for col in user_columns:
                if col in self.filtered_clinical.columns:
                    value = self.filtered_clinical.loc[sample_id, col]
                    if pd.isna(value) or value == '' or value == 'NA':
                        row_data[col] = self.missing_placeholder
                    else:
                        row_data[col] = str(value)
                else:
                    row_data[col] = self.missing_placeholder
            
            results.append(row_data)
        
        # Convert to DataFrame
        self.extracted_data = pd.DataFrame(results)
        print(f"Extracted data for {len(self.extracted_data)} samples")
        
        return True
    
    def save_output(self, output_file):
        """
        (5) Save the extracted data to TSV file
        """
        print(f"(5) Saving output to: {output_file}")
        
        try:
            self.extracted_data.to_csv(output_file, sep='\t', index=False)
            print(f"Output saved successfully")
            
            # Print summary statistics
            print(f"\nOutput Summary:")
            print(f"   Total samples: {len(self.extracted_data)}")
            print(f"   Total columns: {len(self.extracted_data.columns)}")
            
            # Show receptor pattern distribution
            pattern_counts = self.extracted_data['Receptor_Pattern'].value_counts()
            print(f"\nReceptor Pattern Distribution:")
            for pattern, count in pattern_counts.items():
                print(f"   {pattern}: {count} samples")
            
            # Show first few rows
            print(f"\nSample of output (first 3 rows):")
            print(self.extracted_data.head(3).to_string(index=False))
            
            return True
            
        except Exception as e:
            print(f"Error saving output: {e}")
            return False
    
    def list_available_columns(self, search_term=None):
        """
        (5.5) List available columns in the clinical matrix, optionally filtered by search term
        """
        columns = self.clinical_data.columns.tolist()
        
        if search_term:
            filtered_cols = [col for col in columns if search_term.lower() in col.lower()]
            print(f"\n📋 Available columns containing '{search_term}' ({len(filtered_cols)}):")
            for col in sorted(filtered_cols):
                non_empty = (
                    (self.clinical_data[col].notna()) & 
                    (self.clinical_data[col] != '') & 
                    (self.clinical_data[col] != 'NA')
                ).sum()
                print(f"   - {col} ({non_empty} non-empty)")
        else:
            print(f"\n📋 All available columns ({len(columns)}):")
            for col in sorted(columns):
                print(f"   - {col}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract BRCA clinical matrix data with flexible column selection'
    )
    parser.add_argument(
        '-i', '--input', 
        default='BRCA_clinicalMatrix.txt',
        help='Input clinical matrix file (default: BRCA_clinicalMatrix.txt)'
    )
    parser.add_argument(
        '-o', '--output',
        default='brca_extracted_data.tsv',
        help='Output TSV file (default: brca_extracted_data.tsv)'
    )
    parser.add_argument(
        '-c', '--columns',
        nargs='*',
        default=[],
        help='User-defined columns to extract (space-separated, max 10)'
    )
    parser.add_argument(
        '--list-columns',
        action='store_true',
        help='List all available columns and exit'
    )
    parser.add_argument(
        '--search-columns',
        type=str,
        help='Search for columns containing specified term'
    )
    parser.add_argument(
        '--placeholder',
        default='~',
        help='Placeholder for missing values (default: ~)'
    )
    
    args = parser.parse_args()
    
    # Initialize extractor
    extractor = BRCAMatrixExtractor()
    extractor.missing_placeholder = args.placeholder
    
    # Read clinical matrix
    if not extractor.read_clinical_matrix(args.input):
        sys.exit(1)
    
    # Handle column listing/searching
    if args.list_columns:
        extractor.list_available_columns()
        sys.exit(0)
    
    if args.search_columns:
        extractor.list_available_columns(args.search_columns)
        sys.exit(0)
    
    # Validate user-defined columns
    if len(args.columns) > 10:
        print(f"Too many user-defined columns ({len(args.columns)}). Maximum allowed: 10")
        sys.exit(1)
    
    if args.columns:
        valid_columns, invalid_columns = extractor.validate_user_columns(args.columns)
        if invalid_columns:
            print(f"Invalid columns found: {invalid_columns}")
            print(f"Use --search-columns <term> to find similar column names")
            sys.exit(1)
        user_columns = valid_columns
    else:
        print("No user-defined columns specified. Only IHC data will be extracted.")
        user_columns = []
    
    # Filter samples with complete IHC data
    if not extractor.filter_complete_ihc_samples():
        sys.exit(1)
    
    # Extract data
    if not extractor.extract_data(user_columns):
        sys.exit(1)
    
    # Save output
    if not extractor.save_output(args.output):
        sys.exit(1)
    
    print(f"\nAnalysis completed successfully!")
    print(f"Output file: {args.output}")

if __name__ == "__main__":
    main()

"""
Columns recommended for survival analysis: pathologic_stage, Days_to_date_of_Death_nature2012, OS_Time_nature2012, OS_event_nature2012, Vital_Status_nature2012
"""

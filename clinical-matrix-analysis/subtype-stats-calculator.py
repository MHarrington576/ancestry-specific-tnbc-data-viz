#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
from collections import defaultdict
from datetime import datetime

class TSVReceptorAnalyzer:
    
    def __init__(self):
        self.data = None
        self.receptor_patterns = ['+/+/-', '+/-/-', '-/-/+', '-/-/-']
        self.analysis_results = {}
        
    def read_tsv_file(self, input_file):
        """
       (1)  Read TSV file and validate required columns
        """
        print(f"(1) Reading TSV file: {input_file}")
        
        try:
            self.data = pd.read_csv(input_file, sep='\t')
            print(f"Loaded data: {len(self.data)} samples, {len(self.data.columns)} columns")
            
            # Check for required columns
            required_cols = ['SampleID', 'ER_Status', 'PR_Status', 'HER2_Status', 'Receptor_Pattern']
            missing_cols = [col for col in required_cols if col not in self.data.columns]
            
            if missing_cols:
                print(f"Missing required columns: {missing_cols}")
                return False
                
            # Show available columns for analysis
            available_cols = [col for col in self.data.columns if col not in required_cols]
            print(f"Available columns for analysis: {len(available_cols)}")
            for i, col in enumerate(available_cols, 1):
                non_empty = self.data[col].notna().sum()
                print(f"   {i:2d}. {col} ({non_empty} non-empty values)")
                
            return True
            
        except Exception as e:
            print(f"Error reading file: {e}")
            return False
    
    def validate_analysis_columns(self, analysis_columns):
        """
        (2) Validate that analysis columns exist in the data
        """
        print(f"(2) Validating analysis columns...")
        
        available_cols = self.data.columns.tolist()
        valid_columns = []
        invalid_columns = []
        
        for col in analysis_columns:
            if col in available_cols:
                non_empty = self.data[col].notna().sum()
                print(f"   {col}: {non_empty} non-empty values")
                valid_columns.append(col)
            else:
                print(f"   {col}: Column not found")
                invalid_columns.append(col)
        
        if invalid_columns:
            print(f"\n📋 Available columns:")
            for col in available_cols[5:]:  # Skip the first 5 required columns
                print(f"   - {col}")
            
        return valid_columns
    
    def analyze_receptor_patterns(self, analysis_columns):
        """
        (3) Perform analysis of receptor patterns vs specified columns
        """
        print(f"(3) Analyzing receptor patterns vs specified columns...")
        
        # Get receptor pattern distribution
        pattern_counts = self.data['Receptor_Pattern'].value_counts()
        print(f"\nOverall IHC Receptor Pattern Distribution:")
        
        for pattern in self.receptor_patterns:
            count = pattern_counts.get(pattern, 0)
            percentage = (count / len(self.data)) * 100 if len(self.data) > 0 else 0
            print(f"   {pattern:8s}: {count:4d} samples ({percentage:5.1f}%)")
        
        total_samples = pattern_counts.sum()
        print(f"   {'Total':8s}: {total_samples:4d} samples")
        
        # Analyze each specified column
        self.analysis_results = {}
        
        for col in analysis_columns:
            print(f"\nAnalyzing column: {col}")
            
            self.analysis_results[col] = {}
            
            # Get unique values in this column (excluding NaN)
            unique_values = self.data[col].dropna().unique()
            print(f"Unique values in {col}: {len(unique_values)}")
            
            if len(unique_values) > 10:
                print(f"   (Showing first 10): {list(unique_values[:10])}")
            else:
                print(f"   Values: {list(unique_values)}")
            
            # For each receptor pattern, analyze the distribution
            for pattern in self.receptor_patterns:
                pattern_data = self.data[self.data['Receptor_Pattern'] == pattern]
                pattern_total = len(pattern_data)
                
                if pattern_total == 0:
                    continue
                    
                print(f"\nIHC Pattern: {pattern} ({pattern_total} total samples)")
                
                self.analysis_results[col][pattern] = {
                    'total_samples': pattern_total,
                    'value_distribution': {},
                    'non_empty_total': 0
                }
                
                # Count non-empty values for this pattern
                pattern_col_data = pattern_data[col].dropna()
                non_empty_count = len(pattern_col_data)
                self.analysis_results[col][pattern]['non_empty_total'] = non_empty_count
                
                if non_empty_count == 0:
                    print(f"   No data available for {col}")
                    continue
                
                # Count each unique value
                value_counts = pattern_col_data.value_counts()
                
                for value, count in value_counts.items():
                    # Calculate percentage based on total samples in this pattern
                    percentage_of_pattern = (count / pattern_total) * 100
                    # Calculate percentage based on non-empty samples
                    percentage_of_non_empty = (count / non_empty_count) * 100
                    
                    print(f"   {str(value):20s}: {count:3d} samples "
                          f"({percentage_of_pattern:5.1f}% of {pattern} total, "
                          f"{percentage_of_non_empty:5.1f}% of non-empty)")
                    
                    self.analysis_results[col][pattern]['value_distribution'][str(value)] = {
                        'count': count,
                        'percentage_of_pattern_total': percentage_of_pattern,
                        'percentage_of_non_empty': percentage_of_non_empty
                    }
                
                # Show missing data info
                missing_count = pattern_total - non_empty_count
                if missing_count > 0:
                    missing_percentage = (missing_count / pattern_total) * 100
                    print(f"   {'Missing/Empty':20s}: {missing_count:3d} samples "
                          f"({missing_percentage:5.1f}% of {pattern} total)")
    
    def generate_summary_report(self, output_file=None):
        """
        (4) Generate a comprehensive summary report
        """
        print(f"\n(4) ANALYSIS SUMMARY")
        print("=" * 15)
        
        summary_lines = []
        summary_lines.append(f"# TSV Receptor Pattern Analysis Summary")
        summary_lines.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        summary_lines.append(f"# Total samples: {len(self.data)}")
        summary_lines.append(f"# Columns analyzed: {len(self.analysis_results)}")
        summary_lines.append("")
        
        # Overall receptor pattern distribution
        summary_lines.append("## Overall IHC Receptor Pattern Distribution")
        pattern_counts = self.data['Receptor_Pattern'].value_counts()
        
        for pattern in self.receptor_patterns:
            count = pattern_counts.get(pattern, 0)
            percentage = (count / len(self.data)) * 100 if len(self.data) > 0 else 0
            summary_lines.append(f"{pattern}: {count} samples ({percentage:.1f}%)")
        
        summary_lines.append("")
        
        # Detailed analysis for each column
        for col, col_results in self.analysis_results.items():
            summary_lines.append(f"## Analysis: {col}")
            summary_lines.append("")
            
            # Create summary table
            summary_lines.append("| IHC_Pattern | Total_Samples | Non_Empty | Value | Count | % of Pattern | % of Non_Empty |")
            summary_lines.append("|-------------|---------------|-----------|-------|-------|--------------|----------------|")
            
            for pattern in self.receptor_patterns:
                if pattern not in col_results:
                    continue
                    
                pattern_data = col_results[pattern]
                total = pattern_data['total_samples']
                non_empty = pattern_data['non_empty_total']
                
                if not pattern_data['value_distribution']:
                    summary_lines.append(f"| {pattern} | {total} | {non_empty} | No data | - | - | - |")
                    continue
                
                first_row = True
                for value, value_data in pattern_data['value_distribution'].items():
                    count = value_data['count']
                    perc_pattern = value_data['percentage_of_pattern_total']
                    perc_non_empty = value_data['percentage_of_non_empty']
                    
                    if first_row:
                        summary_lines.append(f"| {pattern} | {total} | {non_empty} | {value} | {count} | {perc_pattern:.1f}% | {perc_non_empty:.1f}% |")
                        first_row = False
                    else:
                        summary_lines.append(f"| | | | {value} | {count} | {perc_pattern:.1f}% | {perc_non_empty:.1f}% |")
            
            summary_lines.append("")
        
        # Print to console
        for line in summary_lines:
            print(line)
        
        # Save to file if requested
        if output_file:
            try:
                with open(output_file, 'w') as f:
                    f.write('\n'.join(summary_lines))
                print(f"\nSummary report saved to: {output_file}")
                return True
            except Exception as e:
                print(f"Error saving summary report: {e}")
                return False
        
        return True
    
    def generate_detailed_statistics_file(self, output_file):
        """
        (4.5) Generate detailed statistics as TSV
        """
        print(f"\nGenerating detailed statistics file: {output_file}")
        
        try:
            with open(output_file, 'w') as f:
                f.write("# TCGA-BRCA TSV Receptor Pattern Detailed Statistics\n")
                f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"# Total samples: {len(self.data)}\n")
                f.write(f"# Columns analyzed: {len(self.analysis_results)}\n")
                f.write("#\n")
                
                # Write detailed statistics
                f.write("Column\tIHC_Pattern\tTotal_Samples_in_Pattern\tColumn_Value\tCount\tPercentage_of_Pattern_Total\tPercentage_of_Non_Empty\tNon_Empty_Total\n")
                
                for col, col_results in self.analysis_results.items():
                    for pattern in self.receptor_patterns:
                        if pattern not in col_results:
                            continue
                            
                        pattern_data = col_results[pattern]
                        total = pattern_data['total_samples']
                        non_empty = pattern_data['non_empty_total']
                        
                        if not pattern_data['value_distribution']:
                            f.write(f"{col}\t{pattern}\t{total}\tNo_Data\t0\t0.0\t0.0\t{non_empty}\n")
                            continue
                        
                        for value, value_data in pattern_data['value_distribution'].items():
                            count = value_data['count']
                            perc_pattern = value_data['percentage_of_pattern_total']
                            perc_non_empty = value_data['percentage_of_non_empty']
                            
                            # Clean value for TSV (replace tabs and newlines)
                            clean_value = str(value).replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
                            
                            f.write(f"{col}\t{pattern}\t{total}\t{clean_value}\t{count}\t{perc_pattern:.2f}\t{perc_non_empty:.2f}\t{non_empty}\n")
                
            print("Detailed statistics file generated successfully.")
            return True
            
        except Exception as e:
            print(f"Error generating detailed statistics file: {e}")
            return False

def main():
    parser = argparse.ArgumentParser(
        description='Analyze TSV files for IHC receptor pattern statistics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single column
  python3 %(prog)s -i extracted_data.tsv -c OS_Status
  
  # Analyze multiple columns  
  python3 %(prog)s -i extracted_data.tsv -c OS_Status Vital_Status AJCC_Stage_nature2012
  
  # Generate detailed output files
  python3 %(prog)s -i extracted_data.tsv -c OS_Status -o analysis_summary.md --detailed stats_detailed.tsv
  
  # List available columns first
  python3 %(prog)s -i extracted_data.tsv --list-columns
        """
    )
    
    parser.add_argument('-i', '--input', default='extracted_clinical_data.tsv',
                       help='Input TSV file (default: extracted_clinical_data.tsv)')
    parser.add_argument('-c', '--columns', nargs='+', 
                       help='Column names to analyze (space-separated)')
    parser.add_argument('-o', '--output', 
                       help='Output summary report file (markdown format)')
    parser.add_argument('--detailed',
                       help='Generate detailed statistics TSV file')
    parser.add_argument('--list-columns', action='store_true',
                       help='List available columns and exit')
    
    args = parser.parse_args()
    
    print("IHC Pattern Statistics Calculator")
    print("=" * 15)
    
    # Initialize analyzer
    analyzer = TSVReceptorAnalyzer()
    
    # Read input file
    if not analyzer.read_tsv_file(args.input):
        sys.exit(1)
    
    # List columns if requested
    if args.list_columns:
        print(f"\nAll available columns in {args.input}:")
        required_cols = ['SampleID', 'ER_Status', 'PR_Status', 'HER2_Status', 'Receptor_Pattern']
        available_cols = [col for col in analyzer.data.columns if col not in required_cols]
        
        for i, col in enumerate(available_cols, 1):
            non_empty = analyzer.data[col].notna().sum()
            unique_vals = analyzer.data[col].nunique()
            print(f"   {i:2d}. {col}")
            print(f"       Non-empty values: {non_empty}")
            print(f"       Unique values: {unique_vals}")
            if unique_vals <= 5:
                unique_list = analyzer.data[col].dropna().unique()
                print(f"       Values: {list(unique_list)}")
        sys.exit(0)
    
    # Check if columns specified
    if not args.columns:
        print("No analysis columns specified. Use -c to specify columns or --list-columns to see available options.")
        sys.exit(1)
    
    # Validate analysis columns
    valid_columns = analyzer.validate_analysis_columns(args.columns)
    if not valid_columns:
        print("No valid columns to analyze.")
        sys.exit(1)
    
    # Perform analysis
    analyzer.analyze_receptor_patterns(valid_columns)
    
    # Generate summary report
    analyzer.generate_summary_report(args.output)
    
    # Generate detailed statistics if requested
    if args.detailed:
        analyzer.generate_detailed_statistics_file(args.detailed)
    
    print(f"\nAnalysis complete!")
    print(f"Analyzed {len(analyzer.data)} samples across {len(valid_columns)} columns")
    
    if args.output:
        print(f"Summary report: {args.output}")
    if args.detailed:
        print(f"Detailed statistics: {args.detailed}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: treat_metadata_ENA.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'treat_metadata_ENA.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, ENA_FASTQ_URLS_COLUMNS, TEMPLATE_FINAL_COLUMNS, 
                    VALID_FASTQ_TYPES, VALID_TREATMENTS, OMD_CTK_Exception, 
                    program_header, check_existence_directory_parameter, 
                    check_fastq_PAIRED_patterns, print_list_n_byline, 
                    check_headers, treat_headers_check, rich_text_colored,
                    check_values, treat_values_check, check_na_in_pandas_dataframe,
                    check_duplicates_in_fastq_names, treat_check_fastq_name_type,
                    check_treatment_for_samples, check_rename_samples, check_merge_samples,
                    treat_output_directory_parameter_outfiles, show_advise_legend)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
import pandas as pd

#Program Constants
DEFAULT_ENA_NO_WARNING_COLUMNS = ['sample_accession', 'secondary_sample_accession', 'experiment_accession', 
                                  'run_accession', 'library_name', 'experiment_title', 'experiment_alias', 
                                  'run_alias', 'sample_alias', 'broker_name','sample_title', 'first_public',
                                  'last_updated', 'ENA-FIRST-PUBLIC', 'ENA-LAST-UPDATE', 'first_created']

ENA_COLUMNS_IGNORED = ['nominal_length', 'read_count', 'base_count',  'fastq_bytes', 'fastq_md5', 'fastq_ftp', 
                       'fastq_aspera', 'fastq_galaxy', 'submitted_bytes', 'submitted_md5', 'submitted_ftp', 
                       'submitted_aspera', 'submitted_galaxy', 'submitted_format', 'sra_bytes',
                       'sra_md5', 'sra_ftp', 'sra_aspera', 'sra_galaxy', 'cram_index_ftp',
                       'cram_index_aspera', 'cram_index_galaxy', 'nominal_sdev', 'Read depth']

DEFAULT_TREATED_METADATA_COL_NAMES = ['final_files_sample_name', 'original_files_sample_names', 'treatment_sample_name', 'treatment_fastq_type']

DEFAULT_WARNING_DF_COL_NAMES = ['final_files_sample_name', 'metadata_column_name', 'warning']

WARNING_MESSAGE = 'More than one value detected after metadata combination and this metadata column was not indicated as "no warning" with the "--extra_no_warning_columns" parameter!'

#Program Global Variables
treated_metadata_df = pd.DataFrame()
warnings_df = pd.DataFrame()

#Program Functions
def check_treatment_fastqs_in_metadata(treatment_df, metadata_df, column, color_treatment):
    """
    This function checks if all the provided fastq file names in the Treatment Template
    are found only one time in the Metadata. Differentiates between zero and more than one match.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    column : str
        The ENA metadata column from which the urls were obtained.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If there are Fastqs names in the treatment dataframe that do not have 
        a unique entry in the metadata raises an exception.

    Returns
    -------
    None.

    """
    #Init results lists
    fastqs_warnings_more_than_one_match = []
    fastqs_warnings_zero_matches = []
    
    #Iter and process fastq files info
    for i,r in treatment_df.iterrows():
        temp_fastq_name = r['fastq_file_name']
        temp_fastq_df = metadata_df[metadata_df[column].str.contains(temp_fastq_name)]
        #Treat info depending on temp_fastq_df len()/ Treat warnings/ Else pass
        if len(temp_fastq_df) == 1:
            pass
        elif len(temp_fastq_df) > 1:
            fastqs_warnings_more_than_one_match.append(temp_fastq_name)
        else:
            fastqs_warnings_zero_matches.append(temp_fastq_name)
    
    #If there are warnings raise exception
    if len(fastqs_warnings_more_than_one_match) > 0 or len(fastqs_warnings_zero_matches) > 0:
        if len(fastqs_warnings_more_than_one_match) > 0:
            print(rich_text_colored('\nWarning! The following Fastq files from the Treatment Template have more than one match in the Metadata Table!', 'program_warning', color_treatment))
            print_list_n_byline(fastqs_warnings_more_than_one_match, 5)
        elif len(fastqs_warnings_zero_matches) > 0:
            print(rich_text_colored('\nWarning! The following Fastq files from the Treatment Template do not have a match in the Metadata Table!', 'program_warning', color_treatment))
            print_list_n_byline(fastqs_warnings_zero_matches, 5)
        raise OMD_CTK_Exception('Error! Some of the Fastqs in the Treatment Template do not have a unique match with the Metadata Table (zero or more than one match)!')        


def unique_original_sample_names(sample_treatment_df, fastq_pattern, sep, n_sep):
    """
    This function generates a list with the unique original sample names from
    the treatment template based on "fastq_file_name" column.

    Parameters
    ----------
    sample_treatment_df : pandas dataframe
        The provided filted sample treatment dataframe.
    fastq_pattern : str
        Extension pattern to recognize a Fastq file.
    sep : str
        Sample Name separator.
    n_sep : int
        Sample Name separator appearance.

    Returns
    -------
    unique_original_names : list
        Unique original sample names.

    """
    #Process fastq_file_name (sample_name + rest + fastq_pattern)
    ##Remove fastq_pattern
    original_fastq_names_without_extension = sample_treatment_df['fastq_file_name'].str.replace(fastq_pattern, '', regex = False)
    
    ##Process fastq names without extension to keep original sample_name
    unique_original_names = list(set(original_fastq_names_without_extension.str.split(sep).apply(lambda x:x[:n_sep]).str.join(sep)))
    unique_original_names.sort()
    
    return unique_original_names


def combine_metadata_rows(final_files_sample_name, filtered_metadata_df, interest_metadata_column, no_warnings_metadata_columns):
    """
    This function combines the metadata (unique values, and treat NAs to '')
    for the provided metadata column. It also treats the associated warnings.

    Parameters
    ----------
    final_files_sample_name : str
        The provided final_files_sample_name.
    filtered_metadata_df : pandas dataframe
        Filtered metadata table of interest.
    interest_metadata_column : str
        Name of the metadata column of interest to combine.
    no_warnings_metadata_columns : list
        List of metadata column names that is normal/expected to have
        multiple values. No warning metadata columns.

    Returns
    -------
    result : str
        Returns the combined metadata. If multiple values are found these will
        be separated by semicolons (;).

    """
    
    #Set global variables
    global warnings_df
    
    #1. Get values of column as list
    metadata_col_values = list(filtered_metadata_df[interest_metadata_column])
    
    #2. Treat NAs and convert to str
    metadata_col_values_treated = []
    for i in metadata_col_values:
        if pd.isna(i):
            metadata_col_values_treated.append('')
        else:
            metadata_col_values_treated.append(str(i))
    
    #3. Get list with unique values 
    unique_values_list = list(set(metadata_col_values_treated))
    unique_values_list.sort()
    
    #4. Check warning for column 
    ##If there is more than one value and col_name not in no_warnings_metadata_columns -> write to warnings_df
    if (len(unique_values_list) > 1) and (interest_metadata_column not in no_warnings_metadata_columns):
        #Create new line in warnings_df 
        ##Get final_files_sample_name / DEFAULT_WARNING_DF_COL_NAMES[0]
        ##Get metadata_column_name / DEFAULT_WARNING_DF_COL_NAMES[1]
        ##Get warning message /  DEFAULT_WARNING_DF_COL_NAMES[2] -> WARNING_MESSAGE
        temp_warning = pd.DataFrame([[final_files_sample_name, interest_metadata_column, WARNING_MESSAGE]], columns = DEFAULT_WARNING_DF_COL_NAMES)
        warnings_df = pd.concat([warnings_df, temp_warning], ignore_index = True)
        
    #5.Convert list to formated str and return 
    result = ';'.join(unique_values_list)
    return result


def generic_mode_metadata(sample_treatment_df, treatment_saple_name, metadata_df, ena_download_column, interest_metadata_columns, no_warnings_metadata_columns, fastq_pattern, sep, n_sep):
    """
    This function creates the necessary combined metadata lines for the provided 
    sample in generic cases (rename_mode and merge_mode). Lines will be saved 
    directy in the treated_metadata_df.

    Parameters
    ----------
    sample_treatment_df : pandas dataframe
        The provided filted sample treatment dataframe.
    treatment_saple_name : str
        Sample Name from the Treatment Template.
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    ena_download_column : str
        Provided ENA Download Column from the Metadata Table.
    interest_metadata_columns : list
        List with the metadata columns of interest to treat.
    no_warnings_metadata_columns : list
        List of metadata column names that is normal/expected to have
        multiple values. No warning metadata columns.
    fastq_pattern : str
        Extension pattern to recognize a Fastq file.
    sep : str
        Sample Name separator.
    n_sep : int
        Sample Name separator appearance.

    Returns
    -------
    None.

    """
    #Set global variables
    global treated_metadata_df
    
    #Get list of unique original sample names for treatment template "fastq_file_name"
    original_sample_names = unique_original_sample_names(sample_treatment_df, fastq_pattern, sep, n_sep)
    
    #0. Set result list
    sample_metadata_line = []
    
    #1. Filter metadata table based in original_sample_names (treatment template) and ena_download_column (metadata)
    ##NOTE: We re-add the sep, to avoid unspecific matches. For example, sample1 could match with sample1, sample10, sample100, etc.
    ##NOTE2: We add the fastq_pattern to try que build the fastq name (in some SE cases, there will not be a separator and this is the way that I came up with to take into account this cases)  
    
    ##Previous steps
    #Get original_sample_names + sep
    original_sample_names_sep = [i + sep for i in original_sample_names]
    #Get original_sample_names + fastq_pattern
    original_sample_names_file = [i + fastq_pattern for i in original_sample_names]
    
    ##Filter metadata
    sample_metadata_df = metadata_df[metadata_df[ena_download_column].str.contains('|'.join(original_sample_names_sep + original_sample_names_file))]
    
    #2. Get DEFAULT_TREATED_METADATA_COL_NAMES
    ##Get final_files_sample_name / DEFAULT_TREATED_METADATA_COL_NAMES[0]
    sample_metadata_line.append(treatment_saple_name)
    ##Get original_files_sample_names / DEFAULT_TREATED_METADATA_COL_NAMES[1]
    original_files_sample_names = ';'.join(original_sample_names)
    sample_metadata_line.append(original_files_sample_names)
    ##Get treatment_sample_name / DEFAULT_TREATED_METADATA_COL_NAMES[2]
    sample_metadata_line.append(treatment_saple_name)
    ##Get unique list for treatment_fastq_type / DEFAULT_TREATED_METADATA_COL_NAMES[2]
    unique_fastq_types = list(set(sample_treatment_df['fastq_type']))
    unique_fastq_types.sort()
    sample_metadata_line.append(';'.join(unique_fastq_types))
    
    #3. Combine each metadata column of interest 
    for col in interest_metadata_columns:
        temp_combined_metadata = combine_metadata_rows(treatment_saple_name, sample_metadata_df, col, no_warnings_metadata_columns)
        sample_metadata_line.append(temp_combined_metadata)
    
    #4. Create new line in treated_metadata_df
    new_metadata_line = pd.DataFrame([sample_metadata_line], columns = DEFAULT_TREATED_METADATA_COL_NAMES + interest_metadata_columns)
    treated_metadata_df = pd.concat([treated_metadata_df, new_metadata_line], ignore_index = True)


def copy_only_mode_metadata(sample_treatment_df, treatment_saple_name, metadata_df, ena_download_column, interest_metadata_columns, no_warnings_metadata_columns, fastq_pattern, sep, n_sep):
    """
    This function creates the necessary combined metadata lines for the provided 
    sample in the case of the copy_only_mode. Lines will be saved directy in the
    treated_metadata_df.

    Parameters
    ----------
    sample_treatment_df : pandas dataframe
        The provided filted sample treatment dataframe.
    treatment_saple_name : str
        Sample Name from the Treatment Template.
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    ena_download_column : str
        Provided ENA Download Column from the Metadata Table.
    interest_metadata_columns : list
        List with the metadata columns of interest to treat.
    no_warnings_metadata_columns : list
        List of metadata column names that is normal/expected to have
        multiple values. No warning metadata columns.
    fastq_pattern : str
        Extension pattern to recognize a Fastq file.
    sep : str
        Sample Name separator.
    n_sep : int
        Sample Name separator appearance.

    Returns
    -------
    None.

    """
    #Set global variables
    global treated_metadata_df
    
    #Get list of unique original sample names for treatment template "fastq_file_name"
    original_sample_names = unique_original_sample_names(sample_treatment_df, fastq_pattern, sep, n_sep)
    
    #Treat diferent original sample names
    for sample in original_sample_names:
        #0. Set result list
        temp_sample_metadata_line = []
        
        #1. Filter metadata table based in original_sample_name (treatment template) and ena_download_column (metadata)
        ##NOTE: We re-add the sep, to avoid unspecific matches. For example, sample1 could match with sample1, sample10, sample100, etc.
        ##NOTE2: We add the fastq_pattern to try que build the fastq name (in some SE cases, there will not be a separator and this is the way that I came up with to take into account this cases)  
        
        ##Previous steps
        #Get sample + sep
        sample_plus_sep = sample + sep 
        #Get sample + fastq_pattern
        sample_plus_extension = sample + fastq_pattern 
        
        ##Filter metadata table
        temp_sample_metadata_df = metadata_df[metadata_df[ena_download_column].str.contains('|'.join([sample_plus_sep, sample_plus_extension]))]
        
        #2. Get DEFAULT_TREATED_METADATA_COL_NAMES
        ##Get final_files_sample_name / DEFAULT_TREATED_METADATA_COL_NAMES[0]
        temp_sample_metadata_line.append(sample)
        ##Get original_files_sample_names / DEFAULT_TREATED_METADATA_COL_NAMES[1]
        temp_sample_metadata_line.append(sample)
        ##Get treatment_sample_name / DEFAULT_TREATED_METADATA_COL_NAMES[2]
        temp_sample_metadata_line.append(treatment_saple_name)
        ##Get unique list for treatment_fastq_type / DEFAULT_TREATED_METADATA_COL_NAMES[2]
        unique_fastq_types = list(set(sample_treatment_df['fastq_type']))
        unique_fastq_types.sort()
        temp_sample_metadata_line.append(';'.join(unique_fastq_types))
        
        #3. Combine each metadata column of interest 
        for col in interest_metadata_columns:
            temp_combined_metadata = combine_metadata_rows(sample, temp_sample_metadata_df, col, no_warnings_metadata_columns)
            temp_sample_metadata_line.append(temp_combined_metadata)
        
        #4. Create new line in treated_metadata_df
        temp_metadata_line = pd.DataFrame([temp_sample_metadata_line], columns = DEFAULT_TREATED_METADATA_COL_NAMES + interest_metadata_columns)
        treated_metadata_df = pd.concat([treated_metadata_df, temp_metadata_line], ignore_index = True)


def treat_metadata(treatment_df, unique_samples_list, metadata_df, ena_download_column, interest_metadata_columns, no_warnings_metadata_columns, fastq_pattern, sep, n_sep, output_directory_path, color_treatment):
    """
    This function treats the provided Metadata Table based on the Treatment Template.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    unique_samples_list : list
        List of unique sample_names in Treatment Template.
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    ena_download_column : str
        Provided ENA Download Column from the Metadata Table.
    interest_metadata_columns : list
        List with the metadata columns of interest to treat.
    no_warnings_metadata_columns : list
        List of metadata column names that is normal/expected to have
        multiple values. No warning metadata columns.
    fastq_pattern : str
        Extension pattern to recognize a Fastq file.
    sep : str
        Sample Name separator.
    n_sep : int
        Sample Name separator appearance.
    output_directory_path : str
        The provided Output Directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.

    """
    #Set global variables
    global treated_metadata_df
    global warnings_df
    
    #A)Show main information pre-treatment
    ##Print information Title
    print(rich_text_colored('\nPre-treatment Information:\n', 'section_header', color_treatment))
    ##Show the number of samples
    print('  o', rich_text_colored('Number of Samples in Treatment Template:', 'subsection', color_treatment), len(unique_samples_list))
    ##Show the number of files
    print('  o', rich_text_colored('Number of Lines in the Original Metadata Table:', 'subsection', color_treatment), len(metadata_df))
    
    #B)Get METADATA for each sample_name and treat
    ##Section header
    print(rich_text_colored('\nTreat Metadata:', 'section_header', color_treatment))
    print('\nTreating metadata ...')
    ##Treat metadata per sample
    for sample in unique_samples_list:
        #Filter table by sample
        temp_sample = treatment_df[treatment_df['sample_name'] == sample]
        #Get treatment(get unique with set function)
        ##NOTE: Up to this point we have check there is only one treatment per sample
        temp_treatment = list(temp_sample['treatment'])[0]
        
        #Treat metadata depending on mode used for treating fastqs
        if temp_treatment == 'copy':
            copy_only_mode_metadata(temp_sample, sample, metadata_df, ena_download_column, interest_metadata_columns, no_warnings_metadata_columns, fastq_pattern, sep, n_sep)
        else:
            generic_mode_metadata(temp_sample, sample, metadata_df, ena_download_column, interest_metadata_columns, no_warnings_metadata_columns, fastq_pattern, sep, n_sep)
    
    #C)Show main information post-treatment
    ##Print information Title
    print(rich_text_colored('\nPost-treatment Information:\n', 'section_header', color_treatment))
    ##Show the number of files
    print('  o', rich_text_colored('Number of Lines in the Treated Metadata Table:', 'subsection', color_treatment), len(treated_metadata_df))


#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter metadata_table
    parser.add_argument(
            '-m','--metadata_table', 
            action = 'store',
            required = True,
            help = 'Metadata Table [Expected sep=TABS]. Indicate the path to the Metadata Table file.'
    )
    ##Parameter treatment_file
    parser.add_argument(
            '-t','--treatment_template', 
            action = 'store',
            required = True ,
            help = 'Treatment Template [Expected sep=TABS]. Indicate the path to the Treatment Template file. See Documentation for more information and format details.'
    )
    ##Parameter ena_download_column
    parser.add_argument(
            '-c','--ena_download_column', 
            action = 'store',
            choices = ENA_FASTQ_URLS_COLUMNS,
            required = False,
            default = 'fastq_ftp',
            help = 'ENA Download Column (Optional) [Default:fastq_ftp]. Indicate the ENA Metadata Table column that was used to download Fastq files.'
    )
    ##Parameter extra_no_warning_columns
    parser.add_argument(
            '-e','--extra_no_warning_columns', 
            nargs = '+',
            required = False,
            help = 'Extra No Warning Columns (Optional). Indicate the column names of the Metadata Table that can be safely merged without warning. Provide column names separated by spaces (If a column name has spaces, quote it).'
    )
    ##Parameter sample_name_sep
    parser.add_argument(
            '-sep','--sample_name_sep', 
            action = 'store',
            default = '_',
            required = False,
            help = 'Sample Name separator (Optional). Indicate sample name separator for "fastq_file_name" column in Treatment Template [Default="_"].'
    )
    ##Parameter sample_name_sep_appearance
    parser.add_argument(
            '-n_sep','--sample_name_sep_appearance', 
            action = 'store',
            required = False,
            type = int,
            default = 1,
            help = 'Sample Name separator appearance (Optional). Indicate by which appearance of the separator the file name can be divided in sample_name + rest [Default=1 appearance].'
    )
    ##Parameter fastq_pattern
    parser.add_argument(
            '-p','--fastq_pattern', 
            action = 'store',
            default = '.fastq.gz',
            required = False,
            help = 'Fastq File Pattern (Optional) [Default:".fastq.gz"]. Indicate the pattern to identify Fastq files.'
    )
    ##Parameter r1_patterns
    parser.add_argument(
            '-r1','--r1_pattern', 
            action = 'store',
            default = '_1.fastq.gz',
            required = False,
            help = 'R1 File Pattern (Optional) [Default:"_1.fastq.gz"]. Indicate the pattern to identify R1 PAIRED Fastq files.'
    )
    ##Parameter r2_patterns
    parser.add_argument(
            '-r2','--r2_pattern', 
            action = 'store',
            default = '_2.fastq.gz',
            required = False,
            help = 'R2 File Pattern (Optional) [Default:"_2.fastq.gz"]. Indicate the pattern to identify R2 PAIRED Fastq files.'
    )
    ##Parameter output_directory
    parser.add_argument(
            '-o','--output_directory', 
            action = 'store',
            required = False,
            help = 'Output Directory (Optional). Indicate the path to the Output Directory. Output files will be created in the current directory if not indicated.'
    )
    ##Parameter plain_text
    parser.add_argument(
            '-x', '--plain_text',
            action = 'store_true',
            required = False,
            help = 'Plain Text Mode (Optional). If indicated, it will enable Plain Text mode, and text will appear without colors.'
    )
    ##Parameter version
    parser.add_argument(
            '-v','--version',
            action = 'version',
            version = 'version {} ({})'.format(VERSION, DATE)
    )

    #Process arguments
    args = parser.parse_args()
    metadata_table_path = args.metadata_table
    ena_download_column = args.ena_download_column
    extra_no_warning_cols = args.extra_no_warning_columns
    treatment_file_path = args.treatment_template
    outputdir_path = args.output_directory
    sample_name_separator = args.sample_name_sep
    sep_appearance = args.sample_name_sep_appearance
    fastq_pattern = args.fastq_pattern
    r1_files_pattern = args.r1_pattern
    r2_files_pattern = args.r2_pattern
    plain_text_bool = args.plain_text
    
    #Show Program headers
    print('')
    program_header('#',  __program__, 64, 2, plain_text_bool)
    print('')
    
    #Show Program parameters
    print(rich_text_colored('Program Parameters:', 'section_header', plain_text_bool))
    print(tabulate(vars(args).items(), headers = ['Argument', 'Value'], tablefmt = 'simple_outline'))
    
    #Try/Except block
    try:
        #0)Initial steps and checks
        
        #Set global variables
        global treated_metadata_df
        global warnings_df
        
        #Define headers used and no_warnings_column
        ##Set variables
        headers_used = DEFAULT_ENA_NO_WARNING_COLUMNS + ENA_FASTQ_URLS_COLUMNS
        no_warning_columns = DEFAULT_ENA_NO_WARNING_COLUMNS
        ##If extra_no_warning_cols parameters are given add to headers to be used if they are not already present
        ## And also add to final no_warning_columns list
        if type(extra_no_warning_cols) == list:
            headers_used = list(set(DEFAULT_ENA_NO_WARNING_COLUMNS + ENA_FASTQ_URLS_COLUMNS + extra_no_warning_cols))
            no_warning_columns = list(set(DEFAULT_ENA_NO_WARNING_COLUMNS + extra_no_warning_cols))
        
        #Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #Check PAIRED Fastq patterns
        check_fastq_PAIRED_patterns(fastq_pattern, r1_files_pattern, r2_files_pattern)
        
        #1)Try to load treatment and metadata files
        
        #Section header message
        print(rich_text_colored('\nLoading Files:\n', 'section_header', plain_text_bool))
        
        #Try to load Treatment Template file as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('Treatment Template file:', 'general_text', plain_text_bool))
        print(treatment_file_path)
        ##Load metadata file as pandas df
        treatment_table = pd.read_csv(treatment_file_path, sep='\t')
        
        #Try to load Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('\nMetadata Table file:', 'general_text', plain_text_bool))
        print(metadata_table_path)
        ##Load metadata file as pandas df
        metadata = pd.read_csv(metadata_table_path, sep='\t')
        
        #2)Checks related to the metadata file

        #Do checks and messages
        check_headers_metadata = check_headers(headers_used, metadata)
        frase0_2 = 'Error! Some of the needed headers are not in the Metadata Table!\n Check your metadata file!\n'
        frase1_2 = '\nThe headers needed are:'
        treat_headers_check(headers_used, check_headers_metadata, metadata, frase0_2, frase1_2, plain_text_bool)
        
        #3)Checks related to the treatment file
        
        ##Check headers in file and messages
        check_headers_treatmentfile = check_headers(TEMPLATE_FINAL_COLUMNS, treatment_table)
        frase0_1 = 'Error! Some of the needed headers are not in the Treatment Template!\n Check your treatment file!\n'
        frase1_1 = '\nThe headers needed are:'
        treat_headers_check(TEMPLATE_FINAL_COLUMNS, check_headers_treatmentfile, treatment_table, frase0_1, frase1_1, plain_text_bool)
        
        ##Check that values in "fastq_type" column from Treatment Template are valid and treatment
        check_expected_fastq_type_values = check_values(VALID_FASTQ_TYPES, treatment_table['fastq_type'])
        frase0_2 = 'Error! Some of the values of the "fastq_type" column in the Treatment Template are not valid!\n Check your treatment file!\n'
        frase1_2 = '\nValid values are:'
        treat_values_check(VALID_FASTQ_TYPES, check_expected_fastq_type_values, frase0_2, frase1_2, plain_text_bool)
        
        ##Check that values in "treatment" column from Treatment Template are valid and treatment
        check_expected_treatment_values = check_values(VALID_TREATMENTS, treatment_table['treatment'])
        frase0_3 = 'Error! Some of the values of the "treatment" column in the Treatment Template are not valid!\n Check your treatment file!\n'
        frase1_3 = '\nValid values are:'
        treat_values_check(VALID_TREATMENTS, check_expected_treatment_values, frase0_3, frase1_3, plain_text_bool)
        
        ##Check that there are no NA values in "sample_name" column from Treatment Template
        frase4 = 'Error! Some of the values of the "sample_name" column in the Treatment Template are NAs!\n Check your treatment file!'
        check_na_in_pandas_dataframe(treatment_table['sample_name'], frase4)
        
        ##Check that there are no NA values in "fastq_file_name" column from Treatment Template
        frase5 = 'Error! Some of the values of the "fastq_file_name" column in the Treatment Template are NAs!\n Check your treatment file!'
        check_na_in_pandas_dataframe(treatment_table['fastq_file_name'], frase5)
        
        ##Check that there are not duplicate names in "fastq_file_name" column
        check_duplicates_in_fastq_names(treatment_table, 5)
        
        #4)Checks per file
        
        ##Check that all files in Treatment Template have matching file names and fastq types
        treat_check_fastq_name_type(treatment_table, fastq_pattern, r1_files_pattern, r2_files_pattern, 5)
        
        ##Check that all files in Treatment Templante appear in the Metadata Table
        check_treatment_fastqs_in_metadata(treatment_table, metadata, ena_download_column, plain_text_bool)
        
        #5)Checks per sample
        
        ##Get unique sample names and sort
        unique_sample_names = list(set(treatment_table['sample_name']))
        unique_sample_names.sort()
        
        ##Check mixed treatments per sample
        check_treatment_for_samples(treatment_table, unique_sample_names, 5)
        
        ##Check rename mode files configurations per sample
        check_rename_samples(treatment_table, unique_sample_names, plain_text_bool)
        
        ##Check merge mode files configurations per sample
        check_merge_samples(treatment_table, unique_sample_names, plain_text_bool)
        
        #6)Treat metadata
        
        ##Get list of colnames of interest in metadata table
        interest_metadata_columns = [i for i in list(metadata.columns) if i not in ENA_COLUMNS_IGNORED]
        
        ##Treat metadata
        treat_metadata(treatment_table, unique_sample_names, metadata, ena_download_column, interest_metadata_columns, no_warning_columns, fastq_pattern, sample_name_separator, sep_appearance, outputdir_path, plain_text_bool)
        
        ##Save generated treated_metadata_df
        
        #Previous steps
        ##Get outfile name
        outfile_name = 'treated_' + metadata_table_path.split('/')[-1]
        ##Treat output_directory parameter / Get full output file path
        outputfile_path = treat_output_directory_parameter_outfiles(outfile_name, outputdir_path)
        
        #Show saved file message
        print(rich_text_colored('\nSaving results in file:', 'general_text', plain_text_bool))
        print(outputfile_path)
        
        #Save template
        treated_metadata_df.to_csv(outputfile_path, header = True, index = False, sep = '\t')
        
        #7)Treat warnings
        if len(warnings_df) > 0:
            #Previous steps
            ##Get outfile name
            outreport_name = 'warning_report_4_treated_' + metadata_table_path.split('/')[-1]
            ##Treat output_directory parameter / Get full output file path
            outreport_path = treat_output_directory_parameter_outfiles(outreport_name, outputdir_path)
            
            #Show warning message
            print(rich_text_colored('\nSome of the samples showed warnings when trying to treat their metadata!', 'program_warning', plain_text_bool))
            print(rich_text_colored('\nThis could be due to:','due_to_header', plain_text_bool))
            print('- The presence of different data types in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of quality controlled Fastqs in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of different sequencing technologies', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- Multiple sequencer output Fastqs (runs/lanes) from the same sample', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- Technical replicates from the same sample', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- PAIRED files uploaded as SINGLE Fastq files', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print("- Authors' mishandle, upload errors or others", rich_text_colored('[Dangerous]', 'dangerous', plain_text_bool))
            print(rich_text_colored('\nYou should:', 'you_should_header', plain_text_bool))
            print('- Use the following warning report to manually confirm which is your case')
            print('- Check if in your particular case these metadata combination merges with several values are acceptable or not')
            print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')
            print('- Check the original publication and supplementary tables to get some context')
            print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
            
            #Show saved file message
            print(rich_text_colored('\nSaving report in file:', 'general_text', plain_text_bool))
            print(outreport_path)
            
            #Save template
            warnings_df.to_csv(outreport_path, header = True, index = False, sep = '\t')
            
            #Show advise legend
            show_advise_legend(plain_text_bool)
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
           
    finally:
        #Print empty line for aesthetic purposes
        print('')
                                                                                                   
if __name__ == '__main__':
    main()
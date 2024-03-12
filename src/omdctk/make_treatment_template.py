#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: make_treatment_template.py
@author: Samuel Piquer-Esteban
@date: 5 Mar 2024 

"""

#Set program name
__program__ = 'make_treatment_template.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, ENA_FASTQ_URLS_COLUMNS, TEMPLATE_FINAL_COLUMNS, VALID_MANIFEST_COLUMNS, 
                    OMD_CTK_Exception, program_header, check_existence_directory_parameter, 
                    check_fastq_PAIRED_patterns, rich_text_colored, check_na_in_pandas_dataframe,
                    get_list_fastqs_in_directory, print_list_n_byline,
                    check_headers, treat_headers_check, get_urls_from_ENA_column, 
                    metadata_files_main_information, treat_output_directory_parameter_outfiles)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
import pandas as pd

#Program Constants
DEFAULT_SAMPLE_COLUMNS = ['sample_accession', 'secondary_sample_accession',
                          'run_accession', 'library_name', 'run_alias', 
                          'sample_alias', 'sample_title']

#Program Functions
def check_colname_duplicates(possible_sample_names, color_treatment):
    """
    This function checks if the possible columns for sample names provided match
    the final template column names.

    Parameters
    ----------
    possible_sample_names : list
        List with the possible columns to use as sample name column.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If the possible columns provided match the final template columns names
        raises an exception.

    Returns
    -------
    None.

    """
    #Init results list
    column_duplicates = []
    
    #Check if the possible sample column names match the TEMPLATE_FINAL_COLUMNS
    for possible_name in possible_sample_names:
        if possible_name in TEMPLATE_FINAL_COLUMNS:
            column_duplicates.append(possible_name)
    
    #Raise exception if duplicates are detected
    if len(column_duplicates) > 0:
        #Prepare exception frase
        frase0 = 'Error! Some of the sample name columns provided match with the final template columns names!\n'
        frase1 = '\nPlease, change the following column names in the provided Metadata Table:'
        frase = [frase0, frase1]
        for i in column_duplicates:
            temp_line = ' '.join(['\n-', rich_text_colored(i, 'check_color', color_treatment)])
            frase.append(temp_line)
        #Raise exception
        raise OMD_CTK_Exception(' '.join(frase))


def check_fastqs_uniqueness(dataframe, column, fastqs_in_directory, color_treatment):
    """
    This function checks if the provided Fastqs have a unique match in the 
    Metadata Table.

    Parameters
    ----------
    dataframe : pandas dataframe
        The provided dataframe to use as reference.
        ENA Mode: ENA Metadata Table
        Generic Mode: Generic Manifest Table
    column : str
        The provided column to search the fastq files in reference dataframe.
        ENA Mode: ENA metadata column from which the urls were obtained.
        Generic Mode: "file_name" column in Manifest Table.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If there are Fastqs that do not have a unique entry in the reference dataframe 
        raises an exception.

    Returns
    -------
    None.

    """
    #Init results lists
    fastqs_warnings_more_than_one_match = []
    fastqs_warnings_zero_matches = []
    
    #Iter and process fastq files info
    for fastq in fastqs_in_directory:
        #Get temp fastq information dataframe
        temp_fastq_df = dataframe[dataframe[column].str.contains(fastq)]
        #Treat info depending on temp_fastq_df len()/ Treat warnings/ Else pass
        if len(temp_fastq_df) == 1:
            pass
        elif len(temp_fastq_df) > 1:
            fastqs_warnings_more_than_one_match.append(fastq)
        else:
            fastqs_warnings_zero_matches.append(fastq)
    
    #If there are warnings raise exception
    if len(fastqs_warnings_more_than_one_match) > 0 or len(fastqs_warnings_zero_matches) > 0:
        if len(fastqs_warnings_more_than_one_match) > 0:
            print(rich_text_colored('\nWarning! The following Fastq files in the provided directory have more than one match in the Input Reference File!', 'program_warning', color_treatment))
            print_list_n_byline(fastqs_warnings_more_than_one_match, 5)
        elif len(fastqs_warnings_zero_matches) > 0:
            print(rich_text_colored('\nWarning! The following Fastq files in the provided directory do not have a match in the Input Reference File!', 'program_warning', color_treatment))
            print_list_n_byline(fastqs_warnings_zero_matches, 5)
        raise OMD_CTK_Exception('Error! Some of the Fastqs in the provided directory do not have a unique match with the Input Reference File (zero or more than one match)!')    


def get_fastq_type(fastq_name, r1_pattern, r2_pattern):
    """
    This function gets the fastq_type from fastq_file_name.

    Parameters
    ----------
    fastq_name : str
        Fastq file name in fastqs directory.
    r1_pattern : str
        The provided pair1 pattern.
    r2_pattern : str
        The provided pair2 pattern.

    Returns
    -------
    result : str
        Returns the fastq_type [pair1, pair2 or single] 
        for the provided fastq file name.
    """   
    if fastq_name.endswith(r1_pattern):
        result = 'pair1'
    elif fastq_name.endswith(r2_pattern):
        result = 'pair2'
    else:
        result = 'single'
    return result


def contruct_template_from_ENA(metadata_df, column, fastqs_in_directory, possible_sample_names, r1_files_pattern, r2_files_pattern):
    """
    This function constructs the raw template dataframe when using the ENA mode.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    column : str
        The ENA metadata column from which the urls were obtained.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    possible_sample_names : list
        List with the possible columns to use as sample names.
    r1_files_pattern : str
        The provided R1 File Pattern.
    r2_files_pattern : str
        The provided R2 File Pattern.

    Returns
    -------
    template : pandas dataframe
        The resulting template dataframe.

    """
    #Init empty pandas dataframe
    template=pd.DataFrame()
    
    #Set template columns
    template_columns = ['sample_name'] + possible_sample_names + ['fastq_file_name', 'fastq_type', 'treatment']
    
    #Iter and process fastq files info
    for fastq in fastqs_in_directory:
        #Get temp fastq information dataframe
        temp_fastq_df = metadata_df[metadata_df[column].str.contains(fastq)]
        #Extract info
        temp_row = []
        ##Add empty sample_name
        temp_row.append('')
        ##Add values for each possible_sample_name
        for i in possible_sample_names:
            temp_row.append(temp_fastq_df[i].values[0])
        ##Add fastq_file_name
        temp_row.append(fastq)
        ##Add fastq_type using get_fastq_type()
        temp_row.append(get_fastq_type(fastq, r1_files_pattern, r2_files_pattern))
        ##Add empty treatment
        temp_row.append('')
        #Save info in dataframe
        temp_template = pd.DataFrame([temp_row], columns = template_columns)
        template = pd.concat([template, temp_template], ignore_index = True)
    
    #Return final template
    return template


def contruct_generic_template(manifest_df, fastqs_in_directory, r1_files_pattern, r2_files_pattern):
    """
    This function constructs the raw template dataframe when using the Generic mode.

    Parameters
    ----------
    manifest_df : pandas dataframe
        The provided manifest dataframe.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    r1_files_pattern : str
        The provided R1 File Pattern.
    r2_files_pattern : str
        The provided R2 File Pattern.

    Returns
    -------
    template : pandas dataframe
        The resulting template dataframe.

    """
    #Init empty pandas dataframe
    template=pd.DataFrame()
    
    #Set template columns
    template_columns = ['sample_name', 'fastq_file_name', 'fastq_type', 'treatment']
    
    #Iter and process fastq files info
    for fastq in fastqs_in_directory:
        #Get temp fastq information dataframe
        temp_fastq_df = manifest_df[manifest_df['file_name'].str.contains(fastq)]
        #Extract info
        temp_row = []
        ##Add sample_name from manifest
        temp_row.append(temp_fastq_df['sample_name'].values[0])
        ##Add fastq_file_name
        temp_row.append(fastq)
        ##Add fastq_type using get_fastq_type()
        temp_row.append(get_fastq_type(fastq, r1_files_pattern, r2_files_pattern))
        ##Add empty treatment
        temp_row.append('')
        #Save info in dataframe
        temp_template = pd.DataFrame([temp_row], columns = template_columns)
        template = pd.concat([template, temp_template], ignore_index = True)
    
    #Return final template
    return template

#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter input_file
    parser.add_argument(
            '-i','--input_file', 
            action = 'store',
            required = True,
            help = 'Input Reference File. Indicate the path to the Input Reference File with the information to create raw treatment template.'
    )
    ##Parameter fastqs_directory
    parser.add_argument(
            '-d','--fastqs_directory', 
            action = 'store',
            required = True,
            help = 'Fastqs Directory. Indicate the path to the Fastqs Directory.'
    )
    ##Parameter mode
    parser.add_argument(
            '-s','--mode', 
            action = 'store',
            choices = ['ENA','Generic'],
            required = False,
            default = 'ENA',
            help = 'Execution Mode (Optional) [Default:ENA]. Options: 1) ENA Metadata Table File [Expected sep=TABS] or 2) Generic Manifest Table File [Expected sep=TABS].'
    )
    ##Parameter ena_download_column
    parser.add_argument(
            '-c','--ena_download_column', 
            action = 'store',
            choices = ENA_FASTQ_URLS_COLUMNS,
            required = False,
            default = 'fastq_ftp',
            help = 'ENA Download Column (Optional) [Default:fastq_ftp]. Indicate the ENA Metadata Table column that was used to download Fastq files. This parameter will be skipped if Generic mode is used.'
    ) 
    ##Parameter fastq_pattern
    parser.add_argument(
            '-p','--fastq_pattern', 
            action = 'store',
            default = '.fastq.gz',
            required = False,
            help = 'Fastq File Pattern (Optional) [Default:".fastq.gz"]. Indicate the pattern to identify Fastq files.'
    )
    ##Parameter r1_pattern
    parser.add_argument(
            '-r1','--r1_pattern', 
            action = 'store',
            default = '_1.fastq.gz',
            required = False,
            help = 'R1 File Pattern (Optional) [Default:"_1.fastq.gz"]. Indicate the pattern to identify R1 PAIRED Fastq files.'
    )
    ##Parameter r2_pattern
    parser.add_argument(
            '-r2','--r2_pattern', 
            action = 'store',
            default = '_2.fastq.gz',
            required = False,
            help = 'R2 File Pattern (Optional) [Default:"_2.fastq.gz"]. Indicate the pattern to identify R2 PAIRED Fastq files.'
    )
    ##Parameter extra_sample_columns
    parser.add_argument(
            '-e','--extra_sample_columns', 
            nargs = '+',
            required = False,
            help = 'Extra Sample Columns (Optional). Indicate the column names for the extra possible sample names separated by spaces (If a column name has spaces, quote it). This parameter will be skipped if Generic mode is used.'
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
    input_file_path = args.input_file
    program_mode = args.mode
    ena_download_column = args.ena_download_column
    files_directory = args.fastqs_directory
    outputdir_path = args.output_directory
    extra_sample_names = args.extra_sample_columns
    fastq_pattern = args.fastq_pattern
    r1_files_pattern = args.r1_pattern
    r2_files_pattern = args.r2_pattern
    plain_text_bool = args.plain_text
    #Skip parameters depending on mode
    if program_mode == 'Generic':
        args.ena_download_column = None
        args.extra_sample_columns = None
        
    #Show Program headers
    print('')
    program_header('#',  __program__, 64, 2, plain_text_bool)
    print('')
    
    #Show Program parameters
    print(rich_text_colored('Program Parameters:', 'section_header', plain_text_bool))
    print(tabulate(vars(args).items(), headers = ['Argument', 'Value'], tablefmt = 'simple_outline'))
    
    #Try/Except block
    try:
        #0)Initial checks and previous common preparations
        
        #Check that provided files directory exist
        check_existence_directory_parameter(files_directory, 'Fastqs', '--fastqs_directory')
        
        #Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #Check PAIRED Fastq patterns
        check_fastq_PAIRED_patterns(fastq_pattern, r1_files_pattern, r2_files_pattern)
        
        #Try to get list of Fastq files in the fastqs directory and sort by name
        fastqs_in_directory = get_list_fastqs_in_directory(files_directory, fastq_pattern)
        fastqs_in_directory.sort()
        
        #Previous steps for output file(common to both modes)
        ##Get outfile name
        outfile_name = 'raw_treatment_template_' + input_file_path.split('/')[-1]
        ##Treat output_directory parameter / Get full output file path
        outputfile_path = treat_output_directory_parameter_outfiles(outfile_name, outputdir_path)
        
        #1)Initial Section header message(common to both modes)
        print(rich_text_colored('\nLoading File:\n', 'section_header', plain_text_bool))
        
        #Treat differently depending on program mode
        ##ENA Mode
        if program_mode == 'ENA':
            #1)Try to load metadata file
            
            #Try to load Metadata Table as pandas dataframe
            ##Show loading file message
            print(rich_text_colored('Metadata Table file:', 'general_text', plain_text_bool))
            print(input_file_path)
            ##Load metadata file as pandas df
            metadata = pd.read_csv(input_file_path, sep='\t')
            
            #2)Check headers used in metadata table and set possible_sample_names
            
            #Define headers used and possible_sample_names
            ##Set variables
            possible_sample_names = DEFAULT_SAMPLE_COLUMNS
            headers_used = DEFAULT_SAMPLE_COLUMNS + ENA_FASTQ_URLS_COLUMNS
            ##If extra_column parameters are given add to headers to be used if they are not already present
            if type(extra_sample_names) == list:
                possible_sample_names = list(set(DEFAULT_SAMPLE_COLUMNS + extra_sample_names))
                headers_used = list(set(DEFAULT_SAMPLE_COLUMNS + ENA_FASTQ_URLS_COLUMNS + extra_sample_names))

            #Do checks and messages
            check_headers_metadata = check_headers(headers_used, metadata)
            frase0_2 = 'Error! Some of the needed headers are not in the Metadata Table!\n Check your metadata file!\n'
            frase1_2 = '\nThe headers needed are:'
            treat_headers_check(headers_used, check_headers_metadata, metadata, frase0_2, frase1_2, plain_text_bool)
            
            #3)Check that the possible_sample_names do not have a match with the final template column names 
            check_colname_duplicates(possible_sample_names, plain_text_bool)
            
            #4)Try to get urls from the ena_download_column
            urls = get_urls_from_ENA_column(metadata, ena_download_column)
            
            #5)Check file uniqueness in metadata file
            check_fastqs_uniqueness(metadata, ena_download_column, fastqs_in_directory, plain_text_bool)
            
            #6)Show main info
            #Section header message
            print(rich_text_colored('\nMain Information:', 'section_header', plain_text_bool))
            #Print main information
            metadata_files_main_information(metadata, urls, fastqs_in_directory, plain_text_bool)
            
            #7)Construct raw treatmentfile template
            #Section header message and get template
            print(rich_text_colored('\nCreating Raw Treatment Template:', 'section_header', plain_text_bool))
            template = contruct_template_from_ENA(metadata, ena_download_column, fastqs_in_directory, possible_sample_names, r1_files_pattern, r2_files_pattern)
        
        ##Generic Mode
        else:
            #1)Try to load Manifest Table as pandas dataframe
            ##Show loading file message
            print(rich_text_colored('Manifest Table file:', 'general_text', plain_text_bool))
            print(input_file_path)
            ##Load manifest file as pandas df
            manifest = pd.read_csv(input_file_path, sep='\t')
            
            #2)Check headers used
            ##Check headers in manifest table
            check_headers_manifest = check_headers(VALID_MANIFEST_COLUMNS, manifest)
            frase0_1 = 'Error! Some of the needed headers are not in the Manifest Table!\n Check your manifest file!\n'
            frase1_1 = '\nThe headers needed are:'
            treat_headers_check(VALID_MANIFEST_COLUMNS, check_headers_manifest, manifest, frase0_1, frase1_1, plain_text_bool)
            
            ##Check NAs in 'file_name' column manifest
            frase2 = 'Error! Some of the values of the "file_name" column in the Manifest Table are NAs!\n Check your manifest file!'
            check_na_in_pandas_dataframe(manifest['file_name'], frase2)
            
            ##Check NAs in 'file_name' column manifest
            frase3 = 'Error! Some of the values of the "sample_name" column in the Manifest Table are NAs!\n Check your manifest file!'
            check_na_in_pandas_dataframe(manifest['sample_name'], frase3)
            
            #3)Check uniqueness in  manifest file
            check_fastqs_uniqueness(manifest, 'file_name', fastqs_in_directory, plain_text_bool)
            
            #4)Show main information
            ##Section header message
            print(rich_text_colored('\nMain Information:', 'section_header', plain_text_bool))
            ##Show the number of fastq files in manifest
            print(rich_text_colored('\n1) Number of Fastq files in the provided Manifest Table:', 'subsection', plain_text_bool), len(manifest.index))
            ##Show the number of unique sample_accessions
            print(rich_text_colored('\n2) Number of unique samples for "sample_name" column in Manifest Table:', 'subsection', plain_text_bool), manifest['sample_name'].nunique(dropna = False))
            ##Show the number of files downloaded ending in 'fastq.gz' for files_directory
            print(rich_text_colored('\n3) Number of Fastq files in the provided Fastqs Directory:', 'subsection', plain_text_bool), len(fastqs_in_directory))

            #5)Construct raw treatmentfile template
            #Section header message and get template
            print(rich_text_colored('\nCreating Raw Treatment Template:', 'section_header', plain_text_bool))
            template = contruct_generic_template(manifest, fastqs_in_directory, r1_files_pattern, r2_files_pattern)
        
        #Final Common part
        ##Show saved file message and save result(common to both modes)
        print(rich_text_colored('\nSaving results in file:', 'general_text', plain_text_bool))
        print(outputfile_path)
        #Save template
        template.to_csv(outputfile_path, header = True, index = False, sep = '\t')
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
           
    finally:
        #Print empty line for aesthetic purposes
        print('')
        
if __name__ == '__main__':
    main()    
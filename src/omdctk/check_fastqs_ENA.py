#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: check_fastqs_ENA.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'check_fastqs_ENA.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, ENA_FASTQ_URLS_COLUMNS, program_header, 
                    show_advise_legend, get_urls_from_ENA_column, 
                    print_list_n_byline, metadata_files_main_information,
                    get_list_fastqs_in_directory, rich_text_colored,
                    check_existence_directory_parameter, 
                    check_headers, treat_headers_check)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
from hashlib import md5
import pandas as pd
import os

#Program Constants
HEADERS_USED = ['sample_accession', 'sample_alias', 'library_layout',
                'run_accession', 'fastq_ftp', 'fastq_aspera', 'fastq_galaxy',
                'submitted_ftp', 'submitted_aspera', 'submitted_galaxy',
                'fastq_md5', 'submitted_md5']

#Program Functions
def check_files_in_directory_frase(file_name, file_in_directory, color_treatment):
    """
    This function gets the check_files_in_directory_frase.

    Parameters
    ----------
    file_name : str
        The name of the file.
    file_in_directory : bool
        Bool value to see if the file is in the provided directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool_value : bool
        Bool to treat warnings.
    str_frase : str
        String with the check_files_in_directory_frase.

    """
    if file_in_directory == True:
        bool_value = True
        str_frase = file_name + ' ' + rich_text_colored('True! The file is in the provided directory!', 'true_color', color_treatment)
    else:
        bool_value = False
        str_frase = file_name + ' ' + rich_text_colored('False! The file is not in the provided directory!', 'false_color', color_treatment)
    return bool_value, str_frase

def treat_metadata_on_directory_checks(urls, fastqs_in_directory, color_treatment):
    """
    This function treats the metadata's fastqs checks on the fastqs directory.

    Parameters
    ----------
    urls : list
        List with the fastq urls from the Metadata file.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool
        True: There are warnings.
        False: There are no warnings.

    """
    #Init list results
    directory_checks_warnings = []
    
    #For each url do check
    for i in range(len(urls)):
        ##Get file name
        temp_file_name = urls[i].split('/')[-1]
        ##Check that the file is in fastqs_directory
        file_in_directory = temp_file_name in fastqs_in_directory
        ##Directory_check
        temp_directory_check = check_files_in_directory_frase(temp_file_name, file_in_directory, color_treatment)
        ##If has warning keep file name
        if temp_directory_check[0] == False:
            directory_checks_warnings.append(temp_file_name)
    
    #Treat directory_checks_warnings
    print(rich_text_colored('\n1) Check that expected files from the Metadata Table exist in the provided directory:', 'subsection', color_treatment))
    if len(directory_checks_warnings) == 0:
        print(rich_text_colored('All the expected Fastq files from the Metadata Table are in the provided directory!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        print(rich_text_colored('\nWarning! The following expected Fastq files from the Metadata Table are not in the provided directory!', 'program_warning', color_treatment))
        print_list_n_byline(directory_checks_warnings, 5)
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:','due_to_header', color_treatment))
        print('- The presence of untreated download errors', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm that this is your case')
        print('- If necessary, try to re-download affected files')
        #Return info for lately treat advise messages
        return True
    

def treat_directory_on_metadata_check(urls, fastqs_in_directory, color_treatment):
    """
    This function treats the fastqs directory checks on the metadata's fastqs.

    Parameters
    ----------
    urls : list
        List with the fastq urls from the Metadata file.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool
        True: There are warnings.
        False: There are no warnings.

    """
    ##Init list results
    fastqs_in_metadata = []
    ##For each url check
    for i in range(len(urls)):
        ##Get file name
        temp_file_name = urls[i].split('/')[-1]
        ##Save file name for latter
        fastqs_in_metadata.append(temp_file_name)
    
    #Check fastqs that are in the provided directory but not in the metadata table
    print(rich_text_colored('\n2) Check if there are Fastq files of the provided directory absent in the Metadata Table:', 'subsection', color_treatment))
    ##Get asymmetric difference between list
    fastqs_difference = list(set(fastqs_in_directory) - set(fastqs_in_metadata))
    ##If fastqs_difference is empty show message/ Else show different files
    if len(fastqs_difference) == 0:
        print(rich_text_colored('There are no extra Fastq files in the provided directory!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        print(rich_text_colored('\nWarning! The following Fastq files are not in the Metadata Table!', 'program_warning', color_treatment))
        print_list_n_byline(fastqs_difference, 5)
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:','due_to_header', color_treatment))
        print('- The mixture of different dataset projects', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm that this is your case and act accordingly')
        #Return info for lately treat advise messages
        return True


def treat_multiple_fastq_matches(metadata_df, column, fastqs_in_directory, color_treatment):
    """
    This function treats if fastqs in the provided directory have multiple
    run_accession matches in the metadata.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    column : str
        The ENA metadata column from which the urls were obtained.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool
        True: There are warnings.
        False: There are no warnings.

    """
    #Init results lists
    fastqs_warnings_multiple_match = []
    
    #Iter and process fastq files info
    for fastq in fastqs_in_directory:
        #Get temp fastq information dataframe
        temp_fastq_df = metadata_df[metadata_df[column].str.contains(fastq)]
        #Treat info depending on temp_fastq_df len()/ Treat warnings/ Else pass
        if len(temp_fastq_df) > 1:
            fastqs_warnings_multiple_match.append(fastq)
        else:
            pass
        
    #Check fastqs that are in the provided directory but have multiple matches in the metadata table
    print(rich_text_colored('\n3) Check if there are Fastq files of the provided directory with multiple matches in the Metadata Table:', 'subsection', color_treatment))
    
    ##If fastqs_warnings_multiple_match is empty show message/ Else show different files
    if len(fastqs_warnings_multiple_match) == 0:
        print(rich_text_colored('All Fastq files in the provided directory have a unique run_accession match with the Metadata Table!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        print(rich_text_colored('\nWarning! The following Fastq files have multiple run_accession matches with the Metadata Table!', 'program_warning', color_treatment))
        print_list_n_byline(fastqs_warnings_multiple_match, 5)
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:','due_to_header', color_treatment))
        print('- The presence of duplicated entries in the Metadata Table', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm that this is your case and act accordingly')
        #Return info for lately treat advise messages
        return True
    

def get_md5_from_ENA_column(metadata_df, column, file_name):
    """
    This function gets the md5s for the file in the indicated ENA metadata column.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    column : str
        The ENA metadata column from which the urls were obtained.
    file_name : str
        The Fastq file name

    Returns
    -------
    md5_sums : list
        The list of md5s for the seleccted ENA metadata column.
        
    """
    #Try/except/
    try:
        #Get index for row that contains file_name
        row_index = metadata_df[metadata_df[column].str.contains(file_name)].index.values.item()
        #Get sample fastq urls
        file_match = metadata_df.loc[row_index][column].split(';')
        #Get inner index
        inner_index = [i for i, w in enumerate(file_match) if file_name in w][0]
        #Get md5_sum for fastq file
        if column in ['fastq_ftp', 'fastq_aspera', 'fastq_galaxy']:
            md5_sums = metadata_df.loc[row_index]['fastq_md5'].split(';')[inner_index]
        else:
            md5_sums = metadata_df.loc[row_index]['submitted_md5'].split(';')[inner_index]
    except:
        md5_sums = 'metadata_error'
    return md5_sums


def get_md5_checksum_for_file(path_file):
    """
    This function gets MD5 checksum of a file.

    Parameters
    ----------
    path_file : str
        Path to the file to process.

    Returns
    -------
    digest : str
        The MD5 checksum of the provided file.
        
    """
    try:
        #Read file content
        file = open(path_file, 'rb')
        content = file.read()
        file.close()
    
        #Get md5 checksum
        md5_hash = md5()
        md5_hash.update(content)
        digest = md5_hash.hexdigest()
    except:
        digest = 'file_error'
    return digest


def treat_optional_md5(urls, fastqs_in_directory, metadata_df, column, directory, color_treatment):
    """
    This function treats the optional md5 check.

    Parameters
    ----------
    urls : list
        List with the fastq urls from the Metadata file.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    column : str
        The ENA metadata column from which the urls were obtained.
    directory : str
        Path to the provided directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool
        True: There are warnings.
        False: There are no warnings.

    """
    #Show md5_checks message
    print(rich_text_colored('\n4) Check that expected MD5s from the Metadata Table match calculated MD5s:', 'subsection', color_treatment))
    
    #Get md5_column based on ENA Download column provided
    if column in ['fastq_ftp', 'fastq_aspera', 'fastq_galaxy']:
        md5_column = 'fastq_md5'
    else:
        md5_column = 'submitted_md5'
    
    #Check that the md5_column is not all empty
    if all(metadata_df[md5_column].isna()):
        print(rich_text_colored('\nNone of the MD5 checksums are available for the provided ENA Download Column!', 'program_warning', color_treatment))
        print(rich_text_colored("Skipping MD5s' check, since there is nothig to compare with!", 'program_warning', color_treatment))
        return False
    else:
        #Print message
        print(rich_text_colored('\nThis may take a while...', 'program_warning2', color_treatment))
        
        #Init results lists
        warning_file_not_in_dir = []
        warning_file_md5_not_match = []
        warning_metadata_md5_not_found = []
        warning_md5_file_error = []
        
        ##For each url try to check md5
        for i in range(len(urls)):
            ##Get file name
            temp_file_name = urls[i].split('/')[-1]
            ##Get MD5 checksum from metadata
            temp_metadata_md5 = get_md5_from_ENA_column(metadata_df, column, temp_file_name)
            ##Check that MD5 checksum was found
            if temp_metadata_md5 == 'metadata_error':
                warning_metadata_md5_not_found.append(temp_file_name)
            else:
                ##Check that the file is in fastqs_directory
                file_in_directory = temp_file_name in fastqs_in_directory
                if file_in_directory == False:
                    warning_file_not_in_dir.append(temp_file_name)
                else:
                    #Get full path to file
                    file_path = os.path.join(directory, temp_file_name)
                    #Calcule md5 checksum
                    calculated_md5 = get_md5_checksum_for_file(file_path)
                    #Check if the md5s match or if there was a file_error
                    if calculated_md5 == temp_metadata_md5:
                        pass
                    elif calculated_md5 == 'file_error':
                        warning_md5_file_error.append(temp_file_name)
                    else:
                        warning_file_md5_not_match.append(temp_file_name)
        
        #Calculate total number of warning runs
        n_file_warnings = len(warning_file_not_in_dir + warning_file_md5_not_match + warning_metadata_md5_not_found + warning_md5_file_error)
        
        #Show total number of runs again and total warnings
        print('  o', rich_text_colored('Total number of Fastqs detected in Metadata Table:', 'general_text', color_treatment), len(urls))
        print('  o', rich_text_colored('Total number of warnings for MD5s :', 'general_text', color_treatment), n_file_warnings)
        
        #Show different information depending on result 
        if n_file_warnings == 0:
            print(rich_text_colored('\nAll MD5s present in the Metadata Table match the MD5s calculated for their corresponding Fastq files!', 'acceptable', color_treatment))
            #Return info for lately treat advise messages
            return False
        else:
            ##Show warning types stats
            print(rich_text_colored('\nMD5 Warning Types:', 'subsection2', color_treatment))
            print('  o', rich_text_colored('Number of Fastqs with mismatches between metadata and calculated MD5s:','general_text', color_treatment), len(warning_file_md5_not_match))
            print('  o', rich_text_colored('Number of Fastqs with errors when calculating MD5s:', 'general_text', color_treatment), len(warning_md5_file_error))
            print('  o', rich_text_colored('Number of Fastqs not found in the provided directory:', 'general_text', color_treatment), len(warning_file_not_in_dir))
            print('  o', rich_text_colored('Number of Fastqs with MD5s not found in the Metadata Table:', 'general_text', color_treatment), len(warning_metadata_md5_not_found))
            
            ##Show warning message
            print(rich_text_colored('\nWarning! Some MD5s in the Metadata Table do not match the MD5s calculated for their corresponding Fastq files!', 'program_warning', color_treatment))
            
            #Show run_accession values if there are warnings
            ##warning_file_md5_not_match
            if len(warning_file_md5_not_match) > 0:
                print('\n-', rich_text_colored('Fastqs with mismatches between metadata and calculated MD5s:', 'general_text', color_treatment))
                print_list_n_byline(warning_file_md5_not_match, 5)
            ##warning_md5_file_error
            if len(warning_md5_file_error) > 0:
                print('\n-', rich_text_colored('Fastqs with errors when calculating MD5s:', 'general_text', color_treatment))
                print_list_n_byline(warning_md5_file_error, 5)
            ##warning_file_not_in_dir
            if len(warning_file_not_in_dir) > 0:
                print('\n-', rich_text_colored('Fastqs not found in the provided directory:', 'general_text', color_treatment))
                print_list_n_byline(warning_file_not_in_dir, 5)
            ##warning_metadata_md5_not_found
            if len(warning_metadata_md5_not_found) > 0:
                print('\n-', rich_text_colored('Fastqs with MD5s not found in the Metadata Table:', 'general_text', color_treatment))
                print_list_n_byline(warning_metadata_md5_not_found, 5)
            
            #Show advise messages
            print(rich_text_colored('\nThis could be due to:', 'due_to_header', color_treatment))
            print('- Memory errors when calculating MD5s', rich_text_colored('[Acceptable]', 'acceptable', color_treatment))
            print('- The presence of untreated download errors', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
            print("- Upload errors", rich_text_colored('[Dangerous]', 'dangerous', color_treatment))
            print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
            print('- Manually confirm which is your case')
            print('- If any mismatches, try to re-download affected files')
            print('- If any errors calculating MD5s, try to provide more memory')
            print('- If any, check missing Fastq cases')
            print('- If any, check missing MD5 cases')
            print('- If necessary, contact the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
            
            #Return info for lately treat advise messages
            return True
        
            
#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter metadata_table
    parser.add_argument(
            '-t','--metadata_table', 
            action = 'store',
            required = True,
            help = 'Metadata Table [Expected sep=TABS]. Indicate the path to the Metadata Table file.'
    )
    ##Parameter fastqs_directory
    parser.add_argument(
            '-d','--fastqs_directory', 
            action = 'store',
            required = True,
            help = 'Fastqs Directory. Indicate the path to the Fastqs Directory.'
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
    ##Parameter fastq_pattern
    parser.add_argument(
            '-p','--fastq_pattern', 
            action = 'store',
            default = '.fastq.gz',
            required = False,
            help = 'Fastq File Pattern (Optional) [Default:".fastq.gz"]. Indicate the pattern to identify Fastq files.'
    )
    ##Parameter md5_check
    parser.add_argument(
            '-m', '--md5_check',
            action = 'store_true',
            required = False,
            help = 'MD5 Check (Optional). If indicated, it will enable MD5 Check mode.'
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
    files_directory = args.fastqs_directory
    fastq_pattern = args.fastq_pattern
    md5_option = args.md5_check
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
        #0)Initial checks
        #Check that provided files directory exist
        check_existence_directory_parameter(files_directory, 'Fastqs', '--fastqs_directory')
        
        #Try to get list of fastqs in the provided directory
        fastqs_in_directory = get_list_fastqs_in_directory(files_directory, fastq_pattern)
        
        #1)Try to load metadata file
        
        #Section header message
        print(rich_text_colored('\nLoading File:\n', 'section_header', plain_text_bool))
        
        #Try to load Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('Metadata Table file:', 'general_text', plain_text_bool))
        print(metadata_table_path)
        ##Load metadata file as pandas df
        metadata = pd.read_csv(metadata_table_path, sep='\t')
        
        #2)Check headers used in metadata table
        check_headers_metadata = check_headers(HEADERS_USED, metadata)
        frase0_2 = 'Error! Some of the needed headers are not in the Metadata Table!\n Check your metadata file!\n'
        frase1_2 = '\nThe headers needed are:'
        treat_headers_check(HEADERS_USED, check_headers_metadata, metadata, frase0_2, frase1_2, plain_text_bool)
        
        #3)Try to get urls from the ena_download_column
        urls = get_urls_from_ENA_column(metadata, ena_download_column)
        
        #4)Show Main information
        
        #Section header message
        print(rich_text_colored('\nMain Information:', 'section_header', plain_text_bool))
        
        #Print main information
        metadata_files_main_information(metadata, urls, fastqs_in_directory, plain_text_bool)
        
        #5)Fastq Checks 
        
        #Section header
        print(rich_text_colored("\nFastqs' Checks:", 'section_header', plain_text_bool))
        
        #5.1)Check fastqs in metadata
        check1 = treat_metadata_on_directory_checks(urls, fastqs_in_directory, plain_text_bool)        
        
        #5.2)Check fastqs in directory
        check2 = treat_directory_on_metadata_check(urls, fastqs_in_directory, plain_text_bool)
        
        #5.3)Check multuple fastq matches in metadata
        check3 = treat_multiple_fastq_matches(metadata, ena_download_column, fastqs_in_directory, plain_text_bool)
        
        #5.4)Treat optional md5_option 
        if md5_option == True:
            check4 = treat_optional_md5(urls, fastqs_in_directory, metadata, ena_download_column, files_directory, plain_text_bool)
            functions_warnings = (check1, check2, check3, check4)
        else:
            functions_warnings = (check1, check2, check3)
            
        #6)Show Legend if warnings were detected
        if any(functions_warnings) == True:
            show_advise_legend(plain_text_bool)
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)

    except MemoryError:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), 'MemoryError')
    
    finally:
        #Print empty line for aesthetic purposes
        print('')

if __name__ == '__main__':
    main()
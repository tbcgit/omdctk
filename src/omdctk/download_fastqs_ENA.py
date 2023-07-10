#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: download_fastqs_ENA.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'download_fastqs_ENA.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, ENA_FASTQ_URLS_COLUMNS, OMD_CTK_Exception, 
                    program_header, get_urls_from_ENA_column, rich_text_colored,
                    check_existence_directory_parameter, 
                    check_headers, treat_headers_check,
                    treat_output_directory_parameter,
                    treat_output_directory_parameter_outfiles)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
from itertools import compress
import pandas as pd
import parfive
import os

#Program Constants
METADATA_HEADERS = ['sample_accession', 'sample_alias', 'library_layout',
                    'run_accession', 'fastq_ftp', 'fastq_aspera', 'fastq_galaxy',
                    'submitted_ftp', 'submitted_aspera', 'submitted_galaxy']

#Program Global Variables
urls = []
errors_df = pd.DataFrame()

#Program Functions
def LINKS_mode(inputfile_path, color_treatment):
    """
    This function treats the LINKS mode. Loads the input file, get lines of 
    the provided txt file (download urls) and shows main information.

    Parameters
    ----------
    inputfile_path : str
        Path to the input file.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If all the lines of the inputfile are empty raises an exception.

    Returns
    -------
    urls : list
        List with the urls to be downloaded.

    """
    #Section header message
    print(rich_text_colored('\nLoading Input File:\n', 'section_header', color_treatment))
    
    #Load Metadata Table as pandas dataframe
    ##Show loading file message
    print(rich_text_colored('Links TXT file:', 'general_text', color_treatment))
    print(inputfile_path)
    
    #Load links file
    txt_file = open(inputfile_path, "r")
    urls = txt_file.read().splitlines()
    txt_file.close()
    
    #Process links
    ##Remove empty values in list
    urls = list(filter(None, urls))
    
    #Check if there are empty lines
    if len(urls) == 0:
        raise OMD_CTK_Exception('Error! All lines are empty in the provided Links TXT file!')
    
    #Show info
    ##Print information Title
    print(rich_text_colored('\nMain Information:\n', 'section_header', color_treatment))
    ##Show the number of detected URLS to download
    print(rich_text_colored('Number of URLs to download:', 'general_text', color_treatment), len(urls))
    
    #Return urls
    return urls
    
    
def ENA_mode(inputfile_path, ena_download_column, color_treatment):
    """
    This function treats the ENA mode. Loads the input file as
    pandas dataframe, checks that the needed headers are present, gets download
    urls from the ena_download_column and shows main information.

    Parameters
    ----------
    inputfile_path : str
        Path to the input file.
    ena_download_column : str
        The ENA Download Column.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    urls : list
        List with the urls to be downloaded.

    """
    #Section header message
    print(rich_text_colored('\nLoading Input File:\n', 'section_header', color_treatment))
    
    #Load Metadata Table as pandas dataframe
    ##Show loading file message
    print(rich_text_colored('Metadata Table file:', 'general_text', color_treatment))
    print(inputfile_path)
    ##Try to load the ena table file
    metadata = pd.read_csv(inputfile_path, sep = '\t')
    
    #Check headers and treatment
    check_headers_metadata = check_headers(METADATA_HEADERS, metadata)
    frase0_2 = 'Error! Some of the needed headers are not in the Metadata Table!\n Check your metadata file!\n'
    frase1_2 = '\nThe headers needed are:'
    treat_headers_check(METADATA_HEADERS, check_headers_metadata, metadata, frase0_2, frase1_2, color_treatment)
    
    #Get list of links to download from ENA Download column/ If empty raises exception 
    urls = get_urls_from_ENA_column(metadata, ena_download_column)
    
    #Show stats info
    ##Print information Title
    print(rich_text_colored('\nMain Information:', 'section_header', color_treatment))
    ##Show the number of run_accessions
    print(rich_text_colored('\n1) Number of run_accessions:', 'subsection', color_treatment), len(metadata.index))
    ##Show the number of unique sample_accessions
    print(rich_text_colored('\n2) Number of unique sample_accessions:', 'subsection', color_treatment), metadata['sample_accession'].nunique(dropna = False))
    ##Show the number of unique sample_alias
    print(rich_text_colored('\n3) Number of unique sample_alias:', 'subsection', color_treatment), metadata['sample_alias'].nunique(dropna = False))
    ##Show the number of appereances per library_layout
    print(rich_text_colored('\n4) Appearances per library_layout:', 'subsection', color_treatment))
    print(metadata['library_layout'].value_counts(dropna = False).rename_axis('unique_values').reset_index(name = 'counts').to_string(index = False))
    ##Show the number of detected URLS to download
    print(rich_text_colored('\n5) Number of URLs to download:', 'subsection', color_treatment), len(urls))
    
    #Return urls
    return urls
    

def check_urls_domain_and_treatment(inputfile_type):
    """
    This function checks the urls' domain. It will update the urls list and
    the errors_df if needed depending on the inputfile_type.

    Parameters
    ----------
    inputfile_type : str
        The execution mode ['ENA','LINKS'].

    Raises
    ------
    OMD_CTK_Exception
        If all urls lack a valid domain in the inputfile_type == 'LINKS' 
        raises an exception.

    Returns
    -------
    None.

    """
    #Set global variables
    global urls
    global errors_df
    
    #Check url protocol/ if it does not start by http:// , https:// or ftp://
    url_checks = [i.startswith('http://') or i.startswith('https://') or i.startswith('ftp://') for i in urls]
    
    #If all values are False and inputfile_type == 'LINKS' raise exception
    if (not any(url_checks)) and inputfile_type == 'LINKS':
        raise OMD_CTK_Exception('Error! None of the provided Links start with a valid domain (http://, https://, ftp://)!')
    #Elif all values are False and inputfile_type == 'ENA' modify urls
    #Urls from ENA metadata file lack ftp:// domain
    elif (not any(url_checks)) and inputfile_type == 'ENA':   
        #Add ftp domain
        urls = ['ftp://'+i for i in urls]
    else:
        #Get invalid URLs
        reverse_url_checks = [not i for i in url_checks]
        notvalid_urls = list(compress(urls, reverse_url_checks))
        #Get valid URLs to download
        urls = list(compress(urls, url_checks))
        #Add invalid URLSs to errors dataframe
        for notvalid_url in notvalid_urls:
            temp_error_message = 'URL does not start with a valid domain (http://, https://, ftp://)'
            temp_error = pd.DataFrame([[notvalid_url, temp_error_message]], columns = ['URL', 'Error'])
            errors_df = pd.concat([errors_df, temp_error], ignore_index = True)  


def download_files_with_parfive(urls, parfive_verb, n_max_conn, outputfiles_path, color_treatment):
    """
    This function downloads the provided urls using parfive and updates errors_df
    if needed.

    Parameters
    ----------
    urls : list
        List of urls to download.
    parfive_verb : bool
        Input parameter to treat parfive's verbose mode.
    n_max_conn : int
        Input parameter to indicate number of files to download in parallel.
    outputfiles_path : str
        Output files directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.

    """
    #0)Set global variables
    global errors_df
    
    #1)Section header message 
    print(rich_text_colored('\nDownloading Files (parfive):\n', 'section_header', color_treatment))
    #Show downloading file message
    print(rich_text_colored('Downloader version:', 'general_text', color_treatment))
    ##Get parfive version
    print('parfive')
    print(parfive.__version__)
    print('')
    
    #2)Set env variables/ 
    ##Default timeout to 5 min (300s) 
    os.environ["PARFIVE_SOCK_READ_TIMEOUT"] = "300"
    ##Treat Parfive verbose mode
    if parfive_verb == False:
        os.environ["PARFIVE_HIDE_PROGRESS"] = "True"
        save_file_frase = 'Resulting files saved in:'
    else:
        save_file_frase = '\nResulting files saved in:'
    
    #3)Download files
    ##Init downloader object
    downloader = parfive.Downloader(max_conn = n_max_conn)
    ##Append urls to download queue
    for url in urls: 
        downloader.enqueue_file(url, path = outputfiles_path)
    ##Start Download
    files = downloader.download()
    
    #4)Add errors to erros dataframe if there are errors
    if len(files.errors) != 0:
        for error in files.errors:
            tempURL = error[1]
            temp_error_message = error[2]
            temp_error = pd.DataFrame([[tempURL, temp_error_message]], columns = ['URL', 'Error'])
            errors_df = pd.concat([errors_df, temp_error], ignore_index = True)                    
    
    #5)Files saved message/ if there was any file downloaded
    if len(urls) > len(files.errors):
        #Show saved file message
        print(rich_text_colored(save_file_frase, 'general_text', color_treatment))
        print(outputfiles_path)
    
    
#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter input_file
    parser.add_argument(
            '-i','--input_file', 
            action = 'store',
            required = True,
            help = 'Input File. Indicate the path to the Input File with the information to download Fastq files.'
    )
    ##Parameter mode
    parser.add_argument(
            '-m','--mode', 
            action = 'store',
            choices = ['ENA','LINKS'],
            required = False,
            default = 'ENA',
            help = 'Execution Mode (Optional) [Default:ENA]. Options: 1) ENA Metadata Table File [Expected sep=TABS] or 2) Links TXT File.'
    )
    ##Parameter ena_download_column
    parser.add_argument(
            '-c','--ena_download_column', 
            action = 'store',
            choices = ENA_FASTQ_URLS_COLUMNS,
            required = False,
            default = 'fastq_ftp',
            help = 'ENA Download Column (Optional) [Default:fastq_ftp]. Indicate the ENA Metadata Table column with the links to download. This parameter will be skipped if LINKS mode is used.'
    )
    ##Parameter max_conn
    parser.add_argument(
            '-n','--max_conn', 
            action = 'store',
            required = False,
            type = int,
            default = 5,
            help = 'Max Number of Files in Parallel (Optional:Parfive Parameter) [Default:5]. Indicate the max number of files to be downloaded simultaneously in parallel.'
    )
    ##Parameter parfive_verbose
    parser.add_argument(
            '-p', '--parfive_verbose',
            action = 'store_true',
            required = False,
            default = False,
            help = 'Parfive Verbose (Optional:Parfive Parameter). If indicated it will enable parfive verbose.'
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
    inputfile_path = args.input_file
    inputfile_type = args.mode
    outputdir_path = args.output_directory
    ena_download_column = args.ena_download_column
    n_max_conn = args.max_conn
    parfive_verb = args.parfive_verbose
    plain_text_bool = args.plain_text
    #Skip ena_download_column if LINKS mode is used
    if inputfile_type == 'LINKS':
        args.ENA_download_column = None
    
    #Show Program headers
    print('')
    program_header('#',  __program__, 64, 2, plain_text_bool)
    print('')
    
    #Show Program parameters
    print(rich_text_colored('Program Parameters:', 'section_header', plain_text_bool))
    print(tabulate(vars(args).items(), headers = ['Argument', 'Value'], tablefmt = 'simple_outline'))
    
    #Try/Except block
    try:
        #0)Initial steps
        
        #Set global variables
        global urls
        global errors_df
        
        #Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #Treat provided output directory
        outputdir_path_def = treat_output_directory_parameter(outputdir_path)
        
        #1)Load input file, get urls and show information depending on execution mode
        if inputfile_type == 'ENA':
            urls = ENA_mode(inputfile_path, ena_download_column, plain_text_bool)
        else:
            urls = LINKS_mode(inputfile_path, plain_text_bool)
            
        #2)Check urls domain and modify them if needed
        check_urls_domain_and_treatment(inputfile_type)
        
        #3)Try to download files with parfive
        download_files_with_parfive(urls, parfive_verb, n_max_conn, outputdir_path_def, plain_text_bool)
        
        #4)Treat errors_df
        
        #Previous steps
        ##Get outfile name
        out_name_errors = 'errors_report.tsv'
        ##Treat output_directory parameter / Get full output file path
        outputfile_errors = treat_output_directory_parameter_outfiles(out_name_errors, outputdir_path)
        
        #Save erros dataframe if there are errors
        if len(errors_df) != 0:
            #Show error report message
            print(rich_text_colored('\nError Report Generated!','program_warning', plain_text_bool))
            #Show the number of detected error
            print(rich_text_colored('Total Download Errors Detected:', 'general_text', plain_text_bool), len(errors_df))
            #Show saved file message
            print(rich_text_colored('\nSaving report in file:', 'general_text', plain_text_bool))
            print(outputfile_errors)
            #Save df_merge_ENA
            errors_df.to_csv(outputfile_errors, header = True, index = False, sep = '\t')
    
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
           
    finally:
        #Print empty line for aesthetic purposes
        print('')

if __name__ == '__main__':
    main()
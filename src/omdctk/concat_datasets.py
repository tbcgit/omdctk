#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: concat_datasets.py
@author: Samuel Piquer-Esteban
@date: 21 Feb 2024 

"""

#Set program name
__program__ = 'concat_datasets.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, VALID_REQUIREDNESS, OMD_CTK_Exception, program_header, treat_headers_check,
                    check_headers, check_duplicates_in_variables_dict_column, check_required_variables_in_metadata_table, 
                    check_metadata_table_vars_in_dict, treat_output_directory_parameter_outfiles, check_na_in_pandas_dataframe,
                    check_existence_directory_parameter, treat_values_check, check_values, rich_text_colored)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
import pandas as pd
import os

#Program Constants
VARIABLES_FILE_HEADERS_4_CONCAT = ['variable', 'requiredness']

#Program functions
def treat_output_file_prefix_name(prefix, file_name):
    """
    This function creates the name for the file treating 
    the --output_name_prefix parameter.

    Parameters
    ----------
    prefix : str
        The provided output_name_prefix parameter.
    file_name : str
        The name of the file.

    Returns
    -------
    final_name : str
        The final name of the file.

    """
    if prefix is None:
        final_name = file_name
    else:
        final_name = prefix + '_' + file_name
    return final_name


def get_path_files_in_directory(directory, pattern):
    """
    This function gets a list with the files of the provided 
    directory for the provided pattern.

    Parameters
    ----------
    directory : str
        Path to the directory.
    pattern : str
        Extension pattern to recognize the wanted files.
        
    Raises
    ------
    OMD_CTK_Exception
        If there are no files detected raises an exception.

    Returns
    -------
    wanted_files : list
        List of the wanted files paths.
    
    """
    #Get files in directory
    files = os.listdir(directory)
    
    #Filter to retain only fastq.gz files
    wanted_files = [x for x in files if x.endswith(pattern)]
    
    #Check if fastq files were detected
    if len(wanted_files) == 0:
        frase0 = 'Error! There are no files with pattern "'
        frase1 = '" in the provided Input Directory!'
        raise OMD_CTK_Exception(''.join([frase0, pattern, frase1]))
    
    #Get full paths
    wanted_files_full_path = [os.path.join(directory, x) for x in wanted_files]
    
    return wanted_files_full_path


def get_path_files_in_project_mode(directory, pattern, color_treatment):
    """
    This function gets a list with the files with the provided pattern
    at the inmidiate subdir childs of the provided directory.

    Parameters
    ----------
    directory : str
        Path to the directory.
    pattern : str
        Extension pattern to recognize files of interest.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text
        
    Raises
    ------
    OMD_CTK_Exception
        If there are no final metadata files detected or more than one matching file raises an excepcion.

    Returns
    -------
    mt_files_paths : list
        List of the final metadata files paths.
    
    """
    #Set results list
    mt_files_paths = []
    
    #Get list of child subdirs
    datasets = next(os.walk(directory))[1]
    
    #If there are no child subdirs raise exception
    if len(datasets) == 0:
        raise OMD_CTK_Exception(''.join(['Error! For Search Mode: project. The provided Input Directory has no child subdirectories (datasets)!']))
    else:
        #Show message indicating number of datasets detected
        print('\n-',rich_text_colored('Total Datasets (child directories) detected:', 'general_text', color_treatment), len(datasets))
    #For each dataset check if there is a final metadata file
    for p in datasets:
        #Get child dir path
        temp_child_path = os.path.join(directory, p)
        
        #Get files in child directory
        temp_files = os.listdir(temp_child_path)
    
        #Filter to retain only files with extension
        mt_temp_files = [x for x in temp_files if x.endswith(pattern)]
        
        #Check if files were detected
        ##If there were none raise exception
        if len(mt_temp_files) == 0:
            raise OMD_CTK_Exception(''.join(['Error! For Dataset: ', p, '. There are no files with pattern extension "', pattern,'"!']))
        ##If there were more that 1 ask the user what to do
        elif len(mt_temp_files) > 1:
            
            #Initial message informing that there is more than one wanted matching file in dataset folder
            
            ##Show warning message
            print(rich_text_colored(''.join(['\nWarning! There is more than one file with pattern extension "', pattern,'" for the following Dataset: ', p,'!']), 'program_warning2', color_treatment))
            ##Show files in dataset folder
            print(rich_text_colored('\nThe following files were found:', 'general_text', color_treatment))
            for file in mt_temp_files:
                print(file)
        
            #Little menu to see what the user wants to do
            ##Get initial input
            print(rich_text_colored('\nDo you want to add all the found files anyway?', 'general_text', color_treatment))
            response = input('Provide a valid answer {y, n}:')
            ##While loop only valid answer
            while response not in ('y','n'):
                response = input('Provide a valid answer {y, n}:')
            ##Continue
            else:
                ##If it is okay continue and add all files
                if response == 'y':
                    for file in mt_temp_files:
                        mt_files_paths.append(os.path.join(temp_child_path, file))
                ##Else raise exception
                else:
                    raise OMD_CTK_Exception(''.join(['Error! For Dataset: ', p,'. There were unespected files with pattern extension "', pattern,'"!']))
        ##If there was only one file add normally
        else:
            mt_files_paths.append(os.path.join(temp_child_path, mt_temp_files[0]))
    
    return mt_files_paths


#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter input_directory
    parser.add_argument(
            '-i','--input_directory',
            action = 'store',
            required = True,
            help ='Input Directory. Indicate the path to the Input Directory with the Datasets Metadata files.'
    )
    ##Parameter variables_dictionary
    parser.add_argument(
            '-d','--variables_dictionary',
            action = 'store',
            required = True,
            help ='Variables Dictionary [Expected sep=TABS]. Indicate path to the Variables Dictionary file. See Documentation for more information and format details.'
    )
    ##Parameter search_mode
    parser.add_argument(
            '-s','--search_mode', 
            action = 'store',
            choices = ['simple','project'],
            required = False,
            default = 'simple',
            help = 'Search Mode (Optional) [Default:simple]. Indicate the selected mode to search metadata files for each dataset. Options: 1) Simple (all files are in the provided Input Directory) or 2) Project (the provided Input Directory has one folder for each dataset with their own files).'
    )
    ##Parameter metadata_pattern
    parser.add_argument(
            '-p','--metadata_pattern', 
            action = 'store',
            required = False,
            default = '_metadata_final.tsv',
            help = 'Dataset Metadata File Pattern (Optional) [Default:"_metadata_final.tsv"]. Indicate the pattern to identify Dataset Metadata files.'
    )
    ##Parameter output_file_prefix
    parser.add_argument(
            '-op','--output_name_prefix', 
            action = 'store',
            required = False,
            help = 'Output Name Prefix (Optional). Indicate prefix name for the output files.'
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
    
    #Process common arguments
    args = parser.parse_args()
    input_path = args.input_directory
    search_mode = args.search_mode
    files_pattern = args.metadata_pattern
    outputdir_path = args.output_directory
    output_name_prefix = args.output_name_prefix
    variables_dictionary_path = args.variables_dictionary
    plain_text_bool = args.plain_text
    
    #Show Program headers
    print('')
    program_header('#',  __program__, 64, 2, plain_text_bool)
    print('')
    
    #Show Program parameters
    print(rich_text_colored('Program Parameters:', 'section_header', plain_text_bool))
    print(tabulate(vars(args).items(), headers = ['Argument', 'Value'], tablefmt = 'simple_outline'))
    
    #Check if the files exists and load files/ Try-except-else
    try:
        #0) Initial checks and steps
        
        ##Check that provided input directory exist
        check_existence_directory_parameter(input_path, 'Input', '--input_directory')
    
        ##Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #1) Search files in Input Directory
        #Section header message
        print(rich_text_colored('\nSearching Metadata Files:', 'section_header', plain_text_bool))
        
        ##Get list of metadata files for each dataset depending on selected search mode
        if search_mode == 'project':
            input_files_full = get_path_files_in_project_mode(input_path, files_pattern, plain_text_bool)
        else:
            input_files_full = get_path_files_in_directory(input_path, files_pattern)
        
        ##Show message indicating the number of metadata tables detected
        print('\n-',rich_text_colored('Total Metadata Files detected:', 'general_text', plain_text_bool), len(input_files_full))
        
        #2) Load files and do checks
        #Section header message
        print(rich_text_colored('\nLoading Files:', 'section_header', plain_text_bool))
        
        #2.1.Try to load Variables Dictionary File as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('\nVariables Dictionary file:', 'general_text', plain_text_bool))
        print(variables_dictionary_path)
        ##Try to load the filterfile
        variables_dict_table = pd.read_csv(variables_dictionary_path, sep = '\t', dtype = 'string')
        
        #Do initial checks on Variables Dictionary file for columns of interest
        ##Check headers
        check_headers_variablesfile = check_headers(VARIABLES_FILE_HEADERS_4_CONCAT, variables_dict_table)
        frase0_1 = 'Error! Some of the needed headers are not in the Variables Dictionary!\n Check your variables file!\n'
        frase1_1 = '\nThe headers needed are:'
        treat_headers_check(VARIABLES_FILE_HEADERS_4_CONCAT, check_headers_variablesfile, variables_dict_table, frase0_1, frase1_1, plain_text_bool)
        ##Check duplicates in 'variable' column
        check_duplicates_in_variables_dict_column(variables_dict_table, 'variable', 5)
        ##Check NAs in 'variable' column
        frase2 = 'Error! Some of the values of the "variable" column in the Variables Dictionary are NAs!\n Check your variables file!'
        check_na_in_pandas_dataframe(variables_dict_table['variable'], frase2)
        ##Check valid 'requiredness' column values
        check_requiredness_values = check_values(VALID_REQUIREDNESS, variables_dict_table['requiredness'])
        frase0_3 = 'Error! Some of the values of the "requiredness" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_3 = '\nValid values are:'
        treat_values_check(VALID_REQUIREDNESS, check_requiredness_values, frase0_3, frase1_3, plain_text_bool)
        
        #2.2.Try to load metadata files and do checks
        
        #Prepare variables
        ##Set list with metadata tables paths
        mt_tables = []
        ##Get list of required columns
        required_cols = variables_dict_table[variables_dict_table['requiredness']=='required']['variable']
        ##Get list of all possible columns
        universe_cols = list(variables_dict_table['variable'])
        
        ##Show loading files message
        print(rich_text_colored('\nDatasets Metadata files:', 'general_text', plain_text_bool))
        
        ##Try to read each metadata table and do checks
        for file in input_files_full:
            ##Show metadata file path
            print(file)
            ##Load temp file
            temp_mt = pd.read_csv(file, sep = '\t')
            ##Check required variables
            check_required_variables_in_metadata_table(temp_mt, required_cols, plain_text_bool)
            ##Check that all columns in metadata are in the variables dictionary
            check_metadata_table_vars_in_dict(temp_mt, universe_cols, plain_text_bool)
            ##Save in mt_tables results list
            mt_tables.append(temp_mt)
        
        #3) Concatenate metadata tables
        ##Section header message
        print(rich_text_colored('\nConcatenating Metadata Tables:', 'section_header', plain_text_bool))
        
        ##Concat metadata tables
        result = pd.concat(mt_tables)
        
        ##Show dimentions of final metadata table
        print('\n-', rich_text_colored('Number of final Rows:', 'general_text', plain_text_bool),result.shape[0])
        print('-', rich_text_colored('Number of final Columns:', 'general_text', plain_text_bool),result.shape[1])
        
        ##Save final file
        ###Get final file name
        out_table = treat_output_file_prefix_name(output_name_prefix, 'concatenated_final_metadata.tsv')
        ###Treat output_directory parameter / Get full output file path
        outputfile_table = treat_output_directory_parameter_outfiles(out_table, outputdir_path)
        ###Message
        print(rich_text_colored('\nSaving results in file:', 'general_text', plain_text_bool))
        print(outputfile_table)
        ###Save result dataframe
        result.to_csv(outputfile_table, header = True, index = False, sep = '\t')
    
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
    
    finally:
        #Print empty line for aesthetic purposes
        print('')    

if __name__ == '__main__':
    main()
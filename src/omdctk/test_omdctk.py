#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: test_omdctk.py
@author: Samuel Piquer-Esteban
@date: 12 Mar 2024

"""

#Set program name
__program__ = 'test_omdctk.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, program_header, OMD_CTK_Exception,
                    check_existence_directory_parameter, rich_text_colored)

#Import third-party modules
from pathlib import Path
from argparse import ArgumentParser
from tabulate import tabulate
import subprocess
import shutil
import os

#Program Constants

#0)Set test_info_path
TEST_INFO_PATH = os.path.join(Path(__file__).parent.absolute(), 'test_info')

##Set Datasets Directory names
ENA_DATASET_DIR_NAME = 'ENA_Dataset_Example'
EXTERNAL_DATASET_DIR_NAME = 'External_Dataset_Example'

##Set Child Fastqs Subdirectory names
DOWNLOAD_SUBDIR_NAME = 'downloads'
TREATMENT_SUBDIR_NAME = 'treated_files'

##Check fastqs program log stdout pattern
CHECK_FASTQS_PATTERN = 'Main Information:'

#A) ENA Workflow

##A1)download_metadata_ENA

STUDY_ACCESSION = 'PRJEB10949'

RUN1_OUTPUTFILES = ['PRJEB10949_ENA_browser.tsv', 'PRJEB10949_mg-toolkit.tsv', 'PRJEB10949_ENA_metadata.tsv']

##A2)merge_metadata

ENA_MERGE_COLUMN = 'run_accession'

PUBLICATION_MERGE_COLUMN = 'run_accessions'

PUBLICATION_MTD_FILE_NAME = 'PRJEB10949_publication_example.tsv'

ENA_MTD_FILE_NAME = 'PRJEB10949_ENA_metadata.tsv'

RUN2_OUTPUTFILES = ['merged_PRJEB10949_ENA_metadata.tsv']

##A3)check_metadata_ENA

MERGED_MTD_FILE_NAME = 'merged_PRJEB10949_ENA_metadata.tsv'

CHECK_MTD_REF_LOG = 'PRJEB10949_check_metadata_ENA_reference.log'

CHECK_MTD_LOG = 'PRJEB10949_check_metadata_ENA_stdout.log'

CHECK_MTD_PATTERN = "Runs' Stats:\n"

##A4)filter_metadata

FILTER_FILE_NAME = 'PRJEB10949_filterfile_example.tsv'

RUN4_OUTPUTFILES = ['filtered_merged_PRJEB10949_ENA_metadata.tsv']

##A5)download_fastqs

FILTERED_MTD_FILE_NAME = 'filtered_merged_PRJEB10949_ENA_metadata.tsv'

RUN5_OUTPUTFILES = ['ERR1049851_1.fastq.gz', 'ERR1049851_2.fastq.gz', 'ERR1049852_1.fastq.gz', 'ERR1049852_2.fastq.gz', 'ERR1049853_1.fastq.gz',
                    'ERR1049853_2.fastq.gz', 'ERR1049854_1.fastq.gz', 'ERR1049854_2.fastq.gz', 'ERR1049855_1.fastq.gz', 'ERR1049855_2.fastq.gz',
                    'ERR1049856_1.fastq.gz', 'ERR1049856_2.fastq.gz', 'ERR1049857_1.fastq.gz', 'ERR1049857_2.fastq.gz', 'ERR1049859_1.fastq.gz',
                    'ERR1049859_2.fastq.gz', 'ERR1049861_1.fastq.gz', 'ERR1049861_2.fastq.gz', 'ERR1049911_1.fastq.gz', 'ERR1049911_2.fastq.gz',
                    'ERR1049913_1.fastq.gz', 'ERR1049913_2.fastq.gz', 'ERR1049915_1.fastq.gz', 'ERR1049915_2.fastq.gz', 'ERR1049917_1.fastq.gz',
                    'ERR1049917_2.fastq.gz', 'ERR1049919_1.fastq.gz', 'ERR1049919_2.fastq.gz', 'ERR1049921_1.fastq.gz', 'ERR1049921_2.fastq.gz',
                    'ERR1049923_1.fastq.gz', 'ERR1049923_2.fastq.gz', 'ERR1049924_1.fastq.gz', 'ERR1049924_2.fastq.gz', 'ERR1049925_1.fastq.gz',
                    'ERR1049925_2.fastq.gz', 'ERR1049926_1.fastq.gz', 'ERR1049926_2.fastq.gz', 'ERR1049927_1.fastq.gz', 'ERR1049927_2.fastq.gz',
                    'ERR1049928_1.fastq.gz', 'ERR1049928_2.fastq.gz']

##A6)check_fastqs

CHECK_FASTQS_REF_LOG = 'PRJEB10949_check_fastqs_ENA_reference.log'

CHECK_FASTQS_LOG = 'PRJEB10949_check_fastqs_ENA_stdout.log'

##A7)make_treatment_template

RUN7_OUTPUTFILES = ['raw_treatment_template_filtered_merged_PRJEB10949_ENA_metadata.tsv']

##A8)treat_fastqs

TREATMENT_FILE_NAME = 'treatment_template_filtered_PRJEB10949_merged_metadata_example.tsv'

RUN8_OUTPUTFILES = ['ERR1049851_1.fastq.gz', 'ERR1049851_2.fastq.gz', 'ERR1049852_1.fastq.gz', 'ERR1049852_2.fastq.gz', 'ERR1049853_1.fastq.gz',
                    'ERR1049853_2.fastq.gz', 'Ileum-1_1.fastq.gz', 'Ileum-1_2.fastq.gz', 'Ileum-2_1.fastq.gz', 'Ileum-2_2.fastq.gz', 
                    'Ileum-3_1.fastq.gz', 'Ileum-3_2.fastq.gz', 'Ileum-4_1.fastq.gz', 'Ileum-4_2.fastq.gz', 'Ileum-5_1.fastq.gz', 
                    'Ileum-5_2.fastq.gz', 'Ileum-6_1.fastq.gz', 'Ileum-6_2.fastq.gz', 'MAT-1_1.fastq.gz', 'MAT-1_2.fastq.gz', 'MAT-2_1.fastq.gz',
                    'MAT-2_2.fastq.gz', 'MAT-3_1.fastq.gz', 'MAT-3_2.fastq.gz']

##A9)treat_metadata

RUN9_OUTPUTFILES = ['treated_filtered_merged_PRJEB10949_ENA_metadata.tsv']

#B) External Dataset Workflow

##B1)merge_metadata

GMAIN_MERGE_COLUMN = 'Sample_name'

GPUBLICATION_MERGE_COLUMN = 'sample_id'

GPUBLICATION_MTD_FILE_NAME = 'CRA001372_publication_metadata_example.tsv'

GMAIN_MTD_FILE_NAME = 'CRA001372_main_metadata_example.tsv'

RUNB1_OUTPUTFILES = ['merged_CRA001372_main_metadata_example.tsv']

##B2)filter_metadata

EXT_MERGED_MTD_FILE_NAME = 'merged_CRA001372_main_metadata_example.tsv'

EXT_FILTER_FILE_NAME = 'CRA001372_filterfile_example.tsv'

RUNB2_OUTPUTFILES = ['filtered_merged_CRA001372_main_metadata_example.tsv']

##B3)download_fastqs LINKS mode

EXT_URLS_FILE = 'filtered_CRA001372_URLS_example.txt'

RUNB3_OUTPUTFILES = ['CRR044580.fq.gz','CRR044582.fq.gz','CRR044584.fq.gz','CRR044586.fq.gz','CRR044588.fq.gz',
                     'CRR044590.fq.gz','CRR044592.fq.gz','CRR044594.fq.gz','CRR044596.fq.gz','CRR044598.fq.gz',
                     'CRR044600.fq.gz','CRR044602.fq.gz','CRR044581.fq.gz','CRR044583.fq.gz','CRR044585.fq.gz',
                     'CRR044587.fq.gz','CRR044589.fq.gz','CRR044591.fq.gz','CRR044593.fq.gz','CRR044595.fq.gz',
                     'CRR044597.fq.gz','CRR044599.fq.gz','CRR044601.fq.gz','CRR044603.fq.gz']

##B4)check_fastqs Generic mode

EXT_CHECK_FASTQS_REF_LOG = 'CRA001372_check_fastqs_Generic_reference.log'

EXT_CHECK_FASTQS_LOG = 'CRA001372_check_fastqs_Generic_stdout.log'

EXT_MANIFEST_FILE = 'filtered_manifest_CRA001372_example.tsv'

EXT_FILTERED_MTD_FILE = 'filtered_merged_CRA001372_main_metadata_example.tsv'

##B5)make_treatment_template Generic mode

RUNB5_OUTPUTFILES = ['raw_treatment_template_filtered_manifest_CRA001372_example.tsv']

##B6)treat_fastqs

GTREATMENT_FILE_NAME = 'treatment_template_filtered_CRA001372_example.tsv'

RUNB6_OUTPUTFILES = ['CRR044580.fq.gz','CRR044581.fq.gz','CRR044582.fq.gz','CRR044583.fq.gz','CRR044584.fq.gz',
                     'CRR044585.fq.gz','CRR044586.fq.gz','CRR044587.fq.gz','CRR044588.fq.gz','CRR044589.fq.gz',
                     'CRR044590.fq.gz','CRR044591.fq.gz','Soil2Ha.fq.gz','Soil2Hb.fq.gz','Soil2Hc.fq.gz',
                     'Soil2Hd.fq.gz','Soil2He.fq.gz','Soil2Hf.fq.gz','Soil2La.fq.gz','Soil2Lb.fq.gz',
                     'Soil2Lc.fq.gz','Soil2Ld.fq.gz','Soil2Le.fq.gz','Soil2Lf.fq.gz']

##B7)treat_metadata Generic mode

RUNB7_OUTPUTFILES = ['treated_filtered_merged_CRA001372_main_metadata_example.tsv']

#C) Meta-analysis Workflow

##C1)concat_datasets

CURATED_ENA_MTD_FILE_NAME = 'curated_PRJEB10949_ENA_example_metadata_final.tsv'

CURATED_EXT_MTD_FILE_NAME = 'curated_CRA001372_external_example_metadata_final.tsv'

VAR_DICT_FILE_NAME = 'variables_dictionary_example.tsv'

RUNC1_OUTPUTFILES = ['example_concatenated_final_metadata.tsv']

##C2)check_metadata_values

CHECK_CONCAT_MTD_REF_LOG = 'check_metadata_values_example_reference.log'

CHECK_CONCAT_MTD_LOG = 'check_metadata_values_example_stdout.log'

CONCAT_MTD = 'example_concatenated_final_metadata.tsv'

CHECK_MT_VALS_PATTERN = 'Checking Provided Variables:'

#Program Functions
def check_generic_run(run_object, output_dir, expected_files, n_elements, color_treatment):
    """
    This function checks a generic test for a run program.

    Parameters
    ----------
    run_object : CompletedProcess
        Result of running the program with subprocess.run().
    output_dir : str
        Path to the Output Directory.
    expected_files : list
        List with the expected output files produced by the program.
    n_elements : int
        The number of elements to be printed by line.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If an unsuccessfull exit code is returned or the expected
        output files are not found raises an exception.

    Returns
    -------
    None.

    """
    if run_object.returncode != 0:
        raise OMD_CTK_Exception('Error! Running the program returned an exit code different than 0!')
    else:
        #List output_dir
        files = os.listdir(output_dir)
        #Check if all expected files
        if all([i in files for i in expected_files]):
            print(rich_text_colored('\nResult:', 'general_text', color_treatment), rich_text_colored('\nSuccess! All expected files have been generated!', 'acceptable', color_treatment))
        else:
            print(rich_text_colored('\nResult:', 'general_text', color_treatment), rich_text_colored('\nNot all expected files have been generated!', 'program_warning', color_treatment))
            #Prepare exception frase
            frase0 = 'Error! Not all expected files were found after running the last program!\n'
            frase1 = '\nExpected output files are:\n'
            frase = [frase0, frase1]
            #Split file names
            splited_list = [expected_files[i:i+n_elements] for i in range(0, len(expected_files), n_elements)]
            for i in splited_list:
                #Format depending on if it is the last line
                if i == splited_list[-1]:
                    temp_line = ', '.join(i)
                else:
                    temp_line = ', '.join(i+['\n'])
                #Append line
                frase.append(temp_line)
            #Raise exception
            raise OMD_CTK_Exception(' '.join(frase))


def copy_file(file_name, input_file_path, output_directory_path):
    """
    This function copies a file in the provided Output Directory.

    Parameters
    ----------
    file_name : str
        The file name to be copied.
    input_file_path : str
        The provided input file path.
    output_directory_path : str
        The provided Output Directory.
        
    Raises
    ------
    OMD_CTK_Exception
        If a file could not be successfully copied raises an exception.

    Returns
    -------
    None.
        
    """
    #Set input and output paths
    output_file_path = os.path.join(output_directory_path, file_name)
    
    #Try to copy file in output directory
    try:
        shutil.copyfile(input_file_path, output_file_path)    
    except:
        raise OMD_CTK_Exception('Error! This file could not be successfully copied in the Output Directory!')


def read_reference_log_file(input_file_path):
    """
    This function reads the provided reference log file.

    Parameters
    ----------
    input_file_path : str
        The provided input file path.

    Raises
    ------
    OMD_CTK_Exception
        If the file could not be read raises an exception.

    Returns
    -------
    str
        Content of the reference log file.

    """
    try:
        with open(input_file_path) as my_file:
            return my_file.read()
    except:
        raise OMD_CTK_Exception('Error! The reference log file could not be successfully read!')


def save_log_file(content, output_file_path):
    """
    This function saves the content of the provided log file.

    Parameters
    ----------
    content : str
        Log file content.
    output_file_path : str
        The provided output file path.

    Raises
    ------
    OMD_CTK_Exception
        If the file could not be saved raises an exception.

    Returns
    -------
    None.

    """
    try:
        with open(output_file_path, 'w') as output:
            output.write(content)
    except:
        raise OMD_CTK_Exception('Error! The log file could not be successfully saved!')


def check_generic_log_run(pattern, run_object, stdout_log_name, reference_log, reference_log_name, output_directory_path, color_treatment):
    """
    This function checks if the stdout from a check program matches
    the expected output from the package reference log.

    Parameters
    ----------
    pattern : str
        String pattern to process the log files.
    run_object : CompletedProcess
        Result of running the check program with subprocess.run().
    stdout_log_name : str
        Name for the stdout log file.
    reference_log : str
        Content of the package's reference log from the program.
    reference_log_name : str
        Name for the reference log file.
    output_dir : str
        Path to the Output Directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If incongruencies are detected raises an exception.

    Returns
    -------
    None.

    """
    if run_object.returncode != 0:
        raise OMD_CTK_Exception('Error! Running the program returned an exit code different than 0!')
    else:
        #Process logs
        ##Process run object stdout
        stdout_processed = run_object.stdout[run_object.stdout.index(pattern) : len(run_object.stdout)]
        ##Process reference log
        reference_log_processed = reference_log[reference_log.index(pattern) : len(reference_log)]
        
        #Check if outputs are equal were it matters
        if reference_log_processed == stdout_processed:
            print(rich_text_colored('\nResult:', 'general_text', color_treatment), rich_text_colored('\nSuccess! The program stdout matches the expected output!', 'acceptable', color_treatment))
        else:
            print(rich_text_colored('\nResult:', 'general_text', color_treatment), rich_text_colored('\nThe program stdout does not match the expected output!', 'program_warning', color_treatment))
            
            #Save log files
            ## Ref log
            print('Saving', reference_log_name)
            ref_log_path = os.path.join(output_directory_path, reference_log_name)
            save_log_file(reference_log, ref_log_path)
            ## stdout log
            print('Saving', stdout_log_name)
            stdout_log_path = os.path.join(output_directory_path, stdout_log_name)
            save_log_file(run_object.stdout, stdout_log_path)
            
            #Prepare exception frase
            frase0 = 'Error! The stdout does not match the expected output after running the last program!'
            frase1 = '\n Log files have been saved in Output Directory!'
            frase = [frase0, frase1]
            #Raise exception
            raise OMD_CTK_Exception(' '.join(frase))


def create_subdirectory(subdir_name, output_directory_path):
    """
    This function creates a subdirectory with the provided name
    in the Output Directory.

    Parameters
    ----------
    subdir_name : str
        Name for the subdirectory to be created.
    output_dir : str
        Path to the Output Directory.

    Raises
    ------
    OMD_CTK_Exception
        If the subdirectory coulb not be created raises an exception.

    Returns
    -------
    None.

    """
    try:
        path = os.path.join(output_directory_path, subdir_name)
        os.mkdir(path)
    except:
        raise OMD_CTK_Exception('Error! The subdirectory could not be successfully created!')
        

#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter output_directory
    parser.add_argument(
            '-o','--output_directory', 
            action = 'store',
            required = True,
            help = 'Output Directory. Indicate the path to the Output Directory to save the resulting files.'
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
    outputdir_path = args.output_directory
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
        #0)Initial checks and message
        
        #Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #Message
        print(rich_text_colored('\nThis may take a while...', 'program_warning2', plain_text_bool))
        
        #A) ENA Workflow
        #Section header message
        print(rich_text_colored('\nA) Testing ENA Dataset Workflow:', 'section_header', plain_text_bool))
        #Get full path and create download subdirectory if it does not exit
        ena_dataset_subdir_path = os.path.join(outputdir_path, ENA_DATASET_DIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(ena_dataset_subdir_path) == False:
            print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
            print('Creating ENA Dataset Subdirectory')
            create_subdirectory(ENA_DATASET_DIR_NAME, outputdir_path)
        
        #A1)download_metadata
        
        #Messages
        print('\nA1)', rich_text_colored('Testing download_metadata_ENA.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('download_metadata_ENA -p', STUDY_ACCESSION , '-o', ena_dataset_subdir_path, '--plain_text')
        
        #Run program
        run1 = subprocess.run(['download_metadata_ENA','-p', STUDY_ACCESSION, '-o', ena_dataset_subdir_path, '--plain_text'], stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run1, ena_dataset_subdir_path, RUN1_OUTPUTFILES, 5, plain_text_bool)
        
        #A2)merge_metadata
        
        #Messages
        print('\nA2)', rich_text_colored('Testing merge_metadata.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Extra Metadata Table Example in ENA_Dataset_Example Directory')
        
        #Get Publication Table file path in package and Copy file in output_dir
        publication_table_path = os.path.join(TEST_INFO_PATH, PUBLICATION_MTD_FILE_NAME)
        copy_file(PUBLICATION_MTD_FILE_NAME, publication_table_path, ena_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('merge_metadata -m', ENA_MTD_FILE_NAME, '-mc', ENA_MERGE_COLUMN, '-e', PUBLICATION_MTD_FILE_NAME, '-ec', PUBLICATION_MERGE_COLUMN, '-o', ena_dataset_subdir_path, '--plain_text')
        
        #Run program
        publication_table_path_def = os.path.join(ena_dataset_subdir_path, PUBLICATION_MTD_FILE_NAME)
        ena_table_path = os.path.join(ena_dataset_subdir_path, ENA_MTD_FILE_NAME)
        run2_parameters = ['merge_metadata','-m', ena_table_path, '-mc', ENA_MERGE_COLUMN, '-e', publication_table_path_def, '-ec', PUBLICATION_MERGE_COLUMN, '-o', ena_dataset_subdir_path, '--plain_text']
        run2 = subprocess.run(run2_parameters, stdout = subprocess.DEVNULL)
        
        #Check outputs and run
        check_generic_run(run2, ena_dataset_subdir_path, RUN2_OUTPUTFILES, 5, plain_text_bool)
        
        #A3)check_metadata_ENA
        
        #Messages
        print('\nA3)', rich_text_colored('Testing check_metadata_ENA.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Reading package reference output log file')
        
        check_metadata_log_path = os.path.join(TEST_INFO_PATH, CHECK_MTD_REF_LOG)
        #Get check_metadata_ENA reference log file path in package and read
        check_metadata_reference_log = read_reference_log_file(check_metadata_log_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('check_metadata_ENA -t', MERGED_MTD_FILE_NAME, '--plain_text')
        
        #Run program
        merged_table_path = os.path.join(ena_dataset_subdir_path, MERGED_MTD_FILE_NAME)
        run3 = subprocess.run(['check_metadata_ENA','-t', merged_table_path, '--plain_text'], stdout = subprocess.PIPE, text = True) 

        #Check outputs and run
        check_generic_log_run(CHECK_MTD_PATTERN, run3, CHECK_MTD_LOG, check_metadata_reference_log, CHECK_MTD_REF_LOG, ena_dataset_subdir_path, plain_text_bool)
        
        #A4)filter_metadata
        
        #Messages
        print('\nA4)', rich_text_colored('Testing filter_metadata.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Filter Table Example in ENA_Dataset_Example Directory')
        
        #Get Filter Table file path in package and Copy file in output_dir
        filter_table_path = os.path.join(TEST_INFO_PATH, FILTER_FILE_NAME)
        copy_file(FILTER_FILE_NAME, filter_table_path, ena_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('filter_metadata -t', MERGED_MTD_FILE_NAME, '-f', FILTER_FILE_NAME,'-o', ena_dataset_subdir_path, '--plain_text')
        
        #Run program
        filter_table_path_def = os.path.join(ena_dataset_subdir_path, FILTER_FILE_NAME)
        run4_parameters = ['filter_metadata','-t', merged_table_path, '-f', filter_table_path_def, '-o', ena_dataset_subdir_path, '--plain_text']
        run4 = subprocess.run(run4_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run4, ena_dataset_subdir_path, RUN4_OUTPUTFILES, 5, plain_text_bool)
        
        #A5)download_fastqs
        
        #Message
        print('\nA5)', rich_text_colored('Testing download_fastqs.py in ENA mode:', 'test_color', plain_text_bool))
        
        #Get full path and create download subdirectory if it does not exit
        download_subdir_path = os.path.join(ena_dataset_subdir_path, DOWNLOAD_SUBDIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(download_subdir_path) == False:
            print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
            print('Creating Download Subdirectory in ENA_Dataset_Example Directory')
            create_subdirectory(DOWNLOAD_SUBDIR_NAME, ena_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('download_fastqs -i', FILTERED_MTD_FILE_NAME,'-o', download_subdir_path, '--plain_text')
        
        #Run program
        filtered_table_path = os.path.join(ena_dataset_subdir_path, FILTERED_MTD_FILE_NAME)
        run5_parameters = ['download_fastqs','-i', filtered_table_path, '-o', download_subdir_path, '--plain_text']
        run5 = subprocess.run(run5_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run5, download_subdir_path, RUN5_OUTPUTFILES, 5, plain_text_bool)
        
        #A6)check_fastqs
        
        #Messages
        print('\nA6)', rich_text_colored('Testing check_fastqs.py in ENA mode:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Reading package reference output log file')
        
        #Get check_metadata_ENA reference log file path in package and read
        check_fastqs_log_path = os.path.join(TEST_INFO_PATH, CHECK_FASTQS_REF_LOG)
        check_fastqs_reference_log = read_reference_log_file(check_fastqs_log_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('check_fastqs -t', FILTERED_MTD_FILE_NAME, '-d', download_subdir_path, '--md5_check', '--plain_text')
        
        #Run program
        run6 = subprocess.run(['check_fastqs','-t', filtered_table_path, '-d', download_subdir_path, '--md5_check', '--plain_text'], stdout = subprocess.PIPE, text = True) 

        #Check outputs and run
        check_generic_log_run(CHECK_FASTQS_PATTERN, run6, CHECK_FASTQS_LOG, check_fastqs_reference_log, CHECK_FASTQS_REF_LOG, ena_dataset_subdir_path, plain_text_bool)
        
        #A7)make_treatment_template in ENA mode
        
        #Messages
        print('\nA7)', rich_text_colored('Testing make_treatment_template.py in ENA mode:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('make_treatment_template -i', FILTERED_MTD_FILE_NAME, '-d', download_subdir_path, '--extra_sample_columns sample_column', '-o', ena_dataset_subdir_path, '--plain_text')
        
        #Run program
        run7_parameters = ['make_treatment_template','-i', filtered_table_path, '-d', download_subdir_path, '--extra_sample_columns','sample_column', '-o', ena_dataset_subdir_path, '--plain_text']
        run7 = subprocess.run(run7_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run7, ena_dataset_subdir_path, RUN7_OUTPUTFILES, 5, plain_text_bool)
        
        #A8)treat_fastqs
        
        #Message
        print('\nA8)', rich_text_colored('Testing treat_fastqs.py:', 'test_color', plain_text_bool))
        
        #Get Filter Table file path in package and Copy file in output_dir
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Treatment Template Example in ENA_Dataset_Example Directory')
        treatment_table_path = os.path.join(TEST_INFO_PATH, TREATMENT_FILE_NAME)
        copy_file(TREATMENT_FILE_NAME, treatment_table_path, ena_dataset_subdir_path)
        
        #Get full path and create treated_files subdirectory if it does not exit
        treated_files_subdir_path = os.path.join(ena_dataset_subdir_path, TREATMENT_SUBDIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(treated_files_subdir_path) == False:
            print('Creating Treated Files Subdirectory in ENA_Dataset_Example Directory')
            create_subdirectory(TREATMENT_SUBDIR_NAME, ena_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('treat_fastqs -t', TREATMENT_FILE_NAME, '-i', download_subdir_path, '-o', treated_files_subdir_path, '--plain_text')
        
        #Run program
        treatment_table_path_def = os.path.join(ena_dataset_subdir_path, TREATMENT_FILE_NAME)
        run8_parameters = ['treat_fastqs','-t', treatment_table_path_def, '-i', download_subdir_path, '-o', treated_files_subdir_path, '--plain_text']
        run8 = subprocess.run(run8_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run8, treated_files_subdir_path, RUN8_OUTPUTFILES, 5, plain_text_bool)
        
        ##A9)treat_metadata in ENA mode
        
        #Message
        print('\nA9)', rich_text_colored('Testing treat_metadata.py in ENA mode:', 'test_color', plain_text_bool))
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('treat_metadata -t', TREATMENT_FILE_NAME, '-m', FILTERED_MTD_FILE_NAME, '-o', ena_dataset_subdir_path, '--extra_no_warning_columns Run Sample run_accessions run_label', '--plain_text')
        
        #Run program
        run9_parameters = ['treat_metadata','-t', treatment_table_path_def, '-m', filtered_table_path, '-o', ena_dataset_subdir_path, '--extra_no_warning_columns', 'Run', 'Sample', 'run_accessions', 'run_label', '--plain_text']
        run9 = subprocess.run(run9_parameters, stdout = subprocess.DEVNULL)
        
        #Check outputs and run
        check_generic_run(run9, ena_dataset_subdir_path, RUN9_OUTPUTFILES, 5, plain_text_bool)
        
        #B) External Dataset Workflow
        #Section header message
        print(rich_text_colored('\nB) Testing External Dataset Workflow:', 'section_header', plain_text_bool))
        #Get full path and create download subdirectory if it does not exit
        external_dataset_subdir_path = os.path.join(outputdir_path, EXTERNAL_DATASET_DIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(external_dataset_subdir_path) == False:
            print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
            print('Creating External Dataset Subdirectory')
            create_subdirectory(EXTERNAL_DATASET_DIR_NAME, outputdir_path)
        
        ##B1)merge_metadata
        #Messages
        print('\nB1)', rich_text_colored('Testing merge_metadata.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))

        #Get Generic Main Metadata Table file path in package and Copy file in output_dir
        print('Coping Generic Main Metadata Table Example in External_Dataset_Example Directory')
        gmain_table_path = os.path.join(TEST_INFO_PATH, GMAIN_MTD_FILE_NAME)
        copy_file(GMAIN_MTD_FILE_NAME, gmain_table_path, external_dataset_subdir_path)
        
        #Get Publication Table file path in package and Copy file in output_dir
        print('Coping Generic Extra Metadata Table Example in External_Dataset_Example Directory')
        gpublication_table_path = os.path.join(TEST_INFO_PATH, GPUBLICATION_MTD_FILE_NAME)
        copy_file(GPUBLICATION_MTD_FILE_NAME, gpublication_table_path, external_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('merge_metadata -m', GMAIN_MTD_FILE_NAME, '-mc', GMAIN_MERGE_COLUMN, '-e', GPUBLICATION_MTD_FILE_NAME, '-ec', GPUBLICATION_MERGE_COLUMN, '-o', external_dataset_subdir_path, '-es _publication --plain_text')
        
        #Run program
        gpublication_table_path_def = os.path.join(external_dataset_subdir_path, GPUBLICATION_MTD_FILE_NAME)
        gmain_table_path_def = os.path.join(external_dataset_subdir_path, GMAIN_MTD_FILE_NAME)
        runb1_parameters = ['merge_metadata','-m', gmain_table_path_def, '-mc', GMAIN_MERGE_COLUMN, '-e', gpublication_table_path_def, '-ec', GPUBLICATION_MERGE_COLUMN, '-o', external_dataset_subdir_path, '-es','_publication', '--plain_text']
        runb1 = subprocess.run(runb1_parameters, stdout = subprocess.DEVNULL)
        
        #Check outputs and run
        check_generic_run(runb1, external_dataset_subdir_path, RUNB1_OUTPUTFILES, 5, plain_text_bool)

        ##B2)filter_metadata
        
        #Messages
        print('\nB2)', rich_text_colored('Testing filter_metadata.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Filter Table Example in External_Dataset_Example Directory')
        #Get Filter Table file path in package and Copy file in output_dir
        ext_filter_table_path = os.path.join(TEST_INFO_PATH, EXT_FILTER_FILE_NAME)
        copy_file(EXT_FILTER_FILE_NAME, ext_filter_table_path, external_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('filter_metadata -t', EXT_MERGED_MTD_FILE_NAME, '-f', EXT_FILTER_FILE_NAME,'-o', external_dataset_subdir_path, '--plain_text')
        
        #Run program
        ext_filter_table_path_def = os.path.join(external_dataset_subdir_path, EXT_FILTER_FILE_NAME)
        ext_merged_table_path = os.path.join(external_dataset_subdir_path, EXT_MERGED_MTD_FILE_NAME)
        runb2_parameters = ['filter_metadata','-t', ext_merged_table_path, '-f', ext_filter_table_path_def, '-o', external_dataset_subdir_path, '--plain_text']
        runb2 = subprocess.run(runb2_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(runb2, external_dataset_subdir_path, RUNB2_OUTPUTFILES, 5, plain_text_bool)

        ##B3)download_fastqs LINKS mode
        
        #Message
        print('\nB3)', rich_text_colored('Testing download_fastqs.py in LINKS mode:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        
        #Get full path and create download subdirectory if it does not exit
        ext_download_subdir_path = os.path.join(external_dataset_subdir_path, DOWNLOAD_SUBDIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(ext_download_subdir_path) == False:
            print('Creating Download Subdirectory in External_Dataset_Example Directory')
            create_subdirectory(DOWNLOAD_SUBDIR_NAME, external_dataset_subdir_path)
        
        #Get URLS file path in package and Copy file in output_dir
        print('Coping URLs TXT File Example in External_Dataset_Example Directory')
        ext_urls_file_path = os.path.join(TEST_INFO_PATH, EXT_URLS_FILE)
        copy_file(EXT_URLS_FILE, ext_urls_file_path, external_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('download_fastqs -m LINKS -i', EXT_URLS_FILE,'-o', ext_download_subdir_path, '--plain_text')
        
        #Run program
        ext_urls_file_path_def = os.path.join(external_dataset_subdir_path, EXT_URLS_FILE)
        runb3_parameters = ['download_fastqs','-m','LINKS','-i', ext_urls_file_path_def, '-o', ext_download_subdir_path, '--plain_text']
        runb3 = subprocess.run(runb3_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(runb3, ext_download_subdir_path, RUNB3_OUTPUTFILES, 5, plain_text_bool)

        ##B4)check_fastqs Generic mode
        
        #Messages
        print('\nB4)', rich_text_colored('Testing check_fastqs.py in Generic mode:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Reading package reference output log file')
        
        #Get check_metadata_ENA reference log file path in package and read
        ext_check_fastqs_log_path = os.path.join(TEST_INFO_PATH, EXT_CHECK_FASTQS_REF_LOG)
        ext_check_fastqs_reference_log = read_reference_log_file(ext_check_fastqs_log_path)
        
        #Get manifest file path in package and Copy file in output_dir
        print('Coping Manifest File Example in External_Dataset_Example Directory')
        ext_manifest_file_path = os.path.join(TEST_INFO_PATH, EXT_MANIFEST_FILE)
        copy_file(EXT_MANIFEST_FILE, ext_manifest_file_path, external_dataset_subdir_path)
                
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('check_fastqs -s Generic -t', EXT_FILTERED_MTD_FILE, '-d', ext_download_subdir_path, '-a', EXT_MANIFEST_FILE, '-p', "'.fq.gz'", '--md5_check', '--plain_text')
        
        #Run program
        ext_filtered_table_path = os.path.join(external_dataset_subdir_path, EXT_FILTERED_MTD_FILE)
        ext_manifest_file_path_def = os.path.join(external_dataset_subdir_path, EXT_MANIFEST_FILE)
        runb4 = subprocess.run(['check_fastqs','-s','Generic','-t', ext_filtered_table_path, '-d', ext_download_subdir_path, '-a', ext_manifest_file_path_def, '-p', '.fq.gz', '--md5_check', '--plain_text'], stdout = subprocess.PIPE, text = True) 

        #Check outputs and run
        check_generic_log_run(CHECK_FASTQS_PATTERN, runb4, EXT_CHECK_FASTQS_LOG, ext_check_fastqs_reference_log, EXT_CHECK_FASTQS_REF_LOG, external_dataset_subdir_path, plain_text_bool)

        ##B5)make_treatment_template Generic mode
        #Messages
        print('\nB5)', rich_text_colored('Testing make_treatment_template.py in Generic mode:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('make_treatment_template -s Generic -i', EXT_MANIFEST_FILE, '-d', ext_download_subdir_path, "-p '.fq.gz' -r1 '_1.fq.gz' -r2 '_2.fq.gz'", '-o', external_dataset_subdir_path, '--plain_text')
        
        #Run program
        runb5_parameters = ['make_treatment_template','-s','Generic', '-i', ext_manifest_file_path_def, '-d', ext_download_subdir_path, '-p', '.fq.gz', '-r1', '_1.fq.gz','-r2', '_2.fq.gz','-o', external_dataset_subdir_path, '--plain_text']
        runb5 = subprocess.run(runb5_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(runb5, external_dataset_subdir_path, RUNB5_OUTPUTFILES, 5, plain_text_bool)

        ##B6)treat_fastqs
        
        #Message
        print('\nB6)', rich_text_colored('Testing treat_fastqs.py:', 'test_color', plain_text_bool))
        
        #Get Filter Table file path in package and Copy file in output_dir
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Treatment Template Example in External_Dataset_Example Directory')
        ext_treatment_table_path = os.path.join(TEST_INFO_PATH, GTREATMENT_FILE_NAME)
        copy_file(GTREATMENT_FILE_NAME, ext_treatment_table_path, external_dataset_subdir_path)
        
        #Get full path and create treated_files subdirectory if it does not exit
        ext_treated_files_subdir_path = os.path.join(external_dataset_subdir_path, TREATMENT_SUBDIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(ext_treated_files_subdir_path) == False:
            print('Creating Treated Files Subdirectory in External_Dataset_Example Directory')
            create_subdirectory(TREATMENT_SUBDIR_NAME, external_dataset_subdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('treat_fastqs -t', GTREATMENT_FILE_NAME, '-i', ext_download_subdir_path, "-p '.fq.gz' -r1 '_1.fq.gz' -r2 '_2.fq.gz'", '-o', ext_treated_files_subdir_path, '--plain_text')
        
        #Run program
        ext_treatment_table_path_def = os.path.join(external_dataset_subdir_path, GTREATMENT_FILE_NAME)
        runb6_parameters = ['treat_fastqs','-t', ext_treatment_table_path_def, '-i', ext_download_subdir_path, '-o', ext_treated_files_subdir_path, '-p', '.fq.gz', '-r1', '_1.fq.gz','-r2', '_2.fq.gz','--plain_text']
        runb6 = subprocess.run(runb6_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(runb6, ext_treated_files_subdir_path, RUNB6_OUTPUTFILES, 5, plain_text_bool)

        ##B7)treat_metadata Generic mode
        
        #Message
        print('\nB7)', rich_text_colored('Testing treat_metadata.py in Generic mode:', 'test_color', plain_text_bool))
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('treat_metadata -s Generic -t', GTREATMENT_FILE_NAME, '-m', EXT_FILTERED_MTD_FILE, "-p '.fq.gz' -r1 '_1.fq.gz' -r2 '_2.fq.gz'", '-o', external_dataset_subdir_path, '--plain_text')
        
        #Run program
        runb7_parameters = ['treat_metadata','-s','Generic','-t', ext_treatment_table_path_def, '-m', ext_filtered_table_path, '-o', external_dataset_subdir_path, '-p', '.fq.gz', '-r1', '_1.fq.gz','-r2', '_2.fq.gz', '--plain_text']
        runb7 = subprocess.run(runb7_parameters, stdout = subprocess.DEVNULL)
        
        #Check outputs and run
        check_generic_run(runb7, external_dataset_subdir_path, RUNB7_OUTPUTFILES, 5, plain_text_bool)

        #C) Meta-analysis Workflow
        #Section header message
        print(rich_text_colored('\nC) Testing Multidatasets Programs:', 'section_header', plain_text_bool))

        ##C1)concat_datasets
        #Messages
        print('\nC1)', rich_text_colored('Testing concat_datasets.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))

        #Get Curated ENA Metadata Table file path in package and Copy file in output_dir
        print('Coping Curated ENA Metadata Table Example in Output Directory')
        curated_ENA_table_path = os.path.join(TEST_INFO_PATH, CURATED_ENA_MTD_FILE_NAME)
        copy_file(CURATED_ENA_MTD_FILE_NAME, curated_ENA_table_path, outputdir_path)
        
        #Get Curated External Metadata Table file path in package and Copy file in output_dir
        print('Coping Curated External Metadata Table Example in Output Directory')
        curated_ext_table_path = os.path.join(TEST_INFO_PATH, CURATED_EXT_MTD_FILE_NAME)
        copy_file(CURATED_EXT_MTD_FILE_NAME, curated_ext_table_path, outputdir_path)
        
        #Get Curated External Metadata Table file path in package and Copy file in output_dir
        print('Coping Variables Dictionary Example in Output Directory')
        var_dict_path = os.path.join(TEST_INFO_PATH, VAR_DICT_FILE_NAME)
        copy_file(VAR_DICT_FILE_NAME, var_dict_path, outputdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('concat_datasets -i', outputdir_path, '-d', VAR_DICT_FILE_NAME, "-op 'example'", '-o', outputdir_path, '--plain_text')
        
        #Run program
        var_dict_path_def = os.path.join(outputdir_path, VAR_DICT_FILE_NAME)
        runc1_parameters = ['concat_datasets','-i', outputdir_path, '-d', var_dict_path_def, '-op', 'example','-o', outputdir_path, '--plain_text']
        runc1 = subprocess.run(runc1_parameters, stdout = subprocess.DEVNULL)
        
        #Check outputs and run
        check_generic_run(runc1, outputdir_path, RUNC1_OUTPUTFILES, 5, plain_text_bool)

        ##C2)check_metadata_values
        
        #Messages
        print('\nC2)', rich_text_colored('Testing check_metadata_values.py:', 'test_color', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Reading package reference output log file')
        
        #Get check_metadata_ENA reference log file path in package and read
        check_mt_vals_log_path = os.path.join(TEST_INFO_PATH, CHECK_CONCAT_MTD_REF_LOG)
        check_mt_vals_reference_log = read_reference_log_file(check_mt_vals_log_path)
                
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('check_metadata_values -t', CONCAT_MTD, '-d', VAR_DICT_FILE_NAME, '--plain_text')
        
        #Run program
        concat_mt_table_path = os.path.join(outputdir_path, CONCAT_MTD)
        runc2 = subprocess.run(['check_metadata_values','-t', concat_mt_table_path, '-d', var_dict_path_def, '--plain_text'], stdout = subprocess.PIPE, text = True)

        #Check outputs and run
        check_generic_log_run(CHECK_MT_VALS_PATTERN, runc2, CHECK_CONCAT_MTD_LOG, check_mt_vals_reference_log, CHECK_CONCAT_MTD_REF_LOG, outputdir_path, plain_text_bool)
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
           
    finally:
        #Print empty line for aesthetic purposes
        print('')
        
if __name__ == '__main__':
    main()
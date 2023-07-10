#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: test_omdctk.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

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

##Set test_info_path
TEST_INFO_PATH = os.path.join(Path(__file__).parent.absolute(), 'test_info')

##1)download_metadata_ENA

STUDY_ACCESSION = 'PRJEB10949'

RUN1_OUTPUTFILES = ['PRJEB10949_ENA_browser.tsv', 'PRJEB10949_mg-toolkit.tsv', 'PRJEB10949_ENA_metadata.tsv']

##2)merge_metadata

ENA_MERGE_COLUMN = 'run_accession'

PUBLICATION_MERGE_COLUMN = 'run_accessions'

PUBLICATION_MTD_FILE_NAME = 'PRJEB10949_publication_example.tsv'

ENA_MTD_FILE_NAME = 'PRJEB10949_ENA_metadata.tsv'

RUN2_OUTPUTFILES = ['PRJEB10949_merged_metadata.tsv']

##3)check_metadata_ENA

MERGED_MTD_FILE_NAME = 'PRJEB10949_merged_metadata.tsv'

CHECK_MTD_REF_LOG = 'PRJEB10949_check_metadata_ENA_reference.log'

CHECK_MTD_LOG = 'PRJEB10949_check_metadata_ENA_stdout.log'

CHECK_MTD_PATTERN = "Runs' Stats:\n"

##4)filter_metadata

FILTER_FILE_NAME = 'PRJEB10949_filterfile_example.tsv'

RUN4_OUTPUTFILES = ['filtered_PRJEB10949_merged_metadata.tsv']

##5)download_fastqs_ENA

DOWNLOAD_SUBDIR_NAME = 'downloads'

FILTERED_MTD_FILE_NAME = 'filtered_PRJEB10949_merged_metadata.tsv'

RUN5_OUTPUTFILES = ['ERR1049851_1.fastq.gz', 'ERR1049851_2.fastq.gz', 'ERR1049852_1.fastq.gz', 'ERR1049852_2.fastq.gz', 'ERR1049853_1.fastq.gz',
                    'ERR1049853_2.fastq.gz', 'ERR1049854_1.fastq.gz', 'ERR1049854_2.fastq.gz', 'ERR1049855_1.fastq.gz', 'ERR1049855_2.fastq.gz',
                    'ERR1049856_1.fastq.gz', 'ERR1049856_2.fastq.gz', 'ERR1049857_1.fastq.gz', 'ERR1049857_2.fastq.gz', 'ERR1049859_1.fastq.gz',
                    'ERR1049859_2.fastq.gz', 'ERR1049861_1.fastq.gz', 'ERR1049861_2.fastq.gz', 'ERR1049911_1.fastq.gz', 'ERR1049911_2.fastq.gz',
                    'ERR1049913_1.fastq.gz', 'ERR1049913_2.fastq.gz', 'ERR1049915_1.fastq.gz', 'ERR1049915_2.fastq.gz', 'ERR1049917_1.fastq.gz',
                    'ERR1049917_2.fastq.gz', 'ERR1049919_1.fastq.gz', 'ERR1049919_2.fastq.gz', 'ERR1049921_1.fastq.gz', 'ERR1049921_2.fastq.gz',
                    'ERR1049923_1.fastq.gz', 'ERR1049923_2.fastq.gz', 'ERR1049924_1.fastq.gz', 'ERR1049924_2.fastq.gz', 'ERR1049925_1.fastq.gz',
                    'ERR1049925_2.fastq.gz', 'ERR1049926_1.fastq.gz', 'ERR1049926_2.fastq.gz', 'ERR1049927_1.fastq.gz', 'ERR1049927_2.fastq.gz',
                    'ERR1049928_1.fastq.gz', 'ERR1049928_2.fastq.gz']

##6)check_fastqs_ENA

CHECK_FASTQS_REF_LOG = 'PRJEB10949_check_fastqs_ENA_reference.log'

CHECK_FASTQS_LOG = 'PRJEB10949_check_fastqs_ENA_stdout.log'

CHECK_FASTQS_PATTERN = 'Main Information:'

##7)make_treatment_template_ENA

RUN7_OUTPUTFILES = ['raw_treatment_template_filtered_PRJEB10949_merged_metadata.tsv']

##8)treat_fastqs

TREATMENT_FILE_NAME = 'treatment_template_filtered_PRJEB10949_merged_metadata_example.tsv'

TREATMENT_SUBDIR_NAME = 'treated_files'

RUN8_OUTPUTFILES = ['ERR1049851_1.fastq.gz', 'ERR1049851_2.fastq.gz', 'ERR1049852_1.fastq.gz', 'ERR1049852_2.fastq.gz', 'ERR1049853_1.fastq.gz',
                    'ERR1049853_2.fastq.gz', 'Ileum-1_1.fastq.gz', 'Ileum-1_2.fastq.gz', 'Ileum-2_1.fastq.gz', 'Ileum-2_2.fastq.gz', 
                    'Ileum-3_1.fastq.gz', 'Ileum-3_2.fastq.gz', 'Ileum-4_1.fastq.gz', 'Ileum-4_2.fastq.gz', 'Ileum-5_1.fastq.gz', 
                    'Ileum-5_2.fastq.gz', 'Ileum-6_1.fastq.gz', 'Ileum-6_2.fastq.gz', 'MAT-1_1.fastq.gz', 'MAT-1_2.fastq.gz', 'MAT-2_1.fastq.gz',
                    'MAT-2_2.fastq.gz', 'MAT-3_1.fastq.gz', 'MAT-3_2.fastq.gz']

##9)treat_metadata_ENA

RUN9_OUTPUTFILES = ['treated_filtered_PRJEB10949_merged_metadata.tsv']

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
        
        #1)download_metadata_ENA
        
        #Messages
        print(rich_text_colored('\n1) Testing download_metadata_ENA.py:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('download_metadata_ENA -p', STUDY_ACCESSION , '-o', outputdir_path, '--plain_text')
        
        #Run program
        run1 = subprocess.run(['download_metadata_ENA','-p', STUDY_ACCESSION, '-o', outputdir_path, '--plain_text'], stdout = subprocess.DEVNULL) 
        
        #Check output and run
        check_generic_run(run1, outputdir_path, RUN1_OUTPUTFILES, 5, plain_text_bool)
        
        #2)merge_metadata
        
        #Messages
        print(rich_text_colored('\n2) Testing merge_metadata.py:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Publication Metadata Table Example in Output Directory')
        
        #Get Publication Table file path in package and Copy file in output_dir
        publication_table_path = os.path.join(TEST_INFO_PATH, PUBLICATION_MTD_FILE_NAME)
        copy_file(PUBLICATION_MTD_FILE_NAME, publication_table_path, outputdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('merge_metadata -e', ENA_MTD_FILE_NAME, '-ec', ENA_MERGE_COLUMN, '-p', PUBLICATION_MTD_FILE_NAME, '-pc', PUBLICATION_MERGE_COLUMN, '-o', outputdir_path, '--plain_text')
        
        #Run program
        publication_table_path_def = os.path.join(outputdir_path, PUBLICATION_MTD_FILE_NAME)
        ena_table_path = os.path.join(outputdir_path, ENA_MTD_FILE_NAME)
        run2_parameters = ['merge_metadata','-e', ena_table_path, '-ec', ENA_MERGE_COLUMN, '-p', publication_table_path_def, '-pc', PUBLICATION_MERGE_COLUMN, '-o', outputdir_path, '--plain_text']
        run2 = subprocess.run(run2_parameters, stdout = subprocess.DEVNULL)
        
        #Check outputs and run
        check_generic_run(run2, outputdir_path, RUN2_OUTPUTFILES, 5, plain_text_bool)
        
        #3)check_metadata_ENA
        
        #Messages
        print(rich_text_colored('\n3) Testing check_metadata_ENA.py:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Reading package reference output log file')
        
        check_metadata_log_path = os.path.join(TEST_INFO_PATH, CHECK_MTD_REF_LOG)
        #Get check_metadata_ENA reference log file path in package and read
        check_metadata_reference_log = read_reference_log_file(check_metadata_log_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('check_metadata_ENA -t', MERGED_MTD_FILE_NAME, '--plain_text')
        
        #Run program
        merged_table_path = os.path.join(outputdir_path, MERGED_MTD_FILE_NAME)
        run3 = subprocess.run(['check_metadata_ENA','-t', merged_table_path, '--plain_text'], stdout = subprocess.PIPE, text = True) 

        #Check outputs and run
        check_generic_log_run(CHECK_MTD_PATTERN, run3, CHECK_MTD_LOG, check_metadata_reference_log, CHECK_MTD_REF_LOG, outputdir_path, plain_text_bool)
        
        #4)filter_metadata
        
        #Messages
        print(rich_text_colored('\n4) Testing filter_metadata.py:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Filter Table Example in Output Directory')
        
        #Get Filter Table file path in package and Copy file in output_dir
        filter_table_path = os.path.join(TEST_INFO_PATH, FILTER_FILE_NAME)
        copy_file(FILTER_FILE_NAME, filter_table_path, outputdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('filter_metadata -t', MERGED_MTD_FILE_NAME, '-f', FILTER_FILE_NAME,'-o', outputdir_path, '--plain_text')
        
        #Run program
        filter_table_path_def = os.path.join(outputdir_path, FILTER_FILE_NAME)
        run4_parameters = ['filter_metadata','-t', merged_table_path, '-f', filter_table_path_def, '-o', outputdir_path, '--plain_text']
        run4 = subprocess.run(run4_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run4, outputdir_path, RUN4_OUTPUTFILES, 5, plain_text_bool)
        
        #5)download_fastqs_ENA
        
        #Message
        print(rich_text_colored('\n5) Testing download_fastqs_ENA.py:', 'section_header', plain_text_bool))
        
        #Get full path and create download subdirectory if it does not exit
        download_subdir_path = os.path.join(outputdir_path, DOWNLOAD_SUBDIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(download_subdir_path) == False:
            print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
            print('Creating Download Subdirectory')
            create_subdirectory(DOWNLOAD_SUBDIR_NAME, outputdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('download_fastqs_ENA -i', FILTERED_MTD_FILE_NAME,'-o', download_subdir_path, '--plain_text')
        
        #Run program
        filtered_table_path = os.path.join(outputdir_path, FILTERED_MTD_FILE_NAME)
        run5_parameters = ['download_fastqs_ENA','-i', filtered_table_path, '-o', download_subdir_path, '--plain_text']
        run5 = subprocess.run(run5_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run5, download_subdir_path, RUN5_OUTPUTFILES, 5, plain_text_bool)
        
        #6)check_fastqs_ENA
        
        #Messages
        print(rich_text_colored('\n6) Testing check_fastqs_ENA.py:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Reading package reference output log file')
        
        #Get check_metadata_ENA reference log file path in package and read
        check_fastqs_log_path = os.path.join(TEST_INFO_PATH, CHECK_FASTQS_REF_LOG)
        check_fastqs_reference_log = read_reference_log_file(check_fastqs_log_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('check_fastqs_ENA -t', FILTERED_MTD_FILE_NAME, '-d', download_subdir_path, '--md5_check', '--plain_text')
        
        #Run program
        run6 = subprocess.run(['check_fastqs_ENA','-t', filtered_table_path, '-d', download_subdir_path, '--md5_check', '--plain_text'], stdout = subprocess.PIPE, text = True) 

        #Check outputs and run
        check_generic_log_run(CHECK_FASTQS_PATTERN, run6, CHECK_FASTQS_LOG, check_fastqs_reference_log, CHECK_FASTQS_REF_LOG, outputdir_path, plain_text_bool)
        
        #7)make_treatment_template_ENA
        
        #Messages
        print(rich_text_colored('\n7) Testing make_treatment_template_ENA.py:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('make_treatment_template_ENA -t', FILTERED_MTD_FILE_NAME, '-d', download_subdir_path, '--extra_sample_columns sample_column', '-o', outputdir_path, '--plain_text')
        
        #Run program
        run7_parameters = ['make_treatment_template_ENA','-t', filtered_table_path, '-d', download_subdir_path, '--extra_sample_columns','sample_column', '-o', outputdir_path, '--plain_text']
        run7 = subprocess.run(run7_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run7, outputdir_path, RUN7_OUTPUTFILES, 5, plain_text_bool)
        
        #8)treat_fastqs
        
        #Message
        print(rich_text_colored('\n8) Testing treat_fastqs.py:', 'section_header', plain_text_bool))
        
        #Get Filter Table file path in package and Copy file in output_dir
        print(rich_text_colored('\nPreparation:', 'general_text', plain_text_bool))
        print('Coping Treatment Template Example in Output Directory')
        treatment_table_path = os.path.join(TEST_INFO_PATH, TREATMENT_FILE_NAME)
        copy_file(TREATMENT_FILE_NAME, treatment_table_path, outputdir_path)
        
        #Get full path and create treated_files subdirectory if it does not exit
        treated_files_subdir_path = os.path.join(outputdir_path, TREATMENT_SUBDIR_NAME)
        ##Check subdir not exits create
        if os.path.isdir(treated_files_subdir_path) == False:
            print('Creating Treated Files Subdirectory')
            create_subdirectory(TREATMENT_SUBDIR_NAME, outputdir_path)
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('treat_fastqs -t', TREATMENT_FILE_NAME, '-i', download_subdir_path, '-o', treated_files_subdir_path, '--plain_text')
        
        #Run program
        treatment_table_path_def = os.path.join(outputdir_path, TREATMENT_FILE_NAME)
        run8_parameters = ['treat_fastqs','-t', treatment_table_path_def, '-i', download_subdir_path, '-o', treated_files_subdir_path, '--plain_text']
        run8 = subprocess.run(run8_parameters, stdout = subprocess.DEVNULL)
        
        #Check output and run
        check_generic_run(run8, treated_files_subdir_path, RUN8_OUTPUTFILES, 5, plain_text_bool)
        
        ##9)treat_metadata_ENA
        
        #Message
        print(rich_text_colored('\n9) Testing treat_metadata_ENA.py:', 'section_header', plain_text_bool))
        
        #Print command (shorter version)
        print(rich_text_colored('\nCommand:', 'general_text', plain_text_bool))
        print('treat_metadata_ENA -t', TREATMENT_FILE_NAME, '-m', FILTERED_MTD_FILE_NAME, '-o', outputdir_path, '--extra_no_warning_columns Run Sample run_accessions run_label', '--plain_text')
        
        #Run program
        run9_parameters = ['treat_metadata_ENA','-t', treatment_table_path_def, '-m', filtered_table_path, '-o', outputdir_path, '--extra_no_warning_columns', 'Run', 'Sample', 'run_accessions', 'run_label', '--plain_text']
        run9 = subprocess.run(run9_parameters, stdout = subprocess.DEVNULL)
        
        #Check outputs and run
        check_generic_run(run9, outputdir_path, RUN9_OUTPUTFILES, 5, plain_text_bool)
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
           
    finally:
        #Print empty line for aesthetic purposes
        print('')
        
if __name__ == '__main__':
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: treat_fastqs.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'treat_fastqs.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, TEMPLATE_FINAL_COLUMNS, VALID_FASTQ_TYPES,
                    VALID_TREATMENTS, OMD_CTK_Exception, program_header, 
                    check_existence_directory_parameter, 
                    check_fastq_PAIRED_patterns, print_list_n_byline,
                    get_list_fastqs_in_directory, rich_text_colored,
                    check_headers, treat_headers_check,
                    check_values, treat_values_check, check_na_in_pandas_dataframe,
                    check_duplicates_in_fastq_names, treat_check_fastq_name_type,
                    check_treatment_for_samples, check_rename_samples, check_merge_samples)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
import pandas as pd
import shutil
import os

#Program Functions
def check_directories_are_not_identical(input_dir, output_dir):
    """
    This function checks if the provided directories are the same or not.

    Parameters
    ----------
    input_dir : str
        Path to the Input Directory.
    output_dir : str
        Path to the Output Directory.

    Raises
    ------
    OMD_CTK_Exception
        If both directories are the same raises an exception.

    Returns
    -------
    None.

    """
    same_dir = os.path.samefile(input_dir, output_dir)
    
    if same_dir == True:
        raise OMD_CTK_Exception('Error! The provided Input and Output directories are the same!')


def check_fastq_file_exits(list_fastq_files, fastqs_in_directory, n_elements):
    """
    This function checks if all files in the list of fastq files are
    present in the list of fastq files of the Input Directory.

    Parameters
    ----------
    list_fastq_files : list
        List with the fastq files from the Treatment Template.
    fastqs_in_directory : list
        List of fastqs detected in the Input Directory.
    n_elements : int
        The number of elements to be printed by line.
        
    Raises
    ------
    OMD_CTK_Exception
        If absent Fastq files are detected raises an exception.

    Returns
    -------
    None.
    
    """
    #Get asymmetric difference between list
    fastqs_difference = list(set(list_fastq_files) - set(fastqs_in_directory))
        
    #If difference is not empty raise exception
    if len(fastqs_difference) > 0:
        #Prepare exception frase
        frase0 = 'Error! Some of the Fastq files in the Treatmente Template are not in the Input Directory!\n'
        frase1 = '\nThe following files were absent:\n'
        frase = [frase0, frase1]
        #Split file names
        splited_list = [fastqs_difference[i:i+n_elements] for i in range(0, len(fastqs_difference), n_elements)]
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
        
    
def cat_files(list_inputfiles, outputfile, color_treatment):
    """
    This function cats files together in the specified file.
    It works with both gz compressed and uncompressed Fastq files.

    Parameters
    ----------
    list_inputfiles : list
        The list of paths for the input gzip files to cat.
    outputfile : str
        The path for the resulting output gzip file.
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
    try:
        with open(outputfile,'wb') as output:
            for file in list_inputfiles:
                with open(file,'rb') as inputs:
                    shutil.copyfileobj(inputs, output)
    except:
        print(rich_text_colored('Error! These files could not be successfully concatenated!', 'exception', color_treatment))
        return True
    else:
        return False


def treat_cat_files(sample_df, fastq_type, input_directory_path, output_directory_path, type_pattern, color_treatment):
    """
    This function treats the provided sample fastq type files to cat them together.

    Parameters
    ----------
    sample_df : pandas dataframe
        The provided sample dataframe.
    fastq_type : str
        The fastq type assigned to the fastq file.
    input_directory_path : str
        Path to the Input Directory.
    output_directory_path : str
        Path to the Output Directory.
    type_pattern : str
        The pattern associated to the provided fastq_type.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    result : bool
        True: There are warnings.
        False: There are no warnings.

    """
    #Get sample_name
    sample_name = list(set(sample_df['sample_name']))[0]   
    
    #Get type info
    sample_df_type = sample_df[sample_df['fastq_type'] == fastq_type]
    
    #Get fastq_name and output_path
    fastq_name = sample_name + type_pattern
    output_path = os.path.join(output_directory_path, fastq_name)
    
    #Get list of input files(and sort by name)
    input_files = list(sample_df_type['fastq_file_name'])
    input_files.sort()
    input_paths = [os.path.join(input_directory_path, file) for file in input_files]
    
    #Print cat configuration
    print(input_files, ' > ', fastq_name)
    
    #Cat the gzip files and return bool value
    result = cat_files(input_paths, output_path, color_treatment)
    return result
    

def merge_mode(sample_df, input_directory_path, output_directory_path, fastq_pattern, R1_pattern, R2_pattern, color_treatment):
    """
    This function merges fastq files from the same sample and prints messages.

    Parameters
    ----------
    sample_df : pandas dataframe
        Pandas dataframe with all the information of the sample files.
    input_directory_path : str
        The provided Input Directory.
    output_directory_path : str
        The provided Output Directory.
    fastq_pattern : str
        The Fastq File Pattern provided.
    R1_pattern : str
        The R1 File Pattern provided.
    R2_pattern : str
        The R2 File Pattern provided.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    result : bool
        True: There are warnings.
        False: There are no warnings.
        
    """
    #Init bool results list
    bools_list = []
    
    #List of fastq_types
    temp_list_fastqtypes = list(sample_df['fastq_type'])
    
    #Get temp number of types
    temp_n_pair1 = temp_list_fastqtypes.count('pair1')
    temp_n_pair2 = temp_list_fastqtypes.count('pair2')
    temp_n_single = temp_list_fastqtypes.count('single')   
       
    #Treat different configurations
    ##Paired files with orphan single fastq files with the same number of files
    if temp_n_pair1 > 1 and temp_n_pair2 > 1 and temp_n_single > 1 and temp_n_pair1 == temp_n_pair2 and temp_n_pair1 == temp_n_single:
        #Show message paired
        print(rich_text_colored('Paired Merging:', 'treatment_mode_color', color_treatment))
        #Generate pair1 file and add bool
        catp1 = treat_cat_files(sample_df, 'pair1', input_directory_path, output_directory_path, R1_pattern, color_treatment)
        bools_list.append(catp1)
        #Generate pair2 file and add bool 
        catp2 = treat_cat_files(sample_df, 'pair2', input_directory_path, output_directory_path, R2_pattern, color_treatment)
        bools_list.append(catp2)
        
        #Show message single
        print(rich_text_colored('Single Merging:', 'treatment_mode_color', color_treatment))
        #Generate pair2 file and add bool 
        cat_s = treat_cat_files(sample_df, 'single', input_directory_path, output_directory_path, fastq_pattern, color_treatment)
        bools_list.append(cat_s)
    #More than one single fastq file but no pair1 or pair2 files
    elif temp_n_pair1 == 0 and temp_n_pair2 == 0 and temp_n_single>1:
        #Show message single
        print(rich_text_colored('Single Merging:', 'treatment_mode_color', color_treatment))
        #Generate pair2 file and add bool 
        cat_s = treat_cat_files(sample_df, 'single', input_directory_path, output_directory_path, fastq_pattern, color_treatment)
        bools_list.append(cat_s)
    #More than one fastq file for pair1 and pair2 with the same sumber of files, no single files 
    elif temp_n_single == 0 and temp_n_pair1 > 1 and temp_n_pair2 > 1 and temp_n_pair1 == temp_n_pair2:
        #Show message paired
        print(rich_text_colored('Paired Merging:', 'treatment_mode_color', color_treatment))
        #Generate pair1 file and add bool
        catp1 = treat_cat_files(sample_df, 'pair1', input_directory_path, output_directory_path, R1_pattern, color_treatment)
        bools_list.append(catp1)
        #Generate pair2 file and add bool 
        catp2 = treat_cat_files(sample_df, 'pair2', input_directory_path, output_directory_path, R2_pattern, color_treatment)
        bools_list.append(catp2)
    #Pass because we have ruled out this possibility previusly with check_merge_samples() function
    else:
        pass
    #Return bool value / Check if theare are any warnings with this sample
    return any(bools_list)

        
def rename_mode(sample_name, fastq_file_name, fastq_type, input_directory_path, output_directory_path, fastq_pattern, R1_pattern, R2_pattern, color_treatment):
    """
    This function copies a fastq file with a new name in the provided 
    Output Directory and prints messages.

    Parameters
    ----------
    sample_name : str
        The sample name of the file that will be used as the new file name.
    fastq_file_name : str
        The file name in the provided Input Directory.
    fastq_type : str
        The type of fastq file [single, pair1, pair2].
    input_directory_path : str
        The provided Input Directory.
    output_directory_path : str
        The provided Output Directory.
    fastq_pattern : str
        The Fastq File Pattern provided.
    R1_pattern : str
        The R1 File Pattern provided.
    R2_pattern : str
        The R2 File Pattern provided.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    result : bool
        True: There are warnings.
        False: There are no warnings.

    """
    #Treat single and paired files and get new names
    if fastq_type == 'pair1':
        fastq_name = sample_name + R1_pattern
    elif fastq_type == 'pair2':
        fastq_name = sample_name + R2_pattern
    else:
        fastq_name = sample_name + fastq_pattern
    
    #Show message
    print([fastq_file_name], ' > ', fastq_name)
    
    #Set input and output paths
    output_path = os.path.join(output_directory_path, fastq_name)
    input_path = os.path.join(input_directory_path, fastq_file_name)
    
    #Try to copy file in output directory
    try:
        shutil.copyfile(input_path, output_path)    
    except:
        print(rich_text_colored('Error! This file could not be successfully renamed in the Output Directory!', 'exception', color_treatment))
        return True
    else:
        return False
    
    
def copy_only_mode(fastq_file_name, input_directory_path, output_directory_path, color_treatment):
    """
    This function copies a file in the provided Output Directory.

    Parameters
    ----------
    fastq_file_name : str
        The file name in the provided Input Directory.
    input_directory_path : str
        The provided Input Directory.
    output_directory_path : str
        The provided Output Directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    result : bool
        True: There are warnings.
        False: There are no warnings.
        
    """
    #Show message
    print([fastq_file_name],' > ',fastq_file_name)
    
    #Set input and output paths
    output_path = os.path.join(output_directory_path, fastq_file_name)
    input_path = os.path.join(input_directory_path, fastq_file_name)
    
    #Try to copy file in output directory
    try:
        shutil.copyfile(input_path, output_path)    
    except:
        print(rich_text_colored('Error! This file could not be successfully copied in the Output Directory!', 'exception', color_treatment))
        return True
    else:
        return False


def treat_files(treatment_df, unique_samples_list, input_directory_path, output_directory_path, fastq_pattern, R1_pattern, R2_pattern, color_treatment):
    """
    This function treats the provided Fastq files based on the Treatment Template.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    unique_samples_list : list
        List of unique sample_names in Treatment Template.
    input_directory_path : str
        The provided Input Directory.
    output_directory_path : str
        The provided Output Directory.
    fastq_pattern : str
        The Fastq File Pattern provided.
    R1_pattern : str
        The R1 File Pattern provided.
    R2_pattern : str
        The R2 File Pattern provided.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.

    """
    #Init warning samples
    warning_samples = []
    
    #A)Show main information pre-treatment
    ##Print information Title
    print(rich_text_colored('\nPre-treatment Information:\n', 'section_header', color_treatment))
    ##Show the number of samples
    print('  o', rich_text_colored('Number of Samples in Treatment Template:', 'subsection', color_treatment), len(unique_samples_list))
    ##Show the number of files
    print('  o', rich_text_colored('Number of Fastq files in Treatment Template:', 'subsection', color_treatment), len(treatment_df['fastq_file_name']))   
    
    #B)Get files for each sample_name and treat
    ##Section header
    print(rich_text_colored('\nTreat Fastqs:', 'section_header', color_treatment))
    ##Treat FASTQS per sample
    for sample in unique_samples_list:
        #Filter table by sample
        temp_sample = treatment_df[treatment_df['sample_name'] == sample]
        #Get treatment(get unique with set function)
        ##NOTE: Up to this point we have check there is only one treatment per sample
        temp_treatment = list(temp_sample['treatment'])[0]
        #List of fastq_types
        temp_list_fastqtypes = list(temp_sample['fastq_type'])
        #Get temp number of types
        temp_n_pair1 = temp_list_fastqtypes.count('pair1')
        temp_n_pair2 = temp_list_fastqtypes.count('pair2')
        temp_n_single = temp_list_fastqtypes.count('single')
        
        #Messages                            
        ##Show sample ID
        print(rich_text_colored('\nSample: ', 'treatment_color', color_treatment), sample)
        ##Show Treatments present
        print(rich_text_colored('Treatment:', 'treatment_color', color_treatment), temp_treatment)
        ##Show number of files
        print(rich_text_colored('Number of Fastq files:', 'general_text', color_treatment), len(temp_sample['fastq_file_name']))
        ##Show files configuration
        print(rich_text_colored('Configuration: ', 'general_text', color_treatment) + ''.join(['Number of pair1(s) = ', str(temp_n_pair1), '; Number of pair2(s) = ', str(temp_n_pair2), '; Number of single(s) = ', str(temp_n_single)]))
        
        #Treat different
        if temp_treatment == 'merge':
            #Merge sample's files
            sample_result = merge_mode(temp_sample, input_directory_path, output_directory_path, fastq_pattern, R1_pattern, R2_pattern, color_treatment)
            #If warnings were detected add sample to list
            if sample_result == True:
                warning_samples.append(sample)
        else:
            if temp_treatment == 'rename':
                print(rich_text_colored('Renaming file(s):', 'treatment_mode_color', color_treatment))
            else:
                print(rich_text_colored('Copying file(s):', 'treatment_mode_color', color_treatment))
            #Apply treatments by files of sample
            #Iter files
            for i,r in temp_sample.iterrows():
                #Get temp file values
                temp_file = r['fastq_file_name']
                temp_file_type = r['fastq_type']
                #Treat modes
                if temp_treatment == 'rename':
                    #Change name for file sample name
                    file_result = rename_mode(sample, temp_file, temp_file_type, input_directory_path, output_directory_path, fastq_pattern, R1_pattern, R2_pattern, color_treatment)
                else:
                    #Copy file sample
                    file_result = copy_only_mode(temp_file, input_directory_path, output_directory_path, color_treatment)
                #If warnings were detected for file add sample to list
                if file_result == True:
                    warning_samples.append(sample)
    
    
    #C)Show main information pre-treatment
    ##Previous steps
    fastqs_in_outputdir = get_list_fastqs_in_directory(output_directory_path, fastq_pattern)
    
    ##Print information Title
    print(rich_text_colored('\nPost-treatment Information:\n', 'section_header', color_treatment))
    ##Show the number of files
    print('  o', rich_text_colored('Number of Fastq files in Output Directory:', 'subsection', color_treatment), len(fastqs_in_outputdir))
    
    #D)Show saved file message
    print(rich_text_colored('\nResulting files saved in:', 'general_text', color_treatment))
    print(output_directory_path)
    
    #E)Treat warnings
    ##Remove duplicates
    warning_samples = list(set(warning_samples))
    ##Show message
    if len(warning_samples) > 0:
        print(rich_text_colored('Error! Some of the samples in the Treatment Template showed warnings when trying to treat their files!', 'program_warning', color_treatment))
        print('\nThe following samples present warnings:')
        print_list_n_byline(warning_samples, 5)
        
        
#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter treatment_file
    parser.add_argument(
            '-t','--treatment_template', 
            action = 'store',
            required = True ,
            help = 'Treatment Template [Expected sep=TABS]. Indicate the path to the Treatment Template file. See Documentation for more information and format details.'
    )
    ##Parameter input_directory
    parser.add_argument(
            '-i','--input_directory', 
            action = 'store',
            required = True,
            help = 'Input Directory. Indicate the path to the Input Directory with the Fastq files to treat.'
    )
    ##Parameter output_directory
    parser.add_argument(
            '-o','--output_directory', 
            action = 'store',
            required = True,
            help = 'Output Directory. Indicate the path to the Output Directory to save the resulting treated Fastq files.'
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
    treatment_file_path = args.treatment_template
    inputdir_path = args.input_directory
    outputdir_path = args.output_directory
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
        #0)Initial checks
        
        #Check that provided files directory exist
        check_existence_directory_parameter(inputdir_path, 'Input', '--input_directory')
        
        #Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #Check that directories are not the same
        check_directories_are_not_identical(inputdir_path, outputdir_path)
        
        #Check PAIRED Fastq patterns
        check_fastq_PAIRED_patterns(fastq_pattern, r1_files_pattern, r2_files_pattern)
        
        #Try to get list of fastqs in the provided directory
        fastqs_in_directory = get_list_fastqs_in_directory(inputdir_path, fastq_pattern)
        
        #1)Try to load treatment file
        
        #Section header message
        print(rich_text_colored('\nLoading Files:\n', 'section_header', plain_text_bool))
        
        #Try to load Treatment Template file as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('Treatment Template file:', 'general_text', plain_text_bool))
        print(treatment_file_path)
        ##Load metadata file as pandas df
        treatment_table = pd.read_csv(treatment_file_path, sep='\t')

        #2)Checks related to the treatment file
        
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
                
        #3)Checks per file
        
        ##Check that all files in Treatment Template have matching file names and fastq types
        treat_check_fastq_name_type(treatment_table, fastq_pattern, r1_files_pattern, r2_files_pattern, 5)
        
        ##Check that all files in Treatment Template are in the Input Directory
        list_fastq_files = list(treatment_table['fastq_file_name'])
        check_fastq_file_exits(list_fastq_files, fastqs_in_directory, 5)

        #4)Checks per sample
        
        ##Get unique sample names and sort
        unique_sample_names = list(set(treatment_table['sample_name']))
        unique_sample_names.sort()
        
        ##Check mixed treatments per sample
        check_treatment_for_samples(treatment_table, unique_sample_names, 5)
        
        ##Check rename mode files configurations per sample
        check_rename_samples(treatment_table, unique_sample_names, plain_text_bool)
        
        ##Check merge mode files configurations per sample
        check_merge_samples(treatment_table, unique_sample_names, plain_text_bool)
        
        #5)Treat Fastq files

        ##Treat Fastq files
        treat_files(treatment_table, unique_sample_names, inputdir_path, outputdir_path, fastq_pattern, r1_files_pattern, r2_files_pattern, plain_text_bool)
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
           
    finally:
        #Print empty line for aesthetic purposes
        print('')
                                                                                                   
if __name__ == '__main__':
    main()
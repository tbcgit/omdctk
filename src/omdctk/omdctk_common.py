#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: omdtk_common.py
@author: Samuel Piquer-Esteban
@date:  12 Mar 2024 

"""

#Import third-party modules
from termcolor import colored
import pandas as pd
import os

#Common constants

##Package information
DATE = '12 Mar 2024'
VERSION = 'v1.1.0'

##Package Figlet Ascii (small)
FIGLET = (" ___  __  __ ___      ___              _   _",
          "/ _ \|  \/  |   \    / __|  _ _ _ __ _| |_(_)___ _ _",
          "|(_)|| |\/| | |) |   |(_| || | '_/ _` |  _| / _ \ ' \ ",
          "\___/|_|  |_|___/    \___\_,_|_| \__,_|\__|_\___/_||_|",
          "              _____         _ _   _ _",
          "             |_   _|__  ___| | |_(_) |_",
          "               | |/ _ \/ _ \ | / / |  _|",
          "               |_|\___/\___/_|_\_\_|\__|")

##ENA FASTQ URLs Column Names
ENA_FASTQ_URLS_COLUMNS = ['fastq_ftp', 'fastq_aspera', 'fastq_galaxy',
                      'submitted_ftp', 'submitted_aspera', 'submitted_galaxy']

##Treatment Template check lists
TEMPLATE_FINAL_COLUMNS = ['sample_name', 'fastq_file_name', 'fastq_type', 'treatment']

VALID_FASTQ_TYPES = ['pair1', 'pair2', 'single']

VALID_TREATMENTS = ['merge', 'copy', 'rename']

##Variables Dictionary check lists
VALID_REQUIREDNESS = ['required', 'optional']

##Generic manifest table file headers
VALID_MANIFEST_COLUMNS = ['file_name', 'sample_name','file_md5']

##Colors Dictionary
COLORS_DIC = {'acceptable':'green', 'legend_warning':'magenta', 'dangerous':'red',
              'section_header':'cyan', 'subsection':'yellow', 'general_text':'yellow',
              'subsection2':'magenta', 'subsection3':'green', 'exception':'red', 
              'column_color':'blue', 'pattern_color':'blue', 'filter_color':'green',
              'program_warning':'red', 'program_warning2':'magenta',
              'due_to_header':'magenta', 'you_should_header':'green',
              'true_color':'green', 'false_color':'red', 'check_color':'cyan',
              'treatment_mode_color':'blue', 'treatment_color':'green',
              'figlet_border_color':'blue', 'figlet_package_color':'yellow', 
              'figlet_version_color':'yellow', 'figlet_program_color':'green', 
              'exception':'red','merge_mode_color':'blue', 'merge_table_color':'green',
              'dataset_color':'green', 'variable_color':'yellow', 'test_color':'yellow',
              'dict_column_color':'blue','dict_column_color2':'yellow'}

#Common classes
class OMD_CTK_Exception(Exception):
    pass

#Common function 
def rich_text_colored(text_str, dic_category, color_treatment):
    """
    This function prints the provided text in the color associated to the
    COLORS_DIC category indicated using the rich package.

    Parameters
    ----------
    text_str : str
        The provided text to be color printed.
    dic_category : str
        The category in the COLORS_DIC to get the associated color style.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    result : str
        The resulting string.

    """
    if (dic_category in COLORS_DIC) and (color_treatment == False):
        result = colored(text_str, COLORS_DIC[dic_category], attrs = ['bold'])
    else:
        result = text_str
    return result


def program_header(border_character, program, line_len, start_gap, color_treatment):
    """
    This function prints the header information of the program.

    Parameters
    ----------
    border_character : str
        Border character to construct the outer rectangle.
    program : str
        Name of the program.
    line_len : int
        Total Length of the lines to be printed.
    start_gap : int
        Gap space to leave on the left.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.

    """
    #Print top frame
    print(rich_text_colored(border_character*line_len, 'figlet_border_color', color_treatment))
    print(rich_text_colored(border_character*2 + ' '*(line_len-4) + border_character*2, 'figlet_border_color', color_treatment))
    #Print figlet package name
    for i in FIGLET:
        print(rich_text_colored(border_character*2 + ' '*start_gap, 'figlet_border_color', color_treatment), 
              rich_text_colored(i, 'figlet_package_color', color_treatment), 
              ' '*(line_len - (start_gap+7+len(i))), 
              rich_text_colored(border_character*2, 'figlet_border_color', color_treatment))
    #Print space
    print(rich_text_colored(border_character*2 + ' '*(line_len-4) + border_character*2, 'figlet_border_color', color_treatment))
    #Print program name
    print(rich_text_colored(border_character*2 + ' '*start_gap, 'figlet_border_color', color_treatment), 
          rich_text_colored(program, 'figlet_program_color', color_treatment), 
          ' '*(line_len - (start_gap+7+len(program))), 
          rich_text_colored(border_character*2, 'figlet_border_color', color_treatment))
    #Print package information
    version_frase = ' * {} - {} *'.format(VERSION, DATE)
    print(rich_text_colored(border_character*2 + ' '*start_gap, 'figlet_border_color', color_treatment), 
          rich_text_colored(version_frase, 'figlet_version_color', color_treatment), 
          ' '*(line_len - (start_gap+7+len(version_frase))), 
          rich_text_colored(border_character*2, 'figlet_border_color', color_treatment))
    #Print bottom frame
    print(rich_text_colored(border_character*2 + ' '*(line_len-4) + border_character*2, 'figlet_border_color', color_treatment))
    print(rich_text_colored(border_character*line_len, 'figlet_border_color', color_treatment))


def show_advise_legend(color_treatment):
    """
    This function sets and prints the legend for warning messages.

    Parameters
    ----------
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.

    """
    #Set legend strs
    legend = (rich_text_colored('\nLegend for Advise Messages:', 'section_header', color_treatment), '',
              rich_text_colored('[Acceptable]', 'acceptable', color_treatment) + ' You should be able to continue without much trouble.', '',
              rich_text_colored(' [Dangerous]', 'dangerous', color_treatment) + ' You must be extremely careful. There could be major problems.', '', 
              rich_text_colored('   [Warning]', 'legend_warning', color_treatment) + ' You should be able to continue with some effort,\n             but it could get complicated or even dangerous. Be extra careful. ')
    #Print legend
    for i in legend: print(i)


def check_existence_directory_parameter(dir_path, dir_type, parameter):
    """
    This function will check if the provided directory path exist.

    Parameters
    ----------
    dir_path : str
        Path for the provided directory.
    dir_type : str
        Type of directory (normally, "input" or "output").
    parameter : str
        The name of the parameter.

    Raises
    ------
    OMD_CTK_Exception
        If the provided directory does not exit raises an exception.

    Returns
    -------
    None.

    """
    if dir_path is None:
        pass
    else:
        if os.path.isdir(dir_path)==False:
            frase = ''.join(['Error! The provided ', dir_type, ' Directory path is not a directory!\n Check the provided ', parameter, ' parameter!'])
            raise OMD_CTK_Exception(frase)


def treat_duplicated_outfiles(directory, file_extension, file_name):
    """
    This function generates a unique outputfile name for the provided
    file_name and output directory.

    Parameters
    ----------
    directory : str
        The output directory to check.
    file_extension : str
        The extension of the provided output file_name.
    file_name : str
        The name of the provided output file.

    Returns
    -------
    str
        If file_name is not in directory, returns the same provided file_name.
        Else it will return a new_file_name.

    """
    #Get files in directory
    files = os.listdir(directory)
    #Get copy for file name
    new_file_name = file_name
    #IF/ELSE
    if new_file_name in files:
        #Count
        counter = 1
        #While in list try new number combination
        while new_file_name in files:
            new_file_name = ''.join([file_name.split(file_extension)[0],'(',str(counter),')',file_extension])
            counter += 1
        #Return final file name
        return new_file_name
    else:
        #If file not in directory return the given file_name
        return file_name
       

def treat_output_directory_parameter_outfiles(file_name, outputdir_path):
     """
     This function will treat the output_directory parameter and return
     the full path for the output file_name provided.

     Parameters
     ----------
     file_name : str
         Name for the output file.
     outputdir_path : str
         Output directory path.

     Returns
     -------
     outputfile_path : str
         Path for the output file.

     """
     ##If not path provided create file in current/ else create file in provided path
     if outputdir_path is None:
         #Check if file_name exits in current directory
         #If exits get new name
         out_name = treat_duplicated_outfiles(os.getcwd(), '.tsv', file_name)
         #Get final full path
         outputfile_path = os.path.join(os.getcwd(), out_name)
     else:
         #Check if file_name exits in provided output directory
         #If exits get new name
         out_name = treat_duplicated_outfiles(outputdir_path, '.tsv', file_name)
         #Get final full path
         outputfile_path = os.path.join(outputdir_path, out_name)
     return outputfile_path


def treat_output_directory_parameter(outputdir_path):
     """
     This function will treat the output_directory parameter and return
     the path for the output directory.

     Parameters
     ----------
     outputdir_path : str
         Output directory path.

     Returns
     -------
     output_path : str
         Path for the output directory.

     """
     ##If not path provided create files in current/ else create files in provided path
     if outputdir_path is None:
         #Get final full path
         output_path = os.getcwd()
     else:
         #Get final full path
         output_path = outputdir_path
     return output_path
        
    
def check_headers(list_headers, pandas_df):
    """
    This function checks if a list of headers is in a pandas dataframe.

    Parameters
    ----------
    list_headers : str / list
        The list of headers to check. 
        It can be either a list of one or multiple headers,
        or a str single header. 
        
    pandas_df : pandas dataframe
        The pandas dataframe to check the headers on.

    Returns
    -------
    result : bool
        True: if all provided headers are present in the dataframe.
        False: if any provided header is not present in the dataframe.
        
    """
    result = pd.Series(list_headers).isin(pandas_df.columns).all()
    return result
    
    
def treat_headers_check(headers_used, headers_check, table, frase0, frase1, color_treatment):
    """
    This function treats the results for headers checks on the table.

    Parameters
    ----------
    headers_used : list
        List with the names of the column names to be used.
    headers_check : bool
        The result from check_headers().
        True: if all provided headers are present in the dataframe.
        False: if any provided header is not present in the dataframe.
    table : pandas dataframe
        The pandas dataframe to check the headers on.
    frase0 : str
        The initial frase of the exception which is check dependent.
    frase1 : str
        The headers name frase which is check dependent.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If any of the column names to be used are not present in the 
        provided table raises an exception.

    Returns
    -------
    None.

    """
    if headers_check == False:
        #Prepare exception frase
        frase = [frase0, frase1]
        for i in headers_used:
            if check_headers(i, table) == True:
                temp_line = ' '.join(['\n-', rich_text_colored(i, 'check_color', color_treatment), 
                                      rich_text_colored('Found!', 'true_color', color_treatment)])
            else:
                temp_line = ' '.join(['\n-', rich_text_colored(i, 'check_color', color_treatment), 
                                      rich_text_colored('Not Found!', 'false_color', color_treatment)])
            frase.append(temp_line)
        #Raise exception
        raise OMD_CTK_Exception(' '.join(frase))
    
def check_generic_column_in_metadata(generic_column, metadata_df, parameter):
    """
    This function checks the presence of the provided generic_column in a
    generic metadata table.

    Parameters
    ----------
    generic_column : str
        Generic column in Metadata Table. 
    metadata_df : pandas dataframe
        The provided generic metadata dataframe.
    parameter :str
        The associated parameter.

    Raises
    ------
    OMD_CTK_Exception
        If provided generic_column is absent raises and exception.

    Returns
    -------
    None.

    """
    check = check_headers(generic_column, metadata_df)
    if check == False:
        frase = ''.join(['Error! The provided column is not in the Metadata Table!\n Check your metadata file and --',parameter,' parameter!'])
        #Raise exception
        raise OMD_CTK_Exception(frase)


def print_list_n_byline(list_to_print, n_elements):
    """
    This function prints elements of the provided list 
    n_elements by line.

    Parameters
    ----------
    list_to_print : list
        The provided list to print.
    n_elements : int
        The number of elements to be printed by line.

    Returns
    -------
    None.

    """
    #Generate sublist of n_elements len
    splited_list = [list_to_print[i:i+n_elements] for i in range(0,len(list_to_print),n_elements)]
    #Iter list of list / Each list will be a line
    for i in splited_list:
        #Format depending on if it is the last line
        if i == splited_list[-1]:
            print(*i, sep = ', ', end = '\n')
        else:
            print(*i, sep = ', ', end = ',\n')


def get_urls_from_ENA_column(metadata_df, column):
    """
    This function gets the urls from 
    the indicated ENA metadata column.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    column : str
        The provided ENA metadata column from which to obtain the urls.
    
    Raises
    ------
    OMD_CTK_Exception
        If all cells in the provided ENA download column are empty raises
        an exception.

    Returns
    -------
    urls : list
        The list of urls from the provided ENA metadata column.
        
    """
    #Get urls for the ena_download_column
    ##Set init variables
    urls = []
    check_nans = list(pd.isna(metadata_df[column]))
    #If all nans(empty colum) pass/ Else continue
    if all(check_nans):
        raise OMD_CTK_Exception('Error! All cells are empty in the provided ENA Download Column!\n Check your metadata file!')
    else:
        #Get column content
        fastqs4sample = list(metadata_df[column])
        ##Loop fastqs4sample and process
        for i in fastqs4sample:
            urls += i.split(';')
        ##Remove empty values in list
        urls = list(filter(None, urls))
    return urls


def check_values(list_values, pandas_column):
    """
    This function checks if all values in a pandas_colum are in a list of values.

    Parameters
    ----------
    list_values : list
        The provided list of values used as reference.
    pandas_column : pandas dataframe column
        The pandas dataframe column to check.

    Returns
    -------
    result : bool
        True: all values in the pandas_colum are in the list.
        False: not all values in the pandas_colum are in the list.
        
    """
    result = pandas_column.isin(list_values).all()
    return result


def check_values_inv(pandas_column, list_values):
    """
    This function checks if all values in a list are in a pandas_colum.

    Parameters
    ----------
    pandas_column : pandas dataframe column
        The provided pandas dataframe column used as reference.
    list_values : list
        The provided list of values to check.

    Returns
    -------
    result : bool
        True: all values in the list are in the pandas_colum.
        False: not all values in the list are in the pandas_colum.
        
    """
    result = all(elem in list(pandas_column) for elem in list_values)
    
    return result


def treat_values_check(expected_values, values_check, frase0, frase1, color_treatment):
    """
    This function treats the results for values checks.

    Parameters
    ----------
    expected_values : list
        List with the expected values.
    values_check : bool
        The result from check_values().
        True: if all column values are a valid value in expected_values.
        False: if any column value is not a valid value in expected_values.
    frase0 : str
        The initial frase of the exception which is check dependent.
    frase1 : str
        The expected_values name frase which is check dependent.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If any of the column values are not present in the expected_values
        raises an exception.

    Returns
    -------
    None.

    """
    if values_check == False:
        #Prepare exception frase
        frase = [frase0, frase1]
        for i in expected_values:
            temp_line = ' '.join(['\n-', rich_text_colored(i, 'check_color', color_treatment)])
            frase.append(temp_line)
        #Raise exception
        raise OMD_CTK_Exception(' '.join(frase))


def get_list_fastqs_in_directory(directory, fastq_pattern):
    """
    This function gets a list with the fastq files 
    of the provided directory for the provided fastq_pattern.

    Parameters
    ----------
    directory : str
        Path to the directory.
    fastq_pattern : str
        Extension pattern to recognize a Fastq file.
        
    Raises
    ------
    OMD_CTK_Exception
        If there are no fastq files detected raises an exception.

    Returns
    -------
    fastq_files : list
        List of the fastq files present in the directory.
    
    """
    #Get files in directory
    files = os.listdir(directory)
    
    #Filter to retain only fastq.gz files
    fastq_files = [x for x in files if x.endswith(fastq_pattern)]
    
    #Check if fastq files were detected
    if len(fastq_files) == 0:
        frase0 = 'Error! There are no Fastq files with pattern "'
        frase1 = '" in the provided Fastqs Directory!'
        raise OMD_CTK_Exception(''.join([frase0, fastq_pattern, frase1]))
    
    return fastq_files


def check_fastq_PAIRED_patterns(fastq_pattern, r1_files_pattern, r2_files_pattern):
    """
    This function checks if the provided Fastq PAIRED patterns are compatible
    with the provided Fastq File Pattern.

    Parameters
    ----------
    fastq_pattern : str
        The provided Fastq File Pattern.
    r1_files_pattern : str
        The provided R1 File Pattern.
    r2_files_pattern : str
        The provided R2 File Pattern.

    Raises
    ------
    OMD_CTK_Exception
        If the patterns are not compatible raises an exception.

    Returns
    -------
    None.

    """
    if r1_files_pattern.endswith(fastq_pattern) == False or r2_files_pattern.endswith(fastq_pattern) == False:
        raise OMD_CTK_Exception('Error! The provided Fastq Pattern does not match with some of the PAIRED files patterns!\n Check the provided patterns parameters!')
    elif r1_files_pattern == r2_files_pattern:
        raise OMD_CTK_Exception('Error! The provided PAIRED files patterns are identical!\n Check the provided patterns parameters!')
    elif r1_files_pattern == fastq_pattern or r2_files_pattern == fastq_pattern:
        raise OMD_CTK_Exception('Error! The provided Fastq Pattern is identical to some of the PAIRED files patterns!\n Check the provided patterns parameters!')


def metadata_files_main_information(metadata_df, urls, fastqs_in_directory, color_treatment):
    """
    This function prints the main information for the metadata and downloaded
    Fastq files. Intented to be used in check_fastqs.py and 
    make_treatment_template.py.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
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
    None.

    """    
    ##Show the number of run_accessions
    print(rich_text_colored('\n1) Number of run_accessions:', 'subsection', color_treatment), len(metadata_df.index))
    ##Show the number of unique sample_accessions
    print(rich_text_colored('\n2) Number of unique sample_accessions:', 'subsection', color_treatment), metadata_df['sample_accession'].nunique(dropna = False))
    ##Show the number of unique sample_alias
    print(rich_text_colored('\n3) Number of unique sample_alias:', 'subsection', color_treatment), metadata_df['sample_alias'].nunique(dropna = False))
    ##Show the number of appereances per library_layout
    print(rich_text_colored('\n4) Appearances per library_layout:', 'subsection', color_treatment))
    print(metadata_df['library_layout'].value_counts(dropna = False).rename_axis('unique_values').reset_index(name = 'counts').to_string(index = False))
    ##Show the number of detected URLS to download
    print(rich_text_colored('\n5) Number of URLs expected to be downloaded:', 'subsection', color_treatment), len(urls))
    ##Show the number of files downloaded ending in 'fastq.gz' for files_directory
    print(rich_text_colored('\n6) Number of Fastq files in the provided Fastqs Directory:', 'subsection', color_treatment), len(fastqs_in_directory))


def generic_files_main_information(metadata_df, generic_merge_col, manifest_df, fastqs_in_directory, color_treatment):
    """
    This function print the main information for the generic files and downloaded
    Fastq files. Intented to be used in check_fastqs.py and 
    make_treatment_template.py.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided generic metadata dataframe.
    generic_merge_col : str
        Merge column in Generic Metadata Table. Must be compatible with 
        'sample_name' column from Manifest Table.
    manifest_df : pandas dataframe
        The provided manifest dataframe.
    fastqs_in_directory : list
        List of fastqs detected in the provided directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.

    """
    ##Show the number of rows in generic metadata table
    print(rich_text_colored('\n1) Number of rows in Metadata Table:', 'subsection', color_treatment), len(metadata_df.index))
    ##Show the number of unique sample_accessions
    print(rich_text_colored('\n2) Number of unique samples for provided generic_merge_column in Metadata Table:', 'subsection', color_treatment), metadata_df[generic_merge_col].nunique(dropna = False))
    ##Show the number of fastq files in manifest
    print(rich_text_colored('\n3) Number of Fastq files in the provided Manifest Table:', 'subsection', color_treatment), len(manifest_df.index))
    ##Show the number of unique sample_accessions
    print(rich_text_colored('\n4) Number of unique samples for "sample_name" column in Manifest Table:', 'subsection', color_treatment), manifest_df['sample_name'].nunique(dropna = False))
    ##Show the number of files downloaded ending in 'fastq.gz' for files_directory
    print(rich_text_colored('\n5) Number of Fastq files in the provided Fastqs Directory:', 'subsection', color_treatment), len(fastqs_in_directory))


def check_na_in_pandas_dataframe(pandas_column, frase):
    """
    This function checks if a pandas dataframe column contais any NA.

    Parameters
    ----------
    pandas_column : pandas dataframe column
        Indicate the pandas dataframe column to be checked.
    frase : str
        Error message to be shown if NAs are detected.

    Raises
    ------
    OMD_CTK_Exception
        If any NAs are detected raises an exception.

    Returns
    -------
    None.

    """
    if any(pd.isna(pandas_column)):
        raise OMD_CTK_Exception(frase)


def check_duplicates_in_fastq_names(treatment_df, n_elements):
    """
    This function checks if there are duplicates in the Fastq names
    column of the Treatment Template.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    n_elements : int
        The number of elements to be printed by line.

    Raises
    ------
    OMD_CTK_Exception
        If duplicates are detected raises an exception.

    Returns
    -------
    None.

    """
    #Get duplicates df
    duplicates_df = treatment_df[treatment_df.duplicated(subset=['fastq_file_name'])]
    
    #Treat duplicates
    if len(duplicates_df) > 0:
        #Get list of duplicates
        duplicates = list(duplicates_df['fastq_file_name'])
        #Construct frase
        frase0 = 'Error! Some of the values of the "fastq_file_name" column in the Treatment Template are duplicates!\n Check your treatment file!\n'
        frase1 = '\nThe following files present duplicates:\n'
        frase = [frase0, frase1]
        #Split file names
        splited_list = [duplicates[i:i+n_elements] for i in range(0, len(duplicates), n_elements)]
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


def check_fastq_name_type(fastq_file_name, fastq_type, fastq_pattern, R1_pattern, R2_pattern):
    """
    This function checks if a fastq_file_name ends with the expected pattern
    for its fastq_type. It processes only one input at a time. 

    Parameters
    ----------
    fastq_file_name : str
        The fastq file to check.
    fastq_type : str
        The fastq type assigned to the fastq file.
    fastq_pattern : str
        The Fastq File Pattern provided.
    R1_pattern : str
        The R1 File Pattern provided.
    R2_pattern : str
        The R2 File Pattern provided.

    Returns
    -------
    bool_value : bool
        True: file name and patterns match.
        False: file name and patterns do not match.
    
    """
    if fastq_type == 'pair1' and fastq_file_name.endswith(R1_pattern):
        bool_value = True
    elif fastq_type == 'pair2' and fastq_file_name.endswith(R2_pattern):
        bool_value = True
    elif fastq_type == 'single' and fastq_file_name.endswith(fastq_pattern) and not (fastq_file_name.endswith(R1_pattern) or fastq_file_name.endswith(R2_pattern)):
        bool_value = True
    else:
        bool_value = False
    ##Return results
    return bool_value


def treat_check_fastq_name_type(treatment_df, fastq_pattern, R1_pattern, R2_pattern, n_elements):
    """
    This function treats the check_fastq_name_type results per file.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    fastq_pattern : str
        The Fastq File Pattern provided.
    R1_pattern : str
        The R1 File Pattern provided.
    R2_pattern : str
        The R2 File Pattern provided.
    n_elements : int
        The number of elements to be printed by line.

    Raises
    ------
    OMD_CTK_Exception
        If mismatches between fastq file names and fastq types are detected
        raises an exception.

    Returns
    -------
    None.

    """
    #Init warning files list
    warnings_list = []
    
    #Check each file
    for i,r in treatment_df.iterrows():
        #Get temp values
        temp_file_name = r['fastq_file_name']
        #Get temp type
        temp_type = r['fastq_type']
        #Get temp check
        temp_check = check_fastq_name_type(temp_file_name, temp_type, fastq_pattern, R1_pattern, R2_pattern)
        if temp_check == False:
            warnings_list.append(temp_file_name)
    
    #If warnings_list is not empty raise exception
    if len(warnings_list) > 0:
        #Prepare exception frase
        frase0 = 'Error! Some of the Fastq file names in the Treatmente Template do not match the provided file patterns!\n'
        frase1 = '\nThe following files present mismatches:\n'
        frase = [frase0, frase1]
        #Split file names
        splited_list = [warnings_list[i:i+n_elements] for i in range(0, len(warnings_list), n_elements)]
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


def check_treatment_for_samples(treatment_df, unique_samples_list, n_elements):
    """
    This function checks if the treatment is the same for the files 
    of each sample.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    unique_samples_list : list
        List of unique sample_names in Treatment Template.
    n_elements : int
        The number of elements to be printed by line.
    
    Raises
    ------
    OMD_CTK_Exception
        If mixed treatments per sample are detected raises an exception.

    Returns
    -------
    None.
    
    """
    #Init results list
    warnings_list = []
    
    #Iter by sample
    for sample in unique_samples_list:
        #Filter table by sample
        temp_sample = treatment_df[treatment_df['sample_name'] == sample]
        #Unique treatments with set
        temp_set_treatment = set(temp_sample['treatment'])
        #Check unique values len for treatment
        if len(temp_set_treatment)!=1:
            warnings_list.append(sample)
    
    #If warnings_list is not empty raise exception
    if len(warnings_list) > 0:
        #Prepare exception frase
        frase0 = 'Error! Some of the samples in the Treatment Template show mixed treatments!\n'
        frase1 = '\nThe following samples present mixed treatments:\n'
        frase = [frase0, frase1]
        #Split file names
        splited_list = [warnings_list[i:i+n_elements] for i in range(0, len(warnings_list), n_elements)]
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


def check_rename_samples(treatment_df, unique_samples_list, color_treatment):
    """
    This function checks if the fastq files configurations are
    acceptable for the files of each sample in rename mode.
    Only rename mode samples will be checked.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    unique_samples_list : list
        List of unique sample_names in Treatment Template.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text
    
    Raises
    ------
    OMD_CTK_Exception
        If incompatible fastq files configurations per sample in rename
        mode are detected raises an exception.

    Returns
    -------
    None.
    
    """
    #Init results list
    check_frases = []
    
    #Iter by sample
    for sample in unique_samples_list:
        #Filter table by sample
        temp_sample = treatment_df[treatment_df['sample_name'] == sample]
        #Unique treatments with set
        temp_set_treatment = set(temp_sample['treatment'])
        #List of fastq_types
        temp_list_fastqtypes = list(temp_sample['fastq_type'])
        #Get temp number of types
        temp_n_pair1 = temp_list_fastqtypes.count('pair1')
        temp_n_pair2 = temp_list_fastqtypes.count('pair2')
        temp_n_single = temp_list_fastqtypes.count('single')
        #Treat only rename treatment samples
        if 'rename' not in temp_set_treatment:
            pass
        else:
            #Paired files with orphan single fastq file
            if temp_n_pair1 == 1 and temp_n_pair2 == 1 and temp_n_single == 1:
                pass
            #A unique single fastq file
            elif temp_n_pair1==0 and temp_n_pair2==0 and temp_n_single==1:
                pass
            #A pair of PAIRED fastq files
            elif temp_n_pair1==1 and temp_n_pair2==1 and temp_n_single==0:
                pass
            #Anything else is not an acceptable rename configuration
            else:
                print_strn = '\n- ' + rich_text_colored('Sample Name: ', 'check_color', color_treatment) + sample
                print_strn = print_strn + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + ''.join(['Number of pair1(s) = ', str(temp_n_pair1), '; Number of pair2(s) = ', str(temp_n_pair2), '; Number of single(s) = ', str(temp_n_single)])
                check_frases.append(print_strn)
    
    #If check_frases is not empty raise exception
    if len(check_frases) > 0:
        #Prepare exception frase
        frase0 = 'Error! Some of the samples in the Treatment Template have incompatible file configurations for rename mode!\n'
        frase1 = '\nThe permitted configurations are:'
        frase2 = '\n- ' + rich_text_colored('A pair of PAIRED Fastq files with a unique SINGLE Fastq file', 'acceptable', color_treatment) + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + 'Number of pair1(s) = 1; Number of pair2(s) = 1; Number of single(s) = 1'
        frase3 = '\n- ' + rich_text_colored('A pair of PAIRED Fastq files', 'acceptable', color_treatment) + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + 'Number of pair1(s) = 1; Number of pair2(s) = 1; Number of single(s) = 0'
        frase4 = '\n- ' + rich_text_colored('A unique SINGLE Fastq file', 'acceptable', color_treatment) + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + 'Number of pair1(s) = 0; Number of pair2(s) = 0; Number of single(s) = 1'
        frase5 = '\n\nThe following samples present incompatibilities:'
        frase = [frase0, frase1, frase2, frase3, frase4, frase5]
        final_frase = frase + check_frases
        #Raise exception
        raise OMD_CTK_Exception(' '.join(final_frase))


def check_merge_samples(treatment_df, unique_samples_list, color_treatment):
    """ 
    This function checks if the fastq files configurations are
    acceptable for the files of each sample in merge mode.
    Only merge mode samples will be checked.

    Parameters
    ----------
    treatment_df : pandas dataframe
        The provided treatment dataframe.
    unique_samples_list : list
        List of unique sample_names in Treatment Template.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text
    
    Raises
    ------
    OMD_CTK_Exception
        If incompatible fastq files configurations per sample in merge
        mode are detected raises an exception.

    Returns
    -------
    None.
    
    """
    #Init results list
    check_frases = []
    
    #Iter by sample
    for sample in unique_samples_list:
        #Filter table by sample
        temp_sample = treatment_df[treatment_df['sample_name'] == sample]
        #Unique treatments with set
        temp_set_treatment = set(temp_sample['treatment'])
        #List of fastq_types
        temp_list_fastqtypes = list(temp_sample['fastq_type'])
        #Get temp number of types
        temp_n_pair1 = temp_list_fastqtypes.count('pair1')
        temp_n_pair2 = temp_list_fastqtypes.count('pair2')
        temp_n_single = temp_list_fastqtypes.count('single')
        #Treat only merge treatment samples
        if 'merge' not in temp_set_treatment:
            pass
        else:
            #Paired files with orphan single fastq files with the same number of files and more that one fast each
            if temp_n_pair1 > 1 and temp_n_pair2 > 1 and temp_n_single > 1 and temp_n_pair1 == temp_n_pair2 and temp_n_pair1 == temp_n_single:
                pass
            #More than one single fastq file but no pair1 and no pair2 files
            elif temp_n_pair1 == 0 and temp_n_pair2 == 0 and temp_n_single > 1:
                pass
            #More than one fastq file for pair1 and pair2 with the same number of files, no single files 
            elif temp_n_single == 0 and temp_n_pair1 > 1 and temp_n_pair2 > 1 and temp_n_pair1 == temp_n_pair2:
                pass
            #Anything else is not an acceptable merge configuration
            else:
                print_strn = '\n- ' + rich_text_colored('Sample Name: ', 'check_color', color_treatment) + sample
                print_strn = print_strn + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + ''.join(['Number of pair1(s) = ', str(temp_n_pair1), '; Number of pair2(s) = ', str(temp_n_pair2), '; Number of single(s) = ', str(temp_n_single)])
                check_frases.append(print_strn)
    
    #If check_frases is not empty raise exception
    if len(check_frases) > 0:
        #Prepare exception frase
        frase0 = 'Error! Some of the samples in the Treatment Template have incompatible file configurations for merge mode!\n'
        frase1 = '\nThe permitted configurations are:'
        frase2 = '\n- ' + rich_text_colored('An equal number of PAIRED and SINGLE Fastq files with more than 1 file per Fastq type', 'acceptable', color_treatment) + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + 'Number of pair1(s) > 1; Number of pair2(s) > 1; Number of single(s) > 1;\n                 Number of pair1(s) = Number of pair2(s) = Number of single(s)'
        frase3 = '\n- ' + rich_text_colored('An equal number of PAIRED Fastq files with more than 1 file per Fastq type', 'acceptable', color_treatment) + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + 'Number of pair1(s) > 1; Number of pair2(s) > 1; Number of single(s) = 0;\n                 Number of pair1(s) = Number of pair2(s)'
        frase4 = '\n- ' + rich_text_colored('More than one SINGLE Fastq file', 'acceptable', color_treatment) + rich_text_colored('\n  Configuration: ', 'general_text', color_treatment) + 'Number of pair1(s) = 0; Number of pair2(s) = 0; Number of single(s) > 1'
        frase5 = '\n\nThe following samples present incompatibilities:'
        frase = [frase0, frase1, frase2, frase3, frase4, frase5]
        final_frase = frase + check_frases
        #Raise exception
        raise OMD_CTK_Exception(' '.join(final_frase))


def check_required_variables_in_metadata_table(metadata_df, required_cols, color_treatment):
    """
    This functions checks if all required variables in dictionary are in the provided metadata table.

    Parameters
    ----------
    metadata_df : pandas dataframe
        Metadata Table to check.
    required_cols : pandas series(column)
        Required columns as indicated in Variables Dictionary.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text
    Raises
    ------
    Exception
        If required columns are missing raises an exception.
        Handled by treat_headers_check function.

    Returns
    -------
    None.

    """   
    #Check that the required col_names are present
    requireness = required_cols.isin(list(metadata_df.columns)).all()
    #Get list of required cols
    required_cols_list = list(required_cols)
    if requireness == False:
        ##Check valid 'requiredness' column values
        frase0 = 'Error! Some of the values of the required columns in the Variables Dictionary are not present in the last Metadata Table!\n Check your Variables Dictionary file and/or Metadata Table file!\n'
        frase1 = '\nRequired columns in Variables Dictionary are:'
        ##The exception will be handled by treat_headers_check function
        treat_headers_check(required_cols_list, requireness, metadata_df, frase0, frase1, color_treatment)
        
        
def check_metadata_table_vars_in_dict(metadata_df, universe_cols, color_treatment):
    """
    This functions checks that all metadata columns are present in the variables dictionary.
    
    Parameters
    ----------
    metadata_df : pandas dataframe
        Metadata Table to check.
    universe_cols : list
        List with all the provided variables in Variables Dictionary.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text
    Raises
    ------
    Exception
        If extra columns are present raises an exception.

    Returns
    -------
    None.

    """   
    #Check that all col_names in metadata table exist in variables_dict
    all_present = set(list(metadata_df.columns)).issubset(universe_cols)
    if all_present == False:
        #Prepare exception frase
        frase0 = 'Error! Some of the columns in the last Metadata Table are not present in the Variables Dictionary!\n Check your Variables Dictionary file and/or Metadata Table file!\n'
        frase1 = '\nColumns for this Metadata Table are:'
        frase = [frase0, frase1]
        for col in list(metadata_df.columns):
            if col in universe_cols:
                temp_line = ' '.join(['\n-', rich_text_colored(col, 'check_color', color_treatment), 
                                      rich_text_colored('Found!', 'true_color', color_treatment)])
            else:
                temp_line = ' '.join(['\n-', rich_text_colored(col, 'check_color', color_treatment), 
                                      rich_text_colored('Not Found!', 'false_color', color_treatment)])
            frase.append(temp_line)
        #Raise exception
        raise OMD_CTK_Exception(' '.join(frase))


def check_duplicates_in_variables_dict_column(variable_dict_df, col_name, n_elements):
    """
    This function checks if there are duplicates in the indicated pandas column of the provided variables dictionary dataframe.

    Parameters
    ----------
    variable_dict_df : pandas dataframe
        The provided variables dictionary dataframe.
    col_name : str
        Name of the column to check for duplicates.
    n_elements : int
        The number of elements to be printed by line.

    Raises
    ------
    OMD_CTK_Exception
        If duplicates are detected raises an exception.

    Returns
    -------
    None.

    """
    #Get duplicates df
    duplicates_df = variable_dict_df[variable_dict_df.duplicated(subset=[col_name])]
    
    #Treat duplicates
    if len(duplicates_df) > 0:
        #Get list of duplicates
        duplicates = list(duplicates_df[col_name])
        #Construct frase
        frase0 = ''.join(['Error! Some of the values of the "', col_name,'" column in the Variables Dictionary are duplicates!\n Check your Variables Dictionary file!\n'])
        frase1 = '\nThe following variables present duplicates:\n'
        frase = [frase0, frase1]
        #Split list of duplicates for aesthetics
        splited_list = [duplicates[i:i+n_elements] for i in range(0, len(duplicates), n_elements)]
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


def curly_2_straight_quotation(value):
    """
    This function replaces curly for straight quotation
    in a string.

    Parameters
    ----------
    value : str
        Provided str to replace characters.

    Returns
    -------
    str
        Result str with replaced characters.

    """
    #curly_2_straight_quotation
    return value.replace("“", '"').replace("”", '"').replace("‘", "'").replace("’", "'")


def generic_columns_intersection_checks(common_col_df1, common_col_df2, df1, df2, df1_name, df2_name, color_treatment):
    """
    This function checks intersections between indicated columns of the 
    provided dataframes and shows the corresponding messages if needed.

    Parameters
    ----------
    common_col_df1 : str
        Common intersect column from dataframe1.
    common_col_df2 : str
        Common intersect column from dataframe2.
    df1 : pandas datframe
        The provided dataframe1.
    df2 : pandas datframe
        The provided dataframe2.
    df1_name : str
        The provided pretty dataframe1 name for messages.
    df2_name : str
        The provided pretty dataframe2 name for messages. 
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool
        True: There are intersection warnings.
        False: There are no intersection warnings.

    """
    #Get unique values
    df1_merge_column_values_unique = set(list(df1[common_col_df1]))
    df2_merge_column_values_unique = set(list(df2[common_col_df2]))
    
    #Show stats
    print(rich_text_colored('\nUnique values:', 'subsection2', color_treatment))
    print('  o', rich_text_colored(''.join([df1_name,' Column [',common_col_df1,'] (total unique values):']), 'general_text', color_treatment), len(df1_merge_column_values_unique))
    print('  o', rich_text_colored(''.join([df2_name,' Column [',common_col_df2,'] (total unique values):']), 'general_text', color_treatment), len(df2_merge_column_values_unique))
    
    #Check if all unique values between merge columns are common
    all_df1_in_df2 = df1_merge_column_values_unique.issubset(df2_merge_column_values_unique)
    all_df2_in_df1 = df2_merge_column_values_unique.issubset(df1_merge_column_values_unique)
    if (all_df1_in_df2 and all_df2_in_df1):
        print(rich_text_colored('\nAll unique values are common between the columns provided!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        ##Get lits of non-common values
        df1_difference = list(df1_merge_column_values_unique.difference(df2_merge_column_values_unique))
        df2_difference = list(df2_merge_column_values_unique.difference(df1_merge_column_values_unique))
        
        ##Show intersection stats
        print(rich_text_colored('\nIntersections:', 'subsection2', color_treatment))
        print('  o', rich_text_colored('Number of unique common values between columns:', 'general_text', color_treatment), len(df1_merge_column_values_unique.intersection(df2_merge_column_values_unique)))
        print('  o', rich_text_colored(''.join(['Number of unique specific values for ',df1_name, ' Column:']), 'general_text', color_treatment), len(df1_difference))
        print('  o', rich_text_colored(''.join(['Number of unique specific values for ',df2_name, ' Column:']), 'general_text', color_treatment), len(df2_difference))
        
        ##Show warning message
        print(rich_text_colored('\nWarning! Some values are not common between columns!', 'program_warning', color_treatment))
        
        #Show values if there are warnings
        ##main_difference
        if len(df1_difference) > 0:
            print('\n-', rich_text_colored(''.join(['Specific values for ',df1_name, ' Column:']), 'general_text', color_treatment))
            print_list_n_byline(df1_difference, 5)
        ##extra_difference
        if len(df2_difference) > 0:
            print('\n-', rich_text_colored(''.join(['Specific values for ',df2_name, ' Column:']), 'general_text', color_treatment))
            print_list_n_byline(df2_difference, 5)
        #Return info for lately treat advise messages
        return True
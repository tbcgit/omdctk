#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: check_metadata_ENA.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'check_metadata_ENA.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, ENA_FASTQ_URLS_COLUMNS, program_header, 
                    check_headers, treat_headers_check,
                    print_list_n_byline, get_urls_from_ENA_column,
                    rich_text_colored, show_advise_legend)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
import pandas as pd

#Program Global Variables
headers_used = ['sample_accession','scientific_name','tax_id','instrument_model',
              'instrument_platform','library_layout','library_strategy',
              'library_source','run_accession','fastq_ftp','fastq_aspera',
              'fastq_galaxy','submitted_ftp','submitted_aspera','submitted_galaxy',
              'read_count', 'first_public', 'last_updated']

#Program Functions
def TrueFalseColor(value):
    """
    This function returns different color categories from COLORS_DIC 
    depending on the provided bool values.

    Parameters
    ----------
    value : bool
        Expected values [True, False].

    Returns
    -------
    result : str
        If value True: It returns the str 'green'.
        Else value should be False: It returns the str 'red'.
        
    """
    if value == True:
        result = 'true_color'
    else:
        result = 'false_color'
    return result


def print_run_stats(metadata_df, color_treatment):
    """
    This function obtains and prints the Run Stats checks.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.
    
    """
    #Show the number of run_accessions
    print(rich_text_colored('1) Number of run_accessions:', 'subsection', color_treatment), len(metadata_df.index))
    
    #Show the number of appearances per taxid+scientific_name
    print(rich_text_colored('\n2) Appearances per scientific_name(tax_id):', 'subsection', color_treatment))
    comb1 = metadata_df['scientific_name'] + '(' + metadata_df['tax_id'].astype(str) + ')'
    print(comb1.value_counts(dropna = False).rename_axis('unique_values').reset_index(name = 'counts').to_string(index = False))
    
    #Show the number of appearances per instrument_platfrom+instrument_model
    print(rich_text_colored('\n3) Appearances per instrument_model(instrument_platform):', 'subsection', color_treatment))
    comb2 = metadata_df['instrument_model'] + '(' + metadata_df['instrument_platform'] + ')'
    print(comb2.value_counts(dropna = False).rename_axis('unique_values').reset_index(name = 'counts').to_string(index = False))
    
    #Show the number of appearances per library_layout
    print(rich_text_colored('\n4) Appearances per library_layout:', 'subsection', color_treatment))
    print(metadata_df['library_layout'].value_counts(dropna = False).rename_axis('unique_values').reset_index(name = 'counts').to_string(index = False))
    
    #Show the number of appearances per library_strategy+library_source
    print(rich_text_colored('\n5) Appearances per library_strategy+library_source:', 'subsection', color_treatment))
    comb3 = metadata_df['library_strategy'] + '+' + metadata_df['library_source']
    print(comb3.value_counts(dropna = False).rename_axis('unique_values').reset_index(name = 'counts').to_string(index = False))


def check_PAIRED_SINGLE(metadata_df, ena_download_column, fastq_file_pattern, color_treatment):
    """
    This function checks library_layouts at the run_accessions level. It also 
    prints the check results and shows the corresponding messages if needed.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    ena_download_column : str
        The provided urls column from the metadata table to do checks on.
    fastq_file_pattern : str
        Pattern to identify fastq files.
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
    ##Init lists of run_accessions with warnings
    runs_warnings_PD_NA = []
    runs_warnings_PD_more2 = []
    runs_warnings_PD_only1 = []
    runs_warnings_PD_noFQ = []
    runs_warnings_SG_NA = []
    runs_warnings_SG_more1 = []
    runs_warnings_SG_noFQ = []
    runs_warnings_noSG_noPD = []
    
    #Iter run_accessions
    for i,r in metadata_df.iterrows():
        temp_run = r['run_accession']
        temp_library_layout = r['library_layout']
        temp_fastq_urls = r[ena_download_column]
        #Append to the different warning lists
        if temp_library_layout == 'PAIRED':
            if pd.isna(temp_fastq_urls) == True:
                runs_warnings_PD_NA.append(temp_run)
            else:
                if temp_fastq_urls.count(fastq_file_pattern) == 2:
                    pass
                elif temp_fastq_urls.count(fastq_file_pattern) > 2:
                    runs_warnings_PD_more2.append(temp_run)
                elif temp_fastq_urls.count(fastq_file_pattern) == 1:
                    runs_warnings_PD_only1.append(temp_run)
                else:
                    runs_warnings_PD_noFQ.append(temp_run)
        elif temp_library_layout == 'SINGLE':
            if pd.isna(temp_fastq_urls) == True:
                runs_warnings_SG_NA.append(temp_run)
            else:
                if temp_fastq_urls.count(fastq_file_pattern) == 1:
                    pass
                elif temp_fastq_urls.count(fastq_file_pattern) > 1:
                    runs_warnings_SG_more1.append(temp_run)
                else:
                    runs_warnings_SG_noFQ.append(temp_run)
        else:
            runs_warnings_noSG_noPD.append(temp_run)
    
    #Calculate total number of warning runs
    n_run_warnings = len(runs_warnings_PD_NA + runs_warnings_PD_more2 + runs_warnings_PD_only1 + runs_warnings_PD_noFQ + runs_warnings_SG_NA + runs_warnings_SG_more1 + runs_warnings_SG_noFQ + runs_warnings_noSG_noPD)
    
    #Show total number of runs again and total warnings
    print('-', rich_text_colored('Total number of runs:', 'general_text', color_treatment), len(metadata_df))
    print('-', rich_text_colored('Total number of warning runs:', 'general_text', color_treatment), n_run_warnings)
    
    #Show different information depending on result 
    if n_run_warnings == 0:
        print(rich_text_colored('\nAll library_layouts for ENA Download Column match for all run_accessions!', 'acceptable', color_treatment))
        #Show SINGLE/PAIRED lines depending on if there are any PAIRED/SINGLE
        if 'PAIRED' in list(metadata_df['library_layout']):
            print('-', rich_text_colored('All PAIRED run_accessions have only 2 associated Fastq files', 'acceptable', color_treatment))
        if 'SINGLE' in list(metadata_df['library_layout']):
            print('-', rich_text_colored('All SINGLE run_accessions have only 1 associated Fastq file', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        #Show SINGLE/PAIRED lines depending on if there are any PAIRED/SINGLE
        print(rich_text_colored('\nPAIRED:', 'subsection2', color_treatment))
        print('  o', rich_text_colored('Number of warning runs with empty values in the provided ENA Download Column:', 'general_text', color_treatment), len(runs_warnings_PD_NA))
        print('  o', rich_text_colored('Number of warning runs with more than 2 Fastqs in the provided ENA Download Column:', 'general_text', color_treatment), len(runs_warnings_PD_more2))
        print('  o', rich_text_colored('Number of warning runs with only 1 Fastq in the provided ENA Download Column:', 'general_text', color_treatment), len(runs_warnings_PD_only1))
        print('  o', rich_text_colored('Number of warning runs without Fastq File Pattern in the provided ENA Download Column:', 'general_text', color_treatment), len(runs_warnings_PD_noFQ))
        print(rich_text_colored('SINGLE:', 'subsection2', color_treatment))
        print('  o', rich_text_colored('Number of warning runs with empty values in the provided ENA Download Column:', 'general_text', color_treatment), len(runs_warnings_SG_NA))
        print('  o', rich_text_colored('Number of warning runs with more than 1 Fastq in the provided ENA Download Column:', 'general_text', color_treatment), len(runs_warnings_SG_more1))
        print('  o', rich_text_colored('Number of warning runs without Fastq File Pattern in the provided ENA Download Column:', 'general_text', color_treatment), len(runs_warnings_SG_noFQ))
        print(rich_text_colored('Unexpected library_layouts:', 'subsection2', color_treatment))
        print('  o', rich_text_colored('Number of warning runs with unexpected library_layouts:', 'general_text', color_treatment), len(runs_warnings_noSG_noPD))
        
        #Print warnings message
        print(rich_text_colored('\nWarning! Some runs present unexpected values!', 'program_warning', color_treatment))
        
        #Show run_accession values if there are warnings
        ##PD_NA
        if len(runs_warnings_PD_NA) > 0:
            print('\n-', rich_text_colored('PAIRED runs with empty values in the provided ENA Download Column:', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_PD_NA, 5)
        ##PD_MORE2
        if len(runs_warnings_PD_more2) > 0:
            print('\n-', rich_text_colored('PAIRED runs with more than 2 Fastqs in the provided ENA Download Column:', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_PD_more2, 5)
        ##PD_ONLY1
        if len(runs_warnings_PD_only1) > 0:
            print('\n-', rich_text_colored('PAIRED runs with only 1 Fastq in the provided ENA Download Column:', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_PD_only1, 5)
        ##PD_notFQ
        if len(runs_warnings_PD_noFQ) > 0:
            print('\n-', rich_text_colored('PAIRED runs without Fastq File Pattern in the provided ENA Download Column:', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_PD_noFQ, 5)            
        ##SG_NA
        if len(runs_warnings_SG_NA) > 0:
            print('\n-', rich_text_colored('SINGLE runs with empty values in the provided ENA Download Column:', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_SG_NA, 5)
        ##SG_MORE1
        if len(runs_warnings_SG_more1) > 0:
            print('\n-', rich_text_colored('SINGLE runs with more than 1 Fastq in the provided ENA Download Column:', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_SG_more1, 5)
        ##SG_notFQ
        if len(runs_warnings_SG_noFQ) > 0:
            print('\n-', rich_text_colored('SINGLE runs without Fastq File Pattern in the provided ENA Download Column:', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_SG_noFQ, 5)
        ##notPD_notSG
        if len(runs_warnings_noSG_noPD) > 0:
            print('\n-', rich_text_colored('Runs with unexpected layouts (other than PAIRED, or SINGLE):', 'general_text', color_treatment))
            print_list_n_byline(runs_warnings_noSG_noPD, 5)
        
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:', 'due_to_header', color_treatment))
        print('- The presence of quality controlled Fastqs in the dataset', rich_text_colored('[Acceptable]', 'acceptable', color_treatment))
        print('- The usage of a Fastq file pattern different from the default (".fastq.gz")', rich_text_colored('[Acceptable]', 'acceptable', color_treatment))
        print('- The use of an ENA Download Column that is not available', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print('- The utilization of a Private ENA Project', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm which is your case')
        print('- If any, check the presence of extra Fastqs and the number of reads associated')
        print('- If any, check missing Fastq cases')
        print('- Check the original Fastq names and content (if available)')
        print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')
        print('- Check the original publication and supplementary tables to get some context')
        print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #Return info for lately treat advise messages
        return True


def check_ENA_submitted_files_columns(metadata_df, color_treatment):
    """
    This function checks if the ENA submitted files columns are empty. It also 
    prints the check results and shows the corresponding messages if needed.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
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
    #Check emptyness and show messages
    empty_submitted_ftp = metadata_df['submitted_ftp'].dropna().empty
    print('-', rich_text_colored('Is submitted_ftp available?', TrueFalseColor(empty_submitted_ftp == False), color_treatment), empty_submitted_ftp == False)
    empty_submitted_aspera = metadata_df['submitted_aspera'].dropna().empty
    print('-', rich_text_colored('Is submitted_aspera available?', TrueFalseColor(empty_submitted_aspera == False), color_treatment), empty_submitted_aspera == False)
    empty_submitted_galaxy = metadata_df['submitted_galaxy'].dropna().empty
    print('-', rich_text_colored('Is submitted_galaxy available?', TrueFalseColor(empty_submitted_galaxy == False), color_treatment), empty_submitted_galaxy == False)
    
    if (empty_submitted_ftp and empty_submitted_aspera and empty_submitted_galaxy):  
        print(rich_text_colored('\nNone of the original submitted Fastqs columns are available!', 'program_warning', color_treatment))
        
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:', 'due_to_header', color_treatment))
        print('- That the dataset was not originally uploaded to ENA, but rather to another database', rich_text_colored('[Acceptable]', 'acceptable', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm that this is your case')
        print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')

        #Return info for lately treat advise messages
        return True
    
    else:
        print(rich_text_colored('\nThere are original submitted Fastqs columns available!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    

def match_ENA_submitted_files_columns(column):
    """
    This function matches the provided urls ENA_metadata_column 
    (fastq_ftp, fastq_aspera, fastq_galaxy) with the 
    corresponding submitted urls column.

    Parameters
    ----------
    column : str
        The ENA metadata column for the ENA fastq files.

    Returns
    -------
    submited_files_col_name : str
        The name of the corresponding submited original files.
        
    """
    if column == 'fastq_galaxy':
        submited_files_col_name = 'submitted_galaxy'
    elif column == 'fastq_aspera':
        submited_files_col_name = 'submitted_aspera'
    else:
        submited_files_col_name = 'submitted_ftp'
    return submited_files_col_name


def select_ENA_submitted_files_column(ena_download_column, metadata_df, color_treatment):
    """
    This function selects an ENA submitted files column if possible.

    Parameters
    ----------
    ena_download_column : str
        The ENA download column provided.
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    selected_submitted_col_name : None or str
        If all submitted files columns are empty, returns None.
        Else it returns a not empty selected ENA submitted files column.

    """
    #Set submitted_cols list
    submitted_cols = ['submitted_ftp', 'submitted_aspera', 'submitted_galaxy']
    
    #Get empty checks
    empty_submitted_ftp = metadata_df['submitted_ftp'].dropna().empty
    empty_submitted_aspera = metadata_df['submitted_aspera'].dropna().empty
    empty_submitted_galaxy = metadata_df['submitted_galaxy'].dropna().empty
    
    #Treat different outcomes
    if (empty_submitted_ftp and empty_submitted_aspera and empty_submitted_galaxy):
        print(rich_text_colored('\nNone of the original submitted Fastqs columns are available!', 'program_warning', color_treatment))
        print(rich_text_colored('Duplicated file names can not be searched!', 'program_warning', color_treatment))
        selected_submitted_col_name = None
    
    elif ena_download_column in submitted_cols and metadata_df[ena_download_column].dropna().empty == False:
        print(rich_text_colored('Using the following submitted Fastqs column:', 'column_color', color_treatment), ena_download_column)
        selected_submitted_col_name = ena_download_column
    
    elif ena_download_column in ['fastq_ftp', 'fastq_aspera', 'fastq_galaxy'] and metadata_df[match_ENA_submitted_files_columns(ena_download_column)].dropna().empty == False:
        submitted_match = match_ENA_submitted_files_columns(ena_download_column)
        print(rich_text_colored('Using the following submitted Fastqs column:', 'column_color', color_treatment), submitted_match)
        selected_submitted_col_name = submitted_match
    
    else:
        #For the rest we should at least have one submitted column not empty
        ##Find the first submitted not empty column
        submitted_col_values = [empty_submitted_ftp, empty_submitted_aspera, empty_submitted_galaxy]
        for i in range(len(submitted_col_values)):
            if submitted_col_values[i] == False:
                submitted_col_name = submitted_cols[i]
                break
        ##Indicate submitted column used
        print(rich_text_colored('Using the following submitted Fastqs column:', 'column_color', color_treatment), submitted_col_name)
        ##Set selected_submitted_col_name
        selected_submitted_col_name = submitted_col_name
    
    return selected_submitted_col_name
        
        
def find_duplicated_names_submitted_column(metadata_df, submitted_column, sample_column, color_treatment):
    """
    This function finds duplicated file names for the provided FASTQs 
    submitted column and prints the corresponding messages.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    submitted_column : str
        The ENA metadata submitted column from which to obtain the file names.
     sample_column : str
        The sample_column provided in the -s parameter.
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
    #Get list of urls from the submitted FASTQ files column provided
    urls = get_urls_from_ENA_column(metadata_df, submitted_column)   
    
    #Get file names list
    file_names = [i.split('/')[-1] for i in urls]
    
    #Find duplicates
    seen = set()
    duplicates = list(set([i for i in file_names if i in seen or seen.add(i)]))
    
    #Treat different outcomes
    if len(duplicates) > 0:
        print(rich_text_colored('\nTotal of warning duplicated files:', 'program_warning', color_treatment), len(duplicates))
        
        print('The following duplicates were detected:')
        for dup in duplicates:
            temp_dup_info = metadata_df[metadata_df[submitted_column].str.contains(dup,regex=False)][['sample_accession', sample_column, 'run_accession', 'library_layout', 'read_count', 'first_public', 'last_updated']]
            print(rich_text_colored('\nFile name:', 'program_warning', color_treatment), dup)
            print(rich_text_colored("File's information:", 'general_text', color_treatment))
            print(temp_dup_info.to_string(index = False))
            print(rich_text_colored('Number of run_accessions:', 'general_text', color_treatment), len(temp_dup_info))
            print(rich_text_colored('Number of unique sample_accessions:', 'general_text', color_treatment), temp_dup_info['sample_accession'].nunique(dropna = False))
            message_scol = ''.join(['Number of unique samples for "', sample_column, '":'])
            print(rich_text_colored(message_scol, 'general_text', color_treatment), temp_dup_info[sample_column].nunique(dropna = False))
        
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:', 'due_to_header', color_treatment))
        print('- The utilization of a Private ENA Project', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm which is your case')
        print('- Check these files and confirm that they are indeed the same files')
        print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        print('- Select a representative file for each duplicate')
        
        #Return info for lately treat advise messages
        return True
    
    else:
        print(rich_text_colored('No duplicates were detected!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False


def print_sample_stats(metadata_df, orderID, col_name, color_treatment):
    """
    This function obtains and prints the Sample Stats for the provided column.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    orderID : int
        The order id for the stats analysis messages.
    col_name : str
        The name of the sample_column to get the sample stats.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.
    
    """
    ##1.Show sample_column stats message
    print(rich_text_colored(''.join([str(orderID), ') Stats for "', col_name, '":']),'subsection', color_treatment))
    
    ##1.1.Show Number of unique samples for col_name
    message11 = ''.join(['Number of unique samples for "', col_name, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.1)']), 'subsection2', color_treatment), rich_text_colored(message11, 'subsection3', color_treatment), metadata_df[col_name].nunique(dropna = False))
    
    ##1.2. Show number of samples per scientific_name(tax_id) combination
    message12 = ''.join(['Samples per scientific_name(tax_id) in "', col_name, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.2)']), 'subsection2', color_treatment), rich_text_colored(message12, 'subsection3', color_treatment))
    #Get values by sample
    precomb1 = metadata_df[['scientific_name', 'tax_id', col_name]].value_counts(dropna = False).reset_index(name = 'count')
    comb1 = precomb1['scientific_name'] + '(' + precomb1['tax_id'].astype(str) + ')'
    print(comb1.value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts').to_string(index = False))
    
    ##1.3. Show number of samples per instrument_model(instrument_platform) combination
    message13 = ''.join(['Samples per instrument_model(instrument_platform) in "', col_name, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.3)']), 'subsection2', color_treatment), rich_text_colored(message13,'subsection3', color_treatment))
    #Get values by sample
    precomb2 = metadata_df[['instrument_model', 'instrument_platform', col_name]].value_counts(dropna = False).reset_index(name = 'count')
    comb2 = precomb2['instrument_model'] + '(' + precomb2['instrument_platform'] + ')'
    print(comb2.value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts').to_string(index = False))
    
    ##1.4. Show number of samples per library_layout
    message14 = ''.join(['Samples per library_layout in "', col_name, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.4)']), 'subsection2', color_treatment), rich_text_colored(message14, 'subsection3', color_treatment))
    #Get values by sample
    comb3 = metadata_df[['library_layout', col_name]].value_counts(dropna = False).reset_index(name = 'count')
    print(comb3['library_layout'].value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts').to_string(index = False))
    
    ##1.5. Show number of samples per library_strategy+library_source combination
    message15 = ''.join(['Samples per library_strategy+library_source in "', col_name, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.5)']), 'subsection2', color_treatment), rich_text_colored(message15, 'subsection3', color_treatment))
    #Get values by sample
    precomb4 = metadata_df[['library_strategy', 'library_source', col_name]].value_counts(dropna = False).reset_index(name = 'count')
    comb4 = precomb4['library_strategy'] + '+' + precomb4['library_source']
    print(comb4.value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts').to_string(index = False))
    
    ##1.6.Show groups of samples of col_name by number of run_accessions
    message16 = ''.join(['Groups of samples in "', col_name, '" by number of run_accessions:'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.6)']), 'subsection2', color_treatment), rich_text_colored(message16, 'subsection3', color_treatment))
    ##Generate table of number of run_accessions per samples in col_name
    runs_by_sample = metadata_df[col_name].value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts')
    ##Get groups/ run_accessions per sample_accession
    uniq_count_groups = sorted(list(set(list(runs_by_sample['counts']))))
    ##Set messages
    message12A = rich_text_colored('\nGroup with ', 'general_text', color_treatment)
    message12B = rich_text_colored(' run_accession(s) per sample:', 'general_text', color_treatment)
    message12C = rich_text_colored('Total number of samples in this group:', 'general_text', color_treatment)
    ##Differenciate by number of groups
    if len(uniq_count_groups) == 1:
        ##If there is only one group do not show samples IDs
        print(message12A, uniq_count_groups[0], message12B)
        print('-', message12C, len(runs_by_sample))
    else:    
        ##Show chunks of 5 samples per line if there is more than 1 group
        for group in uniq_count_groups:
            print(message12A, group, message12B)
            #Get temp sample_accessions list
            temp_samples_list = sorted(list(runs_by_sample[runs_by_sample['counts'] == group]['values']))
            print('-', message12C, len(temp_samples_list))
            print('-', rich_text_colored("Showing samples' values:", 'general_text', color_treatment))
            #Split list of run_accessions in chuncks of 5 and show chunks
            print_list_n_byline(temp_samples_list, 5)
    
    
def sample_vs_sample(metadata_df, orderID, suborderID, main_sample_column, other_sample_column, color_treatment):
    """
    This function compares the provided sample_columns in the metadata table.
    It also prints the check results and returns bool for treat advise messages.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    orderID : int
        The order id for the stats analysis messages.
    suborderID : int
        The sub-order id for the stats analysis messages.
    main_sample_column : str
        The name of the sample_column to use as main sample.
    other_sample_column : str
        The name of the sample_column to compare with the main_sample_column.
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
    #Show sample_column check message
    message = ''.join(['Check samples for "', main_sample_column, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.', str(suborderID), ')']), 'subsection2', color_treatment), rich_text_colored(message, 'subsection3', color_treatment))
    
    #Get main_sample_column vs other_sample_column combination
    comb = metadata_df[[main_sample_column, other_sample_column]].value_counts(dropna = False).reset_index(name = 'count')
    
    #Number of samples in other_sample_column per sample in main_sample_column
    df_counts = comb[main_sample_column].value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts')
    
    #Get number of samples with warnings
    n_sample_warnings = len(df_counts[df_counts['counts'] > 1])
    
    #Show total number of samples in main_sample_column and total warnings
    print('-', rich_text_colored('Total number of samples:', 'general_text', color_treatment), metadata_df[main_sample_column].nunique(dropna = False))
    print('-', rich_text_colored('Total number of warning samples:', 'general_text', color_treatment), n_sample_warnings)
    
    #Show different information depending on result 
    if n_sample_warnings == 0:
        message_okey = ''.join(['All samples in "', main_sample_column, '" have only one match with "', other_sample_column, '"!'])
        print(rich_text_colored(message_okey, 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        print(rich_text_colored(''.join(['\nWarning! Some samples in "', main_sample_column, '" have more than one match with "', other_sample_column, '"!']), 'program_warning', color_treatment))
        print('\n-', rich_text_colored("Showing samples' values:", 'general_text', color_treatment))
        #Show sample ids for main_sample_column in chunks of 5
        sample_warnings_list = sorted(list(df_counts[df_counts['counts'] > 1]['values']))
        print_list_n_byline(sample_warnings_list, 5)
        #Return info for lately treat advise messages
        return True
    
    
def sample_vs_column(metadata_df, orderID, suborderID, main_sample_column, columns_list, column_comb_name, color_treatment):
    """
    This function checks the appearances of the provided column in the metadata table 
    with the samples in the provided main_sample_column. It also prints the check 
    results and returns bool for treat the advise messages lately.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    orderID : int
        The order id for the stats analysis messages.
    suborderID : int
        The sub-order id for the stats analysis messages.
    main_sample_column : str
        The name of the sample_column to use as main sample.
    columns_list : list
        The name of the list with the columns of interest to calculate main_sample_column appearances per sample.
    column_comb_name : str
        Name of the column or columns combination str provided, to use in messages.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool
        True: There are sample warnings.
        False: There are no sample warnings.

    """
    #Set updated list to filter dataframe
    list_updated = [main_sample_column] + columns_list
    
    #Show sample_column check message
    message = ''.join(['Check samples for "', main_sample_column, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.', str(suborderID), ')']), 'subsection2', color_treatment), rich_text_colored(message, 'subsection3', color_treatment))
    
    #Get main_sample_column vs other_sample_column combination
    comb = metadata_df[list_updated].value_counts(dropna = False).reset_index(name = 'count')
    
    #Number of samples in other_column_comb per sample in main_sample_column
    df_counts = comb[main_sample_column].value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts')
    
    #Get number of samples with warnings
    n_sample_warnings = len(df_counts[df_counts['counts'] > 1])
    
    #Show total number of samples in main_sample_column and total warnings
    print('-', rich_text_colored('Total number of samples:', 'general_text', color_treatment), metadata_df[main_sample_column].nunique(dropna = False))
    print('-', rich_text_colored('Total number of warning samples:', 'general_text', color_treatment), n_sample_warnings)
    
    #Show different information depending on result 
    if n_sample_warnings == 0:
        message_okey = ''.join(['All samples in "', main_sample_column, '" have only one match per "', column_comb_name, '"!'])
        print(rich_text_colored(message_okey, 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        print(rich_text_colored(''.join(['\nWarning! Some samples in "', main_sample_column, '" present mixed "', column_comb_name, 's"!']), 'program_warning', color_treatment))
        print('\n-', rich_text_colored("Showing samples' values:", 'general_text', color_treatment))
        #Show sample ids for main_sample_column in chunks of 5
        sample_warnings_list = sorted(list(df_counts[df_counts['counts'] > 1]['values']))
        print_list_n_byline(sample_warnings_list, 5)
        #Return info for lately treat advise messages
        return True
    

def sample_vs_library_layout(metadata_df, orderID, suborderID, main_sample_column, color_treatment):
    """
    Modification of sample_vs_column to get extra checks in the case of library_layout.
    This function checks library_layout appearances in the metadata table with 
    with the samples in the provided main_sample_column. It also prints the check 
    results and returns bool for treat the advise messages lately.

    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    orderID : int
        The order id for the stats analysis messages.
    suborderID : int
        The sub-order id for the stats analysis messages.
    main_sample_column : str
        The name of the sample_column to use as main sample.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    bool
        True: There are sample warnings.
        False: There are no sample warnings.

    """
    #Show sample_column check message
    message = ''.join(['Check samples for "', main_sample_column, '":'])
    print(rich_text_colored(''.join(['\n', str(orderID), '.', str(suborderID), ')']), 'subsection2', color_treatment), rich_text_colored(message, 'subsection3', color_treatment))
    
    #Calculate in which samples mean read count per sample+library_layout
    pp = metadata_df[[main_sample_column, 'library_layout', 'read_count']]
    read_means = pp.groupby([main_sample_column, 'library_layout']).mean().rename(columns = {'read_count':'reads_mean'})
    
    #Get list of samples
    sampleIDs_list = sorted(list(set(metadata_df[main_sample_column].tolist())))
    
    #Set init lists
    sample_avg_reads_check_okey_greter = []
    sample_avg_reads_check_okey_equal = []
    sample_avg_reads_check_okey_less = []
    sample_avg_reads_check_rare = []
    
    #Check each sample
    for sampleID in sampleIDs_list:
        ##Get library_layout values for sample
        temp_sample_layouts = set(list(metadata_df.loc[metadata_df[main_sample_column] == sampleID]['library_layout']))
        ##Check if they have both PAIRED AND SINGLE and there are only 2 types/ Then check avg read counts
        if len(temp_sample_layouts) == 2 and ('PAIRED' in temp_sample_layouts) and ('SINGLE' in temp_sample_layouts):
            if read_means.loc[(sampleID, 'PAIRED')]['reads_mean'] > read_means.loc[(sampleID, 'SINGLE')]['reads_mean']:
                sample_avg_reads_check_okey_greter.append(sampleID)
            elif read_means.loc[(sampleID, 'PAIRED')]['reads_mean'] < read_means.loc[(sampleID, 'SINGLE')]['reads_mean']:
                sample_avg_reads_check_okey_less.append(sampleID)
            else:
                sample_avg_reads_check_okey_equal(sampleID)
        #If there is only one layout type and is PAIRED or SINGLE do nothing
        elif len(temp_sample_layouts) == 1 and (('PAIRED' in temp_sample_layouts) or ('SINGLE' in temp_sample_layouts)):
            pass
        ##Put here the weird cases/ Anything else
        else:
            sample_avg_reads_check_rare.append(sampleID)
    
    #Calculate the number of sample warnings
    n_sample_warnings = len(sample_avg_reads_check_okey_greter + sample_avg_reads_check_okey_equal + sample_avg_reads_check_okey_less + sample_avg_reads_check_rare)
    
    #Show total number of samples in main_sample_column and total warnings
    print('-', rich_text_colored('Total number of samples:', 'general_text', color_treatment), metadata_df[main_sample_column].nunique(dropna = False))
    print('-', rich_text_colored('Total number of warning samples:', 'general_text', color_treatment), n_sample_warnings)
    
    #Show different information depending on result 
    if n_sample_warnings == 0:
        message_okey = ''.join(['All samples in "', main_sample_column, '" have only one type of library_layout!'])
        print(rich_text_colored(message_okey, 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        #Print extra stats
        print('  o', rich_text_colored('Number of warning samples with mean read_counts in PAIRED > SINGLE:', 'general_text', color_treatment), len(sample_avg_reads_check_okey_greter))
        print('  o', rich_text_colored('Number of warning samples with mean read_counts in PAIRED < SINGLE:', 'general_text', color_treatment), len(sample_avg_reads_check_okey_less))
        print('  o', rich_text_colored('Number of warning samples with mean read_counts in PAIRED = SINGLE:', 'general_text', color_treatment), len(sample_avg_reads_check_okey_equal))
        print('  o', rich_text_colored('Number of warning samples with unexpected layouts:', 'general_text', color_treatment), len(sample_avg_reads_check_rare))
        
        #Print warning message
        print(rich_text_colored(''.join(['\nWarning! Some samples in "',main_sample_column,'" present mixed "library_layouts"!']),'program_warning', color_treatment))
        
        #If warnings sample subtypes are not empty print each one
        ##Greater PAIRED
        if len(sample_avg_reads_check_okey_greter) > 0:
            print('\n-', rich_text_colored("Samples' values with mean read_counts in PAIRED > SINGLE:", 'general_text', color_treatment))
            print_list_n_byline(sample_avg_reads_check_okey_greter, 5)
        ##Less PAIRED
        if len(sample_avg_reads_check_okey_less) > 0:
            print('\n-', rich_text_colored("Samples' values with mean read_counts in PAIRED < SINGLE:", 'general_text', color_treatment))
            print_list_n_byline(sample_avg_reads_check_okey_less, 5)
        ##Equal PAIRED=SINGLE
        if len(sample_avg_reads_check_okey_equal) > 0:
            print('\n-', rich_text_colored("Samples' values with mean read_counts in PAIRED = SINGLE:", 'general_text', color_treatment))
            print_list_n_byline(sample_avg_reads_check_okey_equal, 5)
        ##Weird cases
        if len(sample_avg_reads_check_rare) > 0:
            print('\n-', rich_text_colored("Samples' values with unexpected layouts (other than PAIRED, or SINGLE):", 'general_text', color_treatment))
            print_list_n_byline(sample_avg_reads_check_rare, 5)
   
        #Return info for lately treat advise messages
        return True
    
    
def extra_stats_column(metadata_df, sample_column, extra_column, color_treatment):
    """
    This function gets and prints extra stats (run and samples) 
    for the provided extra_column.
    
    Parameters
    ----------
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    sample_column : str
        The sample_column provided in the -s parameter.
    extra_column : str
        Name of the extra_column provided to get stats from.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    None.

    """
    #Show extra_column header message   
    message = ''.join(['Extra stats for "', extra_column, '":'])
    print('\n-', rich_text_colored(message, 'subsection3', color_treatment))
    
    #Set black_list of columns that either have been deeply explored or that does not make sense to get stats
    black_list = ['run_accession', 'sample_accession', sample_column, 'library_layout', 'fastq_ftp', 'fastq_aspera', 'fastq_galaxy', 
                  'submitted_ftp', 'submitted_aspera', 'submitted_galaxy', 'read_count', 'base_count', 'fastq_bytes', 'fastq_md5',
                  'submitted_bytes', 'submitted_md5', 'sra_bytes', 'sra_md5', 'sra_ftp', 'sra_aspera', 'sra_galaxy', 'cram_index_ftp',
                  'cram_index_aspera', 'cram_index_galaxy', 'ENA-SPOT-COUNT', 'ENA-BASE-COUNT']
    
    #Check if extra__column is in black_list
    ## If is in black list show message for skip stats
    if extra_column in black_list:
        print(rich_text_colored('\nSkipping extra column!', 'program_warning2', color_treatment))
        print(rich_text_colored('This extra column has already been deeply explored or does not make sense to get stats!', 'program_warning2', color_treatment))
    ## Else get run and sample stats for extra_column
    else:
        #Get run stats for extra_column and print messages
        run_message = ''.join(['\n1) Run_accessions per "', extra_column, '":'])
        print(rich_text_colored(run_message, 'subsection', color_treatment))
        print(metadata_df[extra_column].value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts').to_string(index = False))
        
        #Get Sample stats for sample_accession column and print 
        sacc_message = ''.join(['\n2) Samples per "', extra_column, '" in "sample_accession":'])
        print(rich_text_colored(sacc_message, 'subsection', color_treatment))
        #Get values by sample and print
        combsacc = metadata_df[[extra_column, 'sample_accession']].value_counts(dropna = False).reset_index(name = 'count')
        print(combsacc[extra_column].value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts').to_string(index = False))
        
        #Get Sample stats for provided sample_column and print
        scol_message = ''.join(['\n3) Samples per "', extra_column, '" in "', sample_column, '":'])
        print(rich_text_colored(scol_message, 'subsection', color_treatment))
        #Get values by sample and print
        combscol = metadata_df[[extra_column, sample_column]].value_counts(dropna = False).reset_index(name = 'count')
        print(combscol[extra_column].value_counts(dropna = False).rename_axis('values').reset_index(name = 'counts').to_string(index = False))


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
    ##Parameter ena_download_column
    parser.add_argument(
            '-c','--ena_download_column', 
            action = 'store',
            choices = ENA_FASTQ_URLS_COLUMNS,
            required = False,
            default = 'fastq_ftp',
            help = 'ENA Download Column (Optional) [Default:fastq_ftp]. Indicate the ENA Metadata Table column with the download links.'
    )
    ##Parameter fastq_pattern
    parser.add_argument(
            '-p','--fastq_pattern', 
            action = 'store',
            default = '.fastq.gz',
            required = False,
            help = 'Fastq File Pattern (Optional) [Default:".fastq.gz"]. Indicate the pattern to identify Fastq files.'
    )
    ##Parameter sample_column
    parser.add_argument(
            '-s','--sample_column', 
            action = 'store',
            default = 'sample_alias',
            required = False,
            help = 'Sample Column (Optional) [Default:sample_alias]. Indicate the Metadata Table column to be used as samples.'
    )
    ##Parameter extra_columns_stats
    parser.add_argument(
            '-e','--extra_columns_stats', 
            nargs = '+',
            required = False,
            help = 'Extra Columns Stats (Optional). Indicate the names for the extra columns to check appearances separated by spaces (If a column name has spaces, quote it).'
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
            action='version',
            version='version {} ({})'.format(VERSION, DATE)
    )
    
    #Process arguments
    args = parser.parse_args()
    metadata_table_path = args.metadata_table
    ena_download_column = args.ena_download_column
    fastq_pattern = args.fastq_pattern
    sample_column = args.sample_column
    extra_columns_names = args.extra_columns_stats
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
        #1)Load files
        
        #Section header message
        print(rich_text_colored('\nLoading File:\n', 'section_header', plain_text_bool))
        
        #Try to load Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('Metadata Table file:', 'general_text', plain_text_bool))
        print(metadata_table_path)
        ##Load metadata file as pandas df
        metadata = pd.read_csv(metadata_table_path, sep = '\t')
        
        #2)Prepare and check headers to be used
        #Set global variable
        global headers_used
        #Prepare headers to use
        ##Add sample_column provided if is not already present
        if sample_column not in headers_used:
            headers_used.append(sample_column)
        ##If extra_column parameters are given add to headers to be used if they are not already present
        if type(extra_columns_names) == list:
            for column in extra_columns_names:
                if column not in headers_used:
                    headers_used.append(column)
        
        #Check headers and treatment
        check_headers_metadata = check_headers(headers_used, metadata)
        frase0 = 'Error! Some of the needed headers are not in the Metadata Table!\n Check your metadata file, --sample_column and --extra_columns_stats parameters!\n'
        frase1 = '\nThe headers needed are:'
        treat_headers_check(headers_used, check_headers_metadata, metadata, frase0, frase1, plain_text_bool)
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
    
    else:
        #3) Runs' Stats and Checks
        
        #3.1)Runs' Stats
        
        #Section header message
        print(rich_text_colored("\nRuns' Stats:\n", 'section_header', plain_text_bool))
        
        #Print main stats
        print_run_stats(metadata, plain_text_bool)
        
        #3.2)Runs' Checks
        
        #Section header message
        print(rich_text_colored("\nRuns' Checks:\n", 'section_header', plain_text_bool))
        
        #3.2.1)For each run_accession check if library_layout and ENA_download_column match
        ##Print info
        print(rich_text_colored('1) Check that library_layout and ENA Download Column match:', 'subsection', plain_text_bool))
        print(rich_text_colored('ENA Download Column selected:', 'column_color', plain_text_bool), ena_download_column)
        print(rich_text_colored('Fastq File Pattern provided:', 'pattern_color', plain_text_bool),fastq_pattern)
        ##Get run layout checks and print results
        runs_layout_checks = check_PAIRED_SINGLE(metadata, ena_download_column, fastq_pattern, plain_text_bool)
        
        #3.2.2)Check if original fastqs/names are available(by checking if is empty)
        print(rich_text_colored('\n2) Check if original uploaded Fastqs are available:', 'subsection', plain_text_bool))
        submitted_files_check = check_ENA_submitted_files_columns(metadata, plain_text_bool)
        
        #3.2.3)Check the original files names, searching for duplicated files
        ##Print info
        print(rich_text_colored('\n3) Check if there are duplicated Fastq file names in the original uploaded Fastqs column:', 'subsection', plain_text_bool))
        ##Select ENA submitted files column
        selected_submitted_col_name = select_ENA_submitted_files_column(ena_download_column, metadata, plain_text_bool)
        ##Search duplicated files if possible / Else Do not raise leyend for this.
        if selected_submitted_col_name != None:
            dup_names_submitted_files_check = find_duplicated_names_submitted_column(metadata, selected_submitted_col_name, sample_column, plain_text_bool)            
        else:
            dup_names_submitted_files_check = False
        #4) Samples' Stats and Checks
        
        #4.1)Samples' Stats
        
        #Section header message
        print(rich_text_colored('\nSample Stats:\n','section_header', plain_text_bool))
        
        #Get Sample stats for "sample_accession" and sample_column
        print_sample_stats(metadata, 1, 'sample_accession', plain_text_bool)
        print('')
        print_sample_stats(metadata, 2, sample_column, plain_text_bool)
        
        #4.2)Samples' Checks
        
        #Section header message
        print(rich_text_colored('\nSample Checks:\n', 'section_header', plain_text_bool))
        
        #4.2.1) Check if number run_accessions match number of samples in sample_columns
        ##Check runs equals samples
        print(rich_text_colored('1) Does the number of run_accessions equal the number of samples?', 'subsection', plain_text_bool))
        ##Check if run_accessions equals sample_accessions
        run_sample_accession = len(metadata.index) == metadata['sample_accession'].nunique(dropna = False)
        print('-', rich_text_colored('Number of runs equals samples in "sample_accessions":', TrueFalseColor(run_sample_accession), plain_text_bool), run_sample_accession)
        ##Check if run_accessions equals sample_alias
        run_sample_column = len(metadata.index) == metadata[sample_column].nunique(dropna = False)
        print('-', rich_text_colored(''.join(['Number of runs equals samples in "', sample_column, '":']), TrueFalseColor(run_sample_column), plain_text_bool), run_sample_column)
        ##Show advise messages if necessary
        if (run_sample_accession == False or run_sample_column == False):
            print(rich_text_colored('\nWarning! The number of runs does not match the number of samples!', 'program_warning', plain_text_bool))
            print(rich_text_colored('\nThis could be due to:','due_to_header', plain_text_bool))
            print('- The presence of different data types in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of quality controlled Fastqs in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of different sequencing technologies', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- Multiple sequencer output Fastqs (runs/lanes) from the same sample', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- Technical replicates from the same sample', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- PAIRED files uploaded as SINGLE Fastq files', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', plain_text_bool))
            print(rich_text_colored('\nYou should:', 'you_should_header', plain_text_bool))
            print('- Manually confirm which is your case (the following sample checks may help)')
            print('- Check the extra Fastqs and the number of reads associated')
            print('- Check the original Fastq names and content (if available)')
            print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')
            print('- Check the original publication and supplementary tables to get some context')
            print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #4.2.2) Check samples match between sample_accession and sample_column provided
        messageB22 = ''.join(['\n2) Does the number of samples equal each other between "sample_accession" and "', sample_column, '"?'])
        print(rich_text_colored(messageB22, 'subsection', plain_text_bool))
        ##Check if number of samples in sample_accession equals provided sample_column
        scol_equals_sacc = metadata['sample_accession'].nunique(dropna = False) == metadata[sample_column].nunique(dropna = False)
        messageB22A = ''.join(['Number of samples in "', sample_column, '" equals samples in "sample_accession":'])
        print('-', rich_text_colored(messageB22A, TrueFalseColor(scol_equals_sacc), plain_text_bool), scol_equals_sacc)
        ##Compare sample_accession vs sample_column
        warnings_sacc = sample_vs_sample(metadata, 2, 1, 'sample_accession', sample_column, plain_text_bool)
        ##Compare sample_column vs sample_accession
        warnings_scol = sample_vs_sample(metadata, 2, 2, sample_column, 'sample_accession', plain_text_bool)
        #Show advise message is warnings are detected
        if scol_equals_sacc == False or warnings_sacc == True or warnings_scol == True:
            #Show advise message
            print(rich_text_colored('\nThis could be due to:', 'due_to_header', plain_text_bool))
            print('- The presence of different data types in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of quality controlled Fastqs in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of different sequencing technologies', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- Multiple sequencer output Fastqs (runs/lines) from the same sample', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- Technical replicates from the same sample', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- PAIRED files uploaded as SINGLE Fastq files', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', plain_text_bool))
            print(rich_text_colored('\nYou should:', 'you_should_header', plain_text_bool))
            print('- Manually confirm which is your case')
            print('- Check the sample IDs of the samples with multiple matches')
            print('- Check the original Fastq names and content (if available)')
            print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')
            print('- Check the original publication and supplementary tables to get some context')
            print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #4.2.3) Check if there is more than one library_strategy+library_source per sample
        messageB23 = '\n3) Is there more than one library_strategy+library_source per sample?'
        print(rich_text_colored(messageB23, 'subsection', plain_text_bool))
        ##Check sample_accession column
        warnings_sacc_B23 = sample_vs_column(metadata, 3, 1,'sample_accession', ['library_strategy', 'library_source'], 'library_strategy+library_source', plain_text_bool)
        ##Check provided sample_column
        warnings_scol_B23 = sample_vs_column(metadata, 3, 2, sample_column, ['library_strategy', 'library_source'], 'library_strategy+library_source', plain_text_bool)
        #Show advise message is warnings are detected
        if warnings_sacc_B23 == True or warnings_scol_B23 == True:
            #Show advise message
            print(rich_text_colored('\nThis could be due to:', 'due_to_header', plain_text_bool))
            print('- The presence of different data types in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print(rich_text_colored('\nYou should:', 'you_should_header', plain_text_bool))
            print('- Manually confirm that this is your case')
            print('- Check the samples with mixed library_strategy+library_source')
            print('- Check the original Fastq names and content (if available)')
            print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')
            print('- Check the original publication and supplementary tables to get some context')
            print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #4.2.4) Check if there is more than one scientific_name(tax_id) per sample
        messageB24 = '\n4) Is there more than one scientific_name(tax_id) per sample?'
        print(rich_text_colored(messageB24, 'subsection', plain_text_bool))
        ##Check sample_accession column
        warnings_sacc_B24 = sample_vs_column(metadata, 4, 1, 'sample_accession', ['scientific_name', 'tax_id'], 'scientific_name(tax_id)', plain_text_bool)
        ##Check provided sample_column
        warnings_scol_B24 = sample_vs_column(metadata, 4, 2, sample_column, ['scientific_name', 'tax_id'], 'scientific_name(tax_id)', plain_text_bool)
        #Show advise message is warnings are detected
        if warnings_sacc_B24 == True or warnings_scol_B24 == True:
            #Show advise message
            print(rich_text_colored('\nThis could be due to:', 'due_to_header', plain_text_bool))
            print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', plain_text_bool))
            print(rich_text_colored('You should:', 'you_should_header', plain_text_bool))
            print('- Manually confirm that this is your case')
            print('- Check the samples with mixed scientific_name(tax_id)s')
            print('- Check the original publication and supplementary tables to get some context')
            print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #4.2.5) Check if there is more than one instrument_model(instrument_platform) per sample
        messageB25 = '\n5) Is there more than one instrument_model(instrument_platform) per sample?'
        print(rich_text_colored(messageB25, 'subsection', plain_text_bool))
        ##Check sample_accession column
        warnings_sacc_B25 = sample_vs_column(metadata, 5, 1,'sample_accession', ['instrument_model', 'instrument_platform'], 'instrument_model(instrument_platform)', plain_text_bool)
        ##Check provided sample_column
        warnings_scol_B25 = sample_vs_column(metadata, 5, 2, sample_column, ['instrument_model', 'instrument_platform'], 'instrument_model(instrument_platform)', plain_text_bool)
        #Show advise message is warnings are detected
        if warnings_sacc_B25 == True or warnings_scol_B25 == True:
            #Show advise message
            print(rich_text_colored('\nThis could be due to:', 'due_to_header', plain_text_bool))
            print('- The use of different technologies for sequencing the same sample', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print(rich_text_colored('\nYou should:', 'you_should_header', plain_text_bool))
            print('- Manually confirm that this is your case')
            print('- Check the samples with mixed instrument_model(instrument_platform)s')
            print('- Check the original Fastq names and content (if available)')
            print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')
            print('- Check the original publication and supplementary tables to get some context')
            print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #4.2.6) Check if there is more than one library_layout per sample
        messageB26 = '\n6) Is there more than one library_layout per sample?'
        print(rich_text_colored(messageB26, 'subsection', plain_text_bool))
        ##Check sample_accession column
        warnings_sacc_B26 = sample_vs_library_layout(metadata, 6, 1, 'sample_accession', plain_text_bool)
        ##Check provided sample_column
        warnings_scol_B26 = sample_vs_library_layout(metadata, 6, 2, sample_column, plain_text_bool)
        #Show advise message is warnings are detected
        if warnings_sacc_B26 == True or warnings_scol_B26 == True:
            #Show advise message
            print(rich_text_colored('\nThis could be due to:', 'due_to_header', plain_text_bool))
            print('- The presence of different data types in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of quality controlled Fastqs in the dataset', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- The presence of different sequencing technologies', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', plain_text_bool))
            print(rich_text_colored('\nYou should:','you_should_header', plain_text_bool))
            print('- Manually confirm which is your case (evaluate with the previous sample checks)')
            print('- Check the samples with mixed library_layouts')
            print('- Check the original Fastq names and content (if available)')
            print('- If necessary, search for extra information in the original database (If the dataset was not originally uploaded to ENA, try at Sequence Read Archive, or DNA Data Bank of Japan)')
            print('- Check the original publication and supplementary tables to get some context')
            print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #5)Show Legend if warnings were detected
        
        #Set functions warnings (Trues -> Warning detected)
        functions_warnings = (runs_layout_checks, submitted_files_check, dup_names_submitted_files_check,
                              warnings_sacc, warnings_scol, warnings_sacc_B23, warnings_scol_B23,
                              warnings_sacc_B24, warnings_scol_B24, warnings_sacc_B25, warnings_scol_B25,
                              warnings_sacc_B26, warnings_scol_B26)
        
        #Set main program warnings (Falses -> Warning detected)
        main_program_warnings = (run_sample_accession, run_sample_column, scol_equals_sacc)

        #Legend
        if (False in main_program_warnings) or (any(functions_warnings) == True):
            show_advise_legend(plain_text_bool)
        
        #6)Extra columns stats (Optional parameter)
        ##Check again if there were any extra_columns given
        if type(extra_columns_names) == list:
            #Show section header
            print(rich_text_colored('\nExtra Stats:', 'section_header', plain_text_bool))
            #For each extra_column get stats
            for extra_col in extra_columns_names:
                extra_stats_column(metadata, sample_column, extra_col, plain_text_bool)
    
    finally:
        #Print empty line for aesthetic purposes
        print('')
            
if __name__ == '__main__':
    main()
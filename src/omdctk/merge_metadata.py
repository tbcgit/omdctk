#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: merge_metadata.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'merge_metadata.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, OMD_CTK_Exception, program_header, 
                    check_headers, treat_output_directory_parameter_outfiles, 
                    print_list_n_byline, check_existence_directory_parameter,
                    rich_text_colored, show_advise_legend)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
import pandas as pd

#Program functions
def check_ena_metadata_file_extension(ena_table_path):
    """
    This function checks if the provided ena_table_path has the
    expected file extension for ENA Metadata Table.

    Parameters
    ----------
    ena_table_path : str
        Path for the ENA Metadata Table file.

    Raises
    ------
    OMD_CTK_Exception
        If the file does not have the expected extension raises an exception.
        This is done to avoid problems in determining the study accession from
        the ENA Metadata file name.

    Returns
    -------
    None.

    """
    if ena_table_path.endswith('_ENA_metadata.tsv'):
        pass
    else:
        raise OMD_CTK_Exception('Error! Expected file extension for ENA Metadata Table is "_ENA_metadata.tsv"!')


def treat_merge_columns(ena_merge_col_check, publication_merge_col_check):
    """
    This function treats the merge columns checks.

    Parameters
    ----------
    ena_merge_col_check : bool
        True: if ena_merge_col present in the ena_table.
        False: if ena_merge_col is not present in the ena_table.
    publication_merge_col_check : bool
        True: if publication_merge_col present in the publication_table.
        False: if publication_merge_col is not present in the publication_table.

    Raises
    ------
    OMD_CTK_Exception
        If all or any merge columns are not in their respective metadata tables
        raises an exception.

    Returns
    -------
    None.

    """
    if ena_merge_col_check == True and publication_merge_col_check == True:
        pass
    elif (ena_merge_col_check == True) and (publication_merge_col_check == False):
        raise OMD_CTK_Exception('Error! The provided Publication Merge Column is not in the Publication Metadata Table!\n Check the provided --publication_merge_column parameter!')
    elif (ena_merge_col_check == False) and (publication_merge_col_check == True):
        raise OMD_CTK_Exception('Error! The provided ENA Merge Column is not in the ENA Metadata Table!\n Check the provided --ena_merge_column parameter!')
    else:
        raise OMD_CTK_Exception('Error! Both of the provided merge columns are not in their respective metadata tables!\n Check merge columns parameters!')
        

def merge_columns_intersection_checks(ena_merge_column_values, publication_merge_column_values, color_treatment):
    """
    This function checks intersections between the provided merge column values
    and shows the corresponding messages if needed.

    Parameters
    ----------
    ena_merge_column_values : list
        List of values of the provided merge column for the ENA Metadata Table.
    publication_merge_column_values : list
        List of values of the provided merge column for the Publication Metadata Table.
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
    ena_merge_column_values_unique = set(ena_merge_column_values)
    publication_merge_column_values_unique = set(publication_merge_column_values)
    
    #Show stats
    print(rich_text_colored('\nUnique values:', 'subsection2', color_treatment))
    print('  o', rich_text_colored('ENA Merge Column (total unique values):', 'general_text', color_treatment), len(ena_merge_column_values_unique))
    print('  o', rich_text_colored('Publication Merge Column (total unique values):', 'general_text', color_treatment), len(publication_merge_column_values_unique))
    
    #Check if all unique values between merge columns are common
    all_ena_in_publication = ena_merge_column_values_unique.issubset(publication_merge_column_values_unique)
    all_publication_in_ena = publication_merge_column_values_unique.issubset(ena_merge_column_values_unique)
    if (all_ena_in_publication and all_publication_in_ena):
        print(rich_text_colored('\nAll unique values are common between the merge columns provided!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        ##Get lits of non-common values
        ena_difference = list(ena_merge_column_values_unique.difference(publication_merge_column_values_unique))
        publication_difference = list(publication_merge_column_values_unique.difference(ena_merge_column_values_unique))
        
        ##Show intersection stats
        print(rich_text_colored('\nIntersections:', 'subsection2', color_treatment))
        print('  o', rich_text_colored('Number of unique common values between merge columns:', 'general_text', color_treatment), len(ena_merge_column_values_unique.intersection(publication_merge_column_values_unique)))
        print('  o', rich_text_colored('Number of unique specific values for ENA Merge Column:', 'general_text', color_treatment), len(ena_difference))
        print('  o', rich_text_colored('Number of unique specific values for Publication Merge Column:', 'general_text', color_treatment), len(publication_difference))
        
        ##Show warning message
        print(rich_text_colored('\nWarning! Some values are not common between merge columns!', 'program_warning', color_treatment))
        
        #Show values if there are warnings
        ##ena_difference
        if len(ena_difference) > 0:
            print('\n-', rich_text_colored('Specific values for ENA Merge Column:', 'general_text', color_treatment))
            print_list_n_byline(ena_difference, 5)
        ##publication_difference
        if len(publication_difference) > 0:
            print('\n-', rich_text_colored('Specific values for Publication Merge Column:', 'general_text', color_treatment))
            print_list_n_byline(publication_difference, 5)
        
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:','due_to_header', color_treatment))
        print('- The presence of extra files or samples not used in the publication (values found only in the ENA Merge Column)', rich_text_colored('[Acceptable]', 'acceptable', color_treatment))
        print('- The absense of files or samples from the publication (values found only in the Publication Merge Column)', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm which is your case')
        print('- Check the original publication and supplementary tables to get some context')
        print('- If necessary, contact the authors of the original publication and/or the ENA Support (https://www.ebi.ac.uk/ena/browser/support)')
        
        #Return info for lately treat advise messages
        return True


#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter ena_metadata_table
    parser.add_argument(
            '-e','--ena_metadata_table', 
            action = 'store',
            required = True,
            help= 'ENA Metadata Table [Expected sep=TABS]. Indicate the path to the ENA Metadata Table file.'
    )
    ##Parameter publication_metadata_table
    parser.add_argument(
            '-p','--publication_metadata_table', 
            action = 'store',
            required = True,
            help = 'Publication Metadata Table [Expected sep=TABS]. Indicate the path to the Publication Metadata Table file.'
    )
    ##Parameter ena_merge_column
    parser.add_argument(
            '-ec','--ena_merge_column', 
            action = 'store',
            required = True,
            help = 'ENA Merge Column. ENA Metadata Table column to be used for merging.'
    )
    ##Parameter publication_merge_column
    parser.add_argument(
            '-pc','--publication_merge_column', 
            action = 'store',
            required = True,
            help = 'Publication Merge Column. Publication Metadata Table column to be used for merging.'
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
    ena_table_path = args.ena_metadata_table
    publication_table_path = args.publication_metadata_table
    ena_merge_col = args.ena_merge_column
    publication_merge_col = args.publication_merge_column
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
        #0)Initial checks
        
        #Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #Check ENA Metadata file extension
        check_ena_metadata_file_extension(ena_table_path)
        
        #1)Try to load files
        
        #Section header message
        print(rich_text_colored('\nLoading Files:\n', 'section_header', plain_text_bool))
        
        #Try to load ENA Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('ENA Metadata Table file:', 'general_text', plain_text_bool))
        print(ena_table_path)
        ##Try to load the ena table file
        ena_table = pd.read_csv(ena_table_path, sep = '\t')
        
        #Try to load Publication Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('\nPublication Metadata Table file:', 'general_text', plain_text_bool))
        print(publication_table_path)
        ##Try to load the publication table file
        publication_table = pd.read_csv(publication_table_path, sep = '\t')
        
        #2)Check provided merge columns
        ena_merge_col_check = check_headers(ena_merge_col, ena_table)
        publication_merge_col_check = check_headers(publication_merge_col, publication_table)
        treat_merge_columns(ena_merge_col_check, publication_merge_col_check)
        
        #3)Try to merge both tables
        
        #Section header messages
        print(rich_text_colored('\nCreating Merged Metadata Table:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nCombining tables:', 'general_text', plain_text_bool))
        print('Left Join using ENA Metadata Table as Reference')
        
        #Previous steps
        ##Get project accession from ena_table_path
        study_accession = ena_table_path.split('/')[-1].split('_ENA_metadata.tsv')[0]
        ##Get outfile name
        out_name_merged = study_accession + '_merged_metadata.tsv'
        ##Treat output_directory parameter / Get full output file path
        outputfile_merged = treat_output_directory_parameter_outfiles(out_name_merged, outputdir_path)
        
        #Join both tables by run accession (left join)
        df_merged_table=pd.merge(ena_table, publication_table, left_on = ena_merge_col, right_on = publication_merge_col, how = 'left')
        
        #Show saved file message
        print(rich_text_colored('\nSaving results in file:', 'general_text', plain_text_bool))
        print(outputfile_merged)
        #Save df_merged_table
        df_merged_table.to_csv(outputfile_merged, header = True, index = False, sep = '\t')
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
        
    else:
        #4)Check intersections between merge columns
        
        #Section header message
        print(rich_text_colored('\nMerge Columns Intersection Checks:', 'section_header', plain_text_bool))
        print(rich_text_colored("\nCheck merge columns' unique values:", 'general_text', plain_text_bool))
        print(rich_text_colored('ENA Merge Column selected:','column_color', plain_text_bool), ena_merge_col)
        print(rich_text_colored('Publication Merge Column selected:','column_color', plain_text_bool), publication_merge_col)
        
        #Get list of values for merge columns
        ena_merge_col_values = list(ena_table[ena_merge_col])
        publication_merge_col_values = list(publication_table[publication_merge_col])
        
        #Check intersections
        merge_cols_intersec_check = merge_columns_intersection_checks(ena_merge_col_values, publication_merge_col_values, plain_text_bool)
        
        #Show Legend if warnings were detected
        if merge_cols_intersec_check == True:
            show_advise_legend(plain_text_bool)
        
    finally:
        #Print empty line for aesthetic purposes
        print('')
        
if __name__ == '__main__':
    main()
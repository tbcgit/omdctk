#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: merge_metadata.py
@author: Samuel Piquer-Esteban
@date: 17 Feb 2024 

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
def treat_merge_columns(main_merge_col_check, extra_merge_col_check):
    """
    This function treats the merge columns checks.

    Parameters
    ----------
    main_merge_col_check : bool
        True: if main_merge_col present in the main_table.
        False: if main_merge_col is not present in the main_table.
    extra_merge_col_check : bool
        True: if extra_merge_col present in the extra_table.
        False: if extra_merge_col is not present in the extra_table.

    Raises
    ------
    OMD_CTK_Exception
        If all or any merge columns are not in their respective metadata tables
        raises an exception.

    Returns
    -------
    None.

    """
    if main_merge_col_check == True and extra_merge_col_check == True:
        pass
    elif (main_merge_col_check == True) and (extra_merge_col_check == False):
        raise OMD_CTK_Exception('Error! The provided Extra Metadata Merge Column is not in the Extra Metadata Table!\n Check the provided --extra_merge_column parameter!')
    elif (main_merge_col_check == False) and (extra_merge_col_check == True):
        raise OMD_CTK_Exception('Error! The provided Main Metadata Merge Column is not in the Main Metadata Table!\n Check the provided --main_merge_column parameter!')
    else:
        raise OMD_CTK_Exception('Error! Both of the provided merge columns are not in their respective metadata tables!\n Check merge columns parameters!')
        

def merge_columns_intersection_checks(main_merge_column_values, extra_merge_column_values, color_treatment):
    """
    This function checks intersections between the provided merge column values
    and shows the corresponding messages if needed.

    Parameters
    ----------
    main_merge_column_values : list
        List of values of the provided merge column for the Main Metadata Table.
    extra_merge_column_values : list
        List of values of the provided merge column for the Extra Metadata Table.
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
    main_merge_column_values_unique = set(main_merge_column_values)
    extra_merge_column_values_unique = set(extra_merge_column_values)
    
    #Show stats
    print(rich_text_colored('\nUnique values:', 'subsection2', color_treatment))
    print('  o', rich_text_colored('Main Metadata Merge Column (total unique values):', 'general_text', color_treatment), len(main_merge_column_values_unique))
    print('  o', rich_text_colored('Extra Metadata Merge Column (total unique values):', 'general_text', color_treatment), len(extra_merge_column_values_unique))
    
    #Check if all unique values between merge columns are common
    all_main_in_extra = main_merge_column_values_unique.issubset(extra_merge_column_values_unique)
    all_extra_in_main = extra_merge_column_values_unique.issubset(main_merge_column_values_unique)
    if (all_main_in_extra and all_extra_in_main):
        print(rich_text_colored('\nAll unique values are common between the merge columns provided!', 'acceptable', color_treatment))
        #Return info for lately treat advise messages
        return False
    else:
        ##Get lits of non-common values
        main_difference = list(main_merge_column_values_unique.difference(extra_merge_column_values_unique))
        extra_difference = list(extra_merge_column_values_unique.difference(main_merge_column_values_unique))
        
        ##Show intersection stats
        print(rich_text_colored('\nIntersections:', 'subsection2', color_treatment))
        print('  o', rich_text_colored('Number of unique common values between merge columns:', 'general_text', color_treatment), len(main_merge_column_values_unique.intersection(extra_merge_column_values_unique)))
        print('  o', rich_text_colored('Number of unique specific values for Main Metadata Merge Column:', 'general_text', color_treatment), len(main_difference))
        print('  o', rich_text_colored('Number of unique specific values for Extra Metadata Merge Column:', 'general_text', color_treatment), len(extra_difference))
        
        ##Show warning message
        print(rich_text_colored('\nWarning! Some values are not common between merge columns!', 'program_warning', color_treatment))
        
        #Show values if there are warnings
        ##main_difference
        if len(main_difference) > 0:
            print('\n-', rich_text_colored('Specific values for Main Metadata Merge Column:', 'general_text', color_treatment))
            print_list_n_byline(main_difference, 5)
        ##extra_difference
        if len(extra_difference) > 0:
            print('\n-', rich_text_colored('Specific values for Extra Metadata Merge Column:', 'general_text', color_treatment))
            print_list_n_byline(extra_difference, 5)
        
        #Show advise messages
        print(rich_text_colored('\nThis could be due to:','due_to_header', color_treatment))
        print('- The presence of additional files or unused samples', rich_text_colored('[Acceptable]', 'acceptable', color_treatment))
        print('- The absense of files or samples', rich_text_colored('[Warning]', 'legend_warning', color_treatment))
        print("- Authors' mishandle or upload errors", rich_text_colored('[Dangerous]', 'dangerous', color_treatment))
        print(rich_text_colored('\nYou should:', 'you_should_header', color_treatment))
        print('- Manually confirm which is your case')
        print('- Check the original publication and supplementary tables to get some context')
        print('- If necessary, contact the authors of the original publication')
        
        #Return info for lately treat advise messages
        return True


#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter main_metadata_table
    parser.add_argument(
            '-m','--main_metadata_table', 
            action = 'store',
            required = True,
            help= 'Main Metadata Table [Expected sep=TABS]. Indicate the path to the Main Metadata Table file.'
    )
    ##Parameter extra_metadata_table
    parser.add_argument(
            '-e','--extra_metadata_table', 
            action = 'store',
            required = True,
            help = 'Extra Metadata Table [Expected sep=TABS]. Indicate the path to the Extra Metadata Table file.'
    )
    ##Parameter main_merge_column
    parser.add_argument(
            '-mc','--main_merge_column', 
            action = 'store',
            required = False,
            default = None,
            help = 'Main Metadata Merge Column. Main Metadata Table column to be used for merging. This parameter will be skipped if pandas_merge_mode = cross.'
    )
    ##Parameter extra_merge_column
    parser.add_argument(
            '-ec','--extra_merge_column', 
            action = 'store',
            required = False,
            default = None,
            help = 'Extra Metadata Merge Column. Extra Metadata Table column to be used for merging. This parameter will be skipped if pandas_merge_mode = cross.'
    )
    ##Parameter pandas_merge_mode
    parser.add_argument(
            '-p','--pandas_merge_mode', 
            action = 'store',
            choices = ['left','right','outer','inner','cross'],
            required = False,
            default = 'left',
            help = 'Pandas Merge Mode (Optional) [Default:left]. Indicate the pandas merge mode to be used for merging.'
    )
    ##Parameter main_merge_suffix
    parser.add_argument(
            '-ms','--main_merge_suffix', 
            action = 'store',
            required = False,
            default = '_x',
            help = 'Main Metadata Pandas Merge Suffix (Optional) [Default:"_x"]. Suffix to add to overlapping column names for the Main Metadata columns.'
    )
    ##Parameter extra_merge_suffix
    parser.add_argument(
            '-es','--extra_merge_suffix', 
            action = 'store',
            required = False,
            default = '_y',
            help = 'Extra Metadata Pandas Merge Suffix (Optional) [Default:"_y"]. Suffix to add to overlapping column names for the Extra Metadata columns.'
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
    main_table_path = args.main_metadata_table
    extra_table_path = args.extra_metadata_table
    main_merge_col = args.main_merge_column
    extra_merge_col = args.extra_merge_column
    pd_merge_mode = args.pandas_merge_mode
    main_cols_suffix = args.main_merge_suffix
    extra_cols_suffix = args.extra_merge_suffix
    outputdir_path = args.output_directory
    plain_text_bool = args.plain_text
    #Skip main_merge_column and extra_merge_column if pandas_merge_mode = cross (this is only for aesthetic purposes when showing program headers)
    if pd_merge_mode == 'cross':
        args.main_merge_column = None
        args.extra_merge_column = None
    
    #Parser error for -mc and -ec if any are none in a merge mode different than 'cross'
    if pd_merge_mode != 'cross' and (args.main_merge_column is None or args.extra_merge_column is None):
        parser.error('the following arguments are required: -mc/--main_merge_column, -ec/--extra_merge_column for this particular merge mode')
    
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
        
        #1)Try to load files
        
        #Section header message
        print(rich_text_colored('\nLoading Files:\n', 'section_header', plain_text_bool))
        
        #Try to load Main Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('Main Metadata Table file:', 'general_text', plain_text_bool))
        print(main_table_path)
        ##Try to load the ena table file
        main_table = pd.read_csv(main_table_path, sep = '\t')
        
        #Try to load Extra Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('\nExtra Metadata Table file:', 'general_text', plain_text_bool))
        print(extra_table_path)
        ##Try to load the publication table file
        extra_table = pd.read_csv(extra_table_path, sep = '\t')
        
        #2)Check provided merge columns if pd_merge_mode is not cross
        if pd_merge_mode != 'cross':
                main_merge_col_check = check_headers(main_merge_col, main_table)
                extra_merge_col_check = check_headers(extra_merge_col, extra_table)
                treat_merge_columns(main_merge_col_check, extra_merge_col_check)
        
        #3)Try to merge both tables
        
        #Section header messages
        print(rich_text_colored('\nCreating Merged Metadata Table:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nCombining tables:', 'general_text', plain_text_bool))
        print(rich_text_colored('Pandas Merge Mode:', 'merge_mode_color', plain_text_bool), pd_merge_mode)
        print('  o',rich_text_colored('Left Table (x):', 'merge_table_color', plain_text_bool), 'Main Metadata Table')
        print('  o',rich_text_colored('Right Table (y):', 'merge_table_color', plain_text_bool), 'Extra Metadata Table')
        
        #Previous steps
        ##Get project accession from main_table_path
        main_file_name = main_table_path.split('/')[-1]
        ##Get outfile name
        out_name_merged = 'merged_' + main_file_name
        
        ##Treat output_directory parameter / Get full output file path
        outputfile_merged = treat_output_directory_parameter_outfiles(out_name_merged, outputdir_path)
        
        #Join both tables
        if pd_merge_mode == 'cross':
            df_merged_table = pd.merge(main_table, extra_table, how = pd_merge_mode, suffixes = (main_cols_suffix, extra_cols_suffix))
        else:
            df_merged_table = pd.merge(main_table, extra_table, left_on = main_merge_col, right_on = extra_merge_col, how = pd_merge_mode, suffixes = (main_cols_suffix, extra_cols_suffix))
        
        #Show saved file message
        print(rich_text_colored('\nSaving results in file:', 'general_text', plain_text_bool))
        print(outputfile_merged)
        #Save df_merged_table
        df_merged_table.to_csv(outputfile_merged, header = True, index = False, sep = '\t')
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
        
    else:
        #4)Check intersections between merge columns if pd_merge_mode is not cross
        if pd_merge_mode != 'cross':
            #Section header message
            print(rich_text_colored('\nMerge Columns Intersection Checks:', 'section_header', plain_text_bool))
            print(rich_text_colored("\nCheck merge columns' unique values:", 'general_text', plain_text_bool))
            print(rich_text_colored('Main Metadata Merge Column selected:','column_color', plain_text_bool), main_merge_col)
            print(rich_text_colored('Extra Metadata Merge Column selected:','column_color', plain_text_bool), extra_merge_col)
        
            #Get list of values for merge columns
            main_merge_col_values = list(main_table[main_merge_col])
            extra_merge_col_values = list(extra_table[extra_merge_col])
        
            #Check intersections
            merge_cols_intersec_check = merge_columns_intersection_checks(main_merge_col_values, extra_merge_col_values, plain_text_bool)
        
            #Show Legend if warnings were detected
            if merge_cols_intersec_check == True:
                show_advise_legend(plain_text_bool)
        
    finally:
        #Print empty line for aesthetic purposes
        print('')
        
if __name__ == '__main__':
    main()
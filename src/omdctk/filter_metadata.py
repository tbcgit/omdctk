#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: filter_metadata.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'filter_metadata.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, OMD_CTK_Exception, program_header, 
                    check_existence_directory_parameter, 
                    check_headers, treat_headers_check,
                    check_values, treat_values_check,
                    treat_output_directory_parameter_outfiles,
                    rich_text_colored)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
from ast import literal_eval
import pandas as pd

#Program Constants
FILTERFILE_HEADERS = ['variable', 'values', 'filter_type', 'action', 'NA_treatment']

FILTERFILE_FILTER_TYPES_ACTION_COMB = ['categorical + keep', 'categorical + drop',
                                       'numerical + greater', 'numerical + greater_equal', 
                                       'numerical + equal', 'numerical + less',
                                       'numerical + less_equal']

FILTERFILE_NA_TREATMENT = ['drop', 'keep', 'no']

#Program Functions
def check_filterfile_values_column_format(filterfile_df, color_treatment):
    """
    This function checks the format of the provided values
    in the "values" column of the filter file.

    Parameters
    ----------
    filterfile_df : pandas dataframe
        The provided filterfile dataframe.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If format inconsistencies are detected raises an exception.

    Returns
    -------
    None.

    """
    #Check format(list) for values
    format_checks = []
    for val in filterfile_df['values']:
        try: 
            #Try to strip start and end spaces
            strip_val = val.strip()
            #Check list format
            if strip_val.startswith('[') and strip_val.endswith(']'):
                format_checks.append(True)
            else:
                format_checks.append(False)
        except:
            format_checks.append(False)
    
    #Check empty cells
    NAs_checks = list(filterfile_df['values'].isna())
    
    #Combine checks
    result_bools = []
    result_prints = []
    for i in range(len(format_checks)):
        #Wrong format and Not NA
        if format_checks[i] == False and NAs_checks[i] == False:
            result_bools.append(False)
            print_strn = ' '.join(['-', rich_text_colored(filterfile_df['variable'][i], 'check_color', color_treatment), rich_text_colored('Does not have the expected "[]" format!', 'program_warning', color_treatment)])
            result_prints.append(print_strn)
        #Wrong format and NA
        elif format_checks[i] == False and NAs_checks[i] == True:
            result_bools.append(False)
            print_strn = ' '.join(['-', rich_text_colored(filterfile_df['variable'][i], 'check_color', color_treatment), rich_text_colored('Cell is empty!', 'program_warning', color_treatment)])
            result_prints.append(print_strn)
        #Right format (It will never be NA)
        else:
            result_bools.append(True)
            print_strn = ' '.join(['-', rich_text_colored(filterfile_df['variable'][i], 'check_color', color_treatment), rich_text_colored('Everything seems okay!', 'acceptable', color_treatment)])
            result_prints.append(print_strn)
    
    #Treat final checks
    if False in result_bools:
        #Prepare exception frase
        frase0 = 'Error! Some of the values of the "values" column in the Filter Table are not valid!\n Check your filter file!\n'
        frase1 = '\nShowing checks by variable:'
        frase = [frase0, frase1]
        for i in result_prints:
            temp_line = ' '.join(['\n', i])
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


def process_filterfile_values_and_checks(variable, values, filter_type, metadata_df, color_treatment):
    """
    This function processes the provided values and runs various checks 
    depending on their filter_type. This function processes only one 
    input at a time. 

    Parameters
    ----------
    variable : str
        The provided variable to check(colname in metadata table).
    values : str
        The provided values to check.
    filter_type : str
        the filter_type to be applied.
        Expected values ['categorical','numerical']
    metadata_df : pandas dataframe
        The provided metadata dataframe to check with.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Returns
    -------
    print_strn : str
        Colored string message indicating if everything is okey (green) or
        if there is any warning to check (red).
    bool_value : bool
        Bool value indicating if everything is okey (True) or there is
        any warning to check (False).
    result_value : list
        The list of provided values processed to properly interact with the
        metadata table.
        
    """
    #Pre-process values until get a list of multiple or single values
    ##Convert curly to straight quotation(in case file has been edited in calc, or excel)
    pre_values = curly_2_straight_quotation(values)
    ##Try to evaluate string to convert it to list
    try:
        post_values = literal_eval(pre_values)
    except:
        print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('The provided values can not be converted to list!', 'program_warning', color_treatment)])
        bool_value = False
        result_value = pre_values
    else:
        #Check that the list is not empty
        if len(post_values) == 0:
            print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('The provided list is empty!', 'program_warning', color_treatment)])
            bool_value = False
            result_value = post_values
        else:
            #Type_filter checks
            ##Treat numerical type_filters values
            if filter_type == 'numerical':
                #Get metadata column dtype
                metadata_colum_dtype = metadata_df[variable].dtype
                ##If column in metadata table is numerical continue/ else stop
                if metadata_colum_dtype == int or metadata_colum_dtype == float:
                    ##Numerical filter values must be a unique number
                    if len(post_values) == 1:
                        #Check that the value is int or float
                        if type(post_values[0]) == int or type(post_values[0]) == float:
                            print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('Everything seems okay!', 'acceptable', color_treatment)])
                            bool_value = True
                            result_value = post_values
                        else:
                            print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('The provided value is not numerical!', 'program_warning', color_treatment)])
                            bool_value = False
                            result_value = post_values
                    else:
                        print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('Values for numerical filters must be a unique number!', 'program_warning', color_treatment)])
                        bool_value = False
                        result_value = post_values
                else:
                    print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('Is not a numerical variable in the Metadata Table!', 'program_warning', color_treatment)])
                    bool_value = False
                    result_value = post_values
            ##Treat categorical type_filters values
            else:
                #Check if ALL the provided categorical values are in the corresponding column in metadata table
                bool_value = check_values_inv(metadata_df[variable], post_values)
                result_value = post_values
                #If true continue/ else warning
                if bool_value == True:
                    print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('Everything seems okay!', 'acceptable', color_treatment)])
                else:
                    print_strn = ' '.join(['-', rich_text_colored(variable, 'check_color', color_treatment), rich_text_colored('Some of the provided values are not in the Metadata Table!', 'program_warning', color_treatment)])
            
    return print_strn, bool_value, result_value


def treat_process_filterfile_values_and_checks(filterfile_df, metadata_df, color_treatment):
    """
    This function treats the processed values checks and returns
    the processed information of the filterfile.

    Parameters
    ----------
    filterfile_df : pandas dataframe
        The provided filterfile dataframe to treat.
    metadata_df : pandas dataframe
        The provided metadata dataframe to check with.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If any warning is detected raises an exception.

    Returns
    -------
    variables_4_filtering : list
        List with "variable" column values from filterfile.
    values_4_filtering : list
        List with formated "values" column values from filterfile.
    filtertypes_4_filtering : list
        List with "filter_type" column values from filterfile.
    actions_4_filtering : list
        List with "action" column values from filterfile.
    na_treatment_4_filtering : list
        List with "NA_treatment" column values from filterfile.

    """
    #Check and process values for each variable in filter table
    
    ##Init list for results (internal function results)
    print_4_filter_values_check = []
    bools_4_filter_values_check = []
    ##Init lists for results (results to be returned)
    variables_4_filtering = []
    values_4_filtering = []
    filtertypes_4_filtering = []
    actions_4_filtering = []
    na_treatment_4_filtering = []
    
    ##For each variable
    for i,r in filterfile_df.iterrows():
        ##Get temp values
        temp_variable = r['variable']
        temp_values = r['values']
        temp_filtertype = r['filter_type']
        temp_action = r['action']
        temp_na_treatment = r['NA_treatment']
        
        ##Get temp result check and processed values
        temp_values_check = process_filterfile_values_and_checks(temp_variable, temp_values, temp_filtertype, metadata_df, color_treatment)
        
        ##Append results
        print_4_filter_values_check.append(temp_values_check[0])
        bools_4_filter_values_check.append(temp_values_check[1])
        values_4_filtering.append(temp_values_check[2])
        variables_4_filtering.append(temp_variable)
        filtertypes_4_filtering.append(temp_filtertype)
        actions_4_filtering.append(temp_action)
        na_treatment_4_filtering.append(temp_na_treatment)
    
    ##If/else
    ##If there are inconsistencies exception
    if False in bools_4_filter_values_check:
        #Prepare exception frase
        frase0 = 'Warning! Inconsistencies detected in the "values" column of the Filter Table!\n Check your filter file!\n'
        frase1 = '\nShowing checks by variable:'
        frase = [frase0, frase1]
        for i in range(len(print_4_filter_values_check)):
            frase.append(' '.join(['\n\n' + print_4_filter_values_check[i]]))
            frase.append(' '.join(['\n  Filter Type:', filtertypes_4_filtering[i]]))
            frase.append(' '.join(['\n  Filter Values:', str(values_4_filtering[i])]))
            frase.append(' '.join(['\n  Metadata Column dtype:', str(metadata_df[variables_4_filtering[i]].dtype)]))
        #exception
        raise OMD_CTK_Exception(' '.join(frase))

    ##Else continue an return external results
    else:
        return variables_4_filtering, values_4_filtering, filtertypes_4_filtering, actions_4_filtering, na_treatment_4_filtering
    

def filter_metadata_variable(variable, processed_values, filter_type, filter_applied, NA_treatment, metadata_df_cp):
    """
    This function filters the provided metadata copy with the provided values 
    depending on their filter_type and filter. This function processes only 
    one input at a time.
    
    Parameters
    ----------
    variable : str
        The provided variable to filter the metadata table copy.
    processed_values : list
        The provided list of processed values used to filter the metadata table copy.
    filter_type : str
        The filter_type to be applied.
        Expected values ['categorical', 'numerical']
    filter_applied : str
        The filter action to be applied.
        Expected values for numerical ['equal', 'greater', 'greater_equal', 'less', 'less_equal']
        Expected values for categorical ['keep', 'drop']
    NA_treatment : str
        The NA treatment to be applied.
        Expected values ['keep', 'drop', 'no']
    metadata_df_cp : pandas dataframe
        The copy of the original metadata table to be filtered.

    Returns
    -------
    metadata_df_cp_f : pandas dataframe
        A filtered version of the metadata dataframe copy.
        
    """
    #Construct NA expression
    if NA_treatment == 'keep':
        #' | `variable`.isna()'
        NA_expression = ''.join([' | `',variable, '`.isna()'])
    elif NA_treatment == 'drop':
        #' & (~`variable`.isna())'
        NA_expression = ''.join([' & (~`',variable, '`.isna())'])
    #Else should be 'no'
    else:
        pass
        
    #Construct filter expression
    ##Treat numerical values
    if filter_type == 'numerical':
        if filter_applied == 'greater':
            #'`variable` > processed_values[0]'
            filter_expression = ''.join(['`',variable, '` > ', str(processed_values[0])])
        elif filter_applied == 'greater_equal':
            #'`variable` >= processed_values[0]'
            filter_expression = ''.join(['`', variable, '` >= ', str(processed_values[0])])
        elif filter_applied == 'less':
            #'`variable` < processed_values[0]'
            filter_expression = ''.join(['`', variable, '` < ', str(processed_values[0])])
        elif filter_applied == 'less_equal':
            #'`variable` <= processed_values[0]'
            filter_expression = ''.join(['`', variable, '` <= ', str(processed_values[0])])
        #Else should be equal
        else:
            #'`variable` == processed_values[0]'
            filter_expression = ''.join(['`', variable, '` == ', str(processed_values[0])])
    ##Treat categorical values
    else:
        #Treat different filter options
        if filter_applied == 'keep':
            #Keep rows that match values in processed_values for the variable
            #'`variable`.isin(processed_values)'
            filter_expression = ''.join(['`', variable, '`.isin(', str(processed_values),')'])
        else:
            #Else should be 'drop'
            #Drop rows that match values in processed_values for the variable
            #'(~`variable`.isin(processed_values))'
            filter_expression = ''.join(['(~`', variable, '`.isin(', str(processed_values),'))'])
            
    #Combine expressions and apply filter with query
    if NA_treatment == 'no':
        metadata_df_cp_f = metadata_df_cp.query(filter_expression)
    else:
        metadata_df_cp_f = metadata_df_cp.query(filter_expression + NA_expression)

    #Return filtered table
    return metadata_df_cp_f


def treat_filter_metadata_table(processed_filterfile, metadata_df, color_treatment):
    """
    This function treats the different filters from the filterfile
    and provides the filtered metadata table.

    Parameters
    ----------
    processed_filterfile : tuple of lists
        The treat_process_filterfile_values_and_checks() results.
        Processed filterfile.
    metadata_df : pandas dataframe
        The provided metadata dataframe.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If the filtered metadata table is empty after applying filters raises
        an exception.

    Returns
    -------
    filtered_metadata : pandas dataframe
        The filtered metadata dataframe.

    """   
    ##Init a copy of the metadata table for filtering
    filtered_metadata = metadata_df
    
    ##Filter with the provided variables
    for i in range(len(processed_filterfile[0])):
        ##Get temp values
        temp_variable = processed_filterfile[0][i]
        temp_values = processed_filterfile[1][i]
        temp_filtertype = processed_filterfile[2][i]
        temp_action = processed_filterfile[3][i]
        temp_na_treatment = processed_filterfile[4][i]
        
        ##Print temp filter information
        print(rich_text_colored('\nFilter Variable: ', 'general_text', color_treatment), temp_variable)
        print(rich_text_colored('Values: ', 'general_text', color_treatment), temp_values)
        print(rich_text_colored('Filter Type: ', 'filter_color', color_treatment), temp_filtertype)
        print(rich_text_colored('Filter Action: ', 'filter_color', color_treatment), temp_action)
        print(rich_text_colored('NA Treatment: ', 'filter_color', color_treatment), temp_na_treatment)
        
        ##Apply filter
        filtered_metadata = filter_metadata_variable(temp_variable, temp_values, temp_filtertype, temp_action, temp_na_treatment, filtered_metadata)
        
        ##Check if the metadata copy is empty after filtering
        if filtered_metadata.empty == True:
            #exception
            raise OMD_CTK_Exception('Warning! The Metadata Table became empty after this last filter!\n Check your filter file and variables!')
            
    ##Return filtered metadata table if no exceptions
    return filtered_metadata

     
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
    ##Parameter filter_table
    parser.add_argument(
            '-f','--filter_table', 
            action = 'store',
            required = True ,
            help = 'Filter Table [Expected sep=TABS]. Indicate the path to the Filter Table file. See Documentation for more information and format details.'
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
    metadata_table_path = args.metadata_table
    filterfile_path = args.filter_table
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
        
        #1)Try to load files
        
        #Section header message
        print(rich_text_colored('\nLoading Files:\n', 'section_header', plain_text_bool))
        
        #Try to load Metadata Table as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('Metadata Table file:', 'general_text', plain_text_bool))
        print(metadata_table_path)
        ##Try to load the metadata table file
        metadata = pd.read_csv(metadata_table_path, sep = '\t')
        
        #Try to load Filter File as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('\nFilter Table file:', 'general_text', plain_text_bool))
        print(filterfile_path)
        ##Try to load the ena table file
        filterfile = pd.read_csv(filterfile_path, sep = '\t', dtype = 'string')
        ##Check headers and treatment
        check_headers_filterfile = check_headers(FILTERFILE_HEADERS, filterfile)
        frase0_3 = 'Error! Some of the needed headers are not in the Filter Table!\n Check your filter file!\n'
        frase1_3 = '\nThe headers needed are:'
        treat_headers_check(FILTERFILE_HEADERS, check_headers_filterfile, filterfile, frase0_3, frase1_3, plain_text_bool)
        
        #2)Check that values in "variable" column from Filter File are headers in Metadata Table
        
        #Get list of "variable" column values in filter file
        provided_variables = list(filterfile['variable'])
        
        #Check headers in metadata table and treatment
        check_provided_variables = check_headers(provided_variables, metadata)
        frase0_4 = 'Error! Some of the provided variables in the Filter Table do not match the Metadata Table!\n Check your filter file and/or metadata file!\n'
        frase1_4 = '\nProvided variables are:'
        treat_headers_check(provided_variables, check_provided_variables, metadata, frase0_4, frase1_4, plain_text_bool)
        
        #3)Check that "filter_type" and "action" column combinations are valid
        
        #Get "filter_type" and "action" column combinations in filter file
        type_sum_filter = filterfile['filter_type'] + ' + ' + filterfile['action']
        
        #Check that combinations are valid and treatment
        check_expected_type_sum_filter_values = check_values(FILTERFILE_FILTER_TYPES_ACTION_COMB, type_sum_filter)
        frase0_5 = 'Error! Some of the "filter_type" to "action" combinations in the Filter Table are not valid!\n Check your filter file!\n'
        frase1_5 = '\nValid combinations are:'
        treat_values_check(FILTERFILE_FILTER_TYPES_ACTION_COMB, check_expected_type_sum_filter_values, frase0_5, frase1_5, plain_text_bool)
        
        #4)Check that values in "NA_treatment" column from Filter File are valid and treatment
        check_expected_NA_treatment_values = check_values(FILTERFILE_NA_TREATMENT, filterfile['NA_treatment'])
        frase0_6 = 'Error! Some of the values of the "NA_treatment" column in the Filter Table are not valid!\n Check your filter file!\n'
        frase1_6 = '\nValid values are:'
        treat_values_check(FILTERFILE_NA_TREATMENT, check_expected_NA_treatment_values, frase0_6, frase1_6, plain_text_bool)
        
        #5)Check that values in "values" column from Filter File has the right format and treatment
        check_filterfile_values_column_format(filterfile, plain_text_bool)
        
        #6)Check and process values in "values" column from Filter File. And treatment of checks.
        processed_filterfile = treat_process_filterfile_values_and_checks(filterfile, metadata, plain_text_bool)
        
        #7)Filtering metadata table
        
        #Previous steps
        ##Get outfile name
        outfile_name = 'filtered_'+metadata_table_path.split('/')[-1]
        ##Treat output_directory parameter / Get full output file path
        outputfile_path = treat_output_directory_parameter_outfiles(outfile_name, outputdir_path)
        
        #Show section header
        print(rich_text_colored('\nFiltering Metadata Table:', 'section_header', plain_text_bool))
        
        #Get filtered metadata table
        filtered_metadata = treat_filter_metadata_table(processed_filterfile, metadata, plain_text_bool)
        
        #Show saved file message
        print(rich_text_colored('\nSaving results in file:', 'general_text', plain_text_bool))
        print(outputfile_path)
        #Save filtered_metadata
        filtered_metadata.to_csv(outputfile_path, header = True, index = False, sep = '\t')
    
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
           
    finally:
        #Print empty line for aesthetic purposes
        print('')

if __name__ == '__main__':
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: check_metadata_values.py
@author: Samuel Piquer-Esteban
@date: 1 Mar 2024 

"""

#Set program name
__program__ = 'check_metadata_values.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, VALID_REQUIREDNESS, check_values_inv, print_list_n_byline, 
                    show_advise_legend, program_header, treat_headers_check, curly_2_straight_quotation, 
                    check_headers, check_duplicates_in_variables_dict_column, 
                    check_required_variables_in_metadata_table, check_metadata_table_vars_in_dict, 
                    check_na_in_pandas_dataframe, treat_values_check, check_values, rich_text_colored)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
from ast import literal_eval
import pandas as pd

#Program Constants
VARIABLES_FILE_HEADERS = ['variable', 'requiredness', 'class_type','uniqueness_within_variable', 'check_uniqueness_between_variables', 'variables_for_uniqueness_between', 'NAs_allowed', 'allowed_values_treatment', 'allowed_values']

VALID_CLASS_TYPE = ['character', 'numeric']

VALID_UNIQUENESS = ['unique', 'nonunique']

VALID_CHECK_UNIQUENESS_BETWEEN_OPTIONS = ['yes', 'no']

VALID_UNIQUENESS_BETWEEN_VARIABLES = ['the "none" string', 'a non-empty variables list formatted as a Python list (only variables present in the "variable" column different than the variable been evaluated are allowed)']

VALID_UNIQUENESS_BETWEEN_VARIABLES_COMB = ['yes + a non-empty variables list formatted as a Python list (only variables present in the "variable" column different than the variable been evaluated are allowed)', 
                                           'no + the "none" string']

VALID_NAS_ALLOWED = ['yes', 'no']

VALID_ALLOWED_VALUES_TREATMENT = ['wholeset', 'subset','range', 'any']

VALID_ALLOWED_VALUES = ['the "any" string', 'a non-empty list formatted as a Python list']

VALID_CLASSTYPE_TREATMENT_COMB = ['character + wholeset', 'character + subset', 'character + any','numeric + wholeset', 
                                  'numeric + subset', 'numeric + any', 'numeric + range']

VALID_REQUIREDNESS_NAS_COMB = ['required + NAs_allowed:no', 'optional + NAs_allowed:yes', 'optional + NAs_allowed:no']

VALID_TREATMENT_ALOWED_VALUES_COMB = ['wholeset + a non-empty list formatted as a Python list', 'subset + a non-empty list formatted as a Python list', 
                                      'range + a non-empty numeric [min,max] list formatted as a Python list', 'any + the "any" string']

#Program functions
def define_allowed_values_column_type(value):
    """
    This function checks the provided string from allowed values column and
    defines the type of value provided.

    Parameters
    ----------
    value : str
        String from allowed_values column.

    Returns
    -------
    result : str
        any -> any string indicating that any value is allowed
        NA_value -> the cell is empty
        not_list -> the cell is filled but not a Python list
        empty_list -> the cell is filled with an empty Python list
        filled_list -> the cell is filled with a non-empty Python list
        filled_minmax_numeric_list -> the cell the cell is filled with a non-empty Python nuemric list of len=2 and [min,max]

    """
    if value == 'any':
        result = 'any'
    elif pd.isna(value):
        result = 'NA_value'
    else:
        #Try to convert
        ##Try to evaluate string to convert it to list
        try:
            evaluated_values = literal_eval(curly_2_straight_quotation(value))
        except:
            result = 'not_list'
        else:
            #Check if is list
            if not (type(evaluated_values) == list):
                result = 'not_list'
            else:
                #Check if list is not empty
                if len(evaluated_values) == 0:
                    result = 'empty_list'
                else:
                    #Check if all elements are numeric
                    if all((type(x) == int or type(x) == float) for x in evaluated_values):
                        #Check len=2 for 1st range criteria(numeric len=2 list)
                        if len(evaluated_values) == 2:
                            #Check evaluated_values[0]<evaluated_values[1] for 2nd range criteria([min,max] range list)
                            if evaluated_values[0] < evaluated_values[1]:
                                result = 'filled_minmax_numeric_list'
                            else:
                                result = 'filled_list'
                        else:
                            result = 'filled_list'
                    else:
                        result = 'filled_list'
    #Return allowed value type
    return result


def check_allowed_values_column(variables_dict_table):
    """
    This function checks if all cells in allowed_values column from
    Variables Dictionary file are valid values.

    Parameters
    ----------
    variables_dict_table : pandas dataframe
        The provided variables dictionary dataframe.

    Returns
    -------
    bool
        True: all values are valid.
        False: not all values are valid.

    """
    
    #Get list of provided allowed values
    allowed_values = list(variables_dict_table['allowed_values'])
    
    #Init check_results values(bools)
    results = []
    
    #Check each provided cell value
    for val in allowed_values:
        #Get temp value type
        temp_value_type = define_allowed_values_column_type(val)
        
        #Set value based on cell_type
        if temp_value_type == 'any' or temp_value_type =='filled_list' or temp_value_type =='filled_minmax_numeric_list':
            results.append(True)
        else:
            results.append(False)

    #Check if there are False values and return bool
    return all(results)


def check_allowed_values_treatment_comb(variables_dict_table):
    """
    This function checks the allowed_values to allowed_values treatment combinations.

    Parameters
    ----------
    variables_dict_table : pandas dataframe
        The provided variables dictionary dataframe.

    Returns
    -------
    bool
        True: all values are valid.
        False: not all values are valid.

    """  
    #Init results for check combinations are valid (bools)
    result = []
    
    #Iter variables dict
    for i,r in variables_dict_table.iterrows():
        #Get temp values
        temp_treatment = r['allowed_values_treatment']
        temp_allowed_values = r['allowed_values']
        
        #Get temp value type
        temp_value_type = define_allowed_values_column_type(temp_allowed_values)
        
        #Set value based on combination
        if temp_treatment == 'any' and temp_value_type == 'any':
            result.append(True)
        elif temp_treatment == 'wholeset' and (temp_value_type == 'filled_list' or temp_value_type == 'filled_minmax_numeric_list'):
            result.append(True)
        elif temp_treatment == 'subset' and (temp_value_type == 'filled_list' or temp_value_type == 'filled_minmax_numeric_list'):
            result.append(True)
        elif temp_treatment == 'range' and temp_value_type == 'filled_minmax_numeric_list':
            result.append(True)
        else:
            result.append(False)
                        
    #Check if there are False values and return bool
    return all(result)


def define_variables_for_uniqueness_between_type(variables_dict_table, variables_cell, variable):
    """
    This function checks the provided string from 'variables_for_uniqueness_between' column and
    defines the type of value provided.

    Parameters
    ----------
    variables_dict_table : pandas dataframe
        The provided variables dictionary dataframe.
    variables_cell : str
        The provided temp cell value from variables_for_uniqueness_between column.
    variable : str
        The current associated variable from 'variable' column

    Returns
    -------
    result : str
        none -> none string indicating that no variable is provided
        NA_value -> the cell is empty
        not_list -> the cell is filled but not a Python list
        empty_list -> the cell is filled with an empty Python list
        invalid_list -> the cell is filled but does not meet criteria
        valid_list -> the cell is filled fullfills criteria(provided variables are in 'variable' column and is not the current variable been evaluated)

    """
    #Check the temp cell value
    if variables_cell == 'none':
        result = 'none'
    elif pd.isna(variables_cell):
        result = 'NA_value'
    else:
        #Try to convert
        ##Try to evaluate string to convert it to list
        try:
            evaluated_values = literal_eval(curly_2_straight_quotation(variables_cell))
        except:
            result = 'not_list'
        else:
            #Check if is list
            if not (type(evaluated_values) == list):
                result = 'not_list'
            else:
                #Check list is not empty
                if len(evaluated_values) == 0:
                    result = 'empty_list'
                else:
                    #Check if evaluated_values from temp_cell variables are all present in 'variable' column from dict
                    if not check_values_inv(variables_dict_table['variable'], evaluated_values):
                        result = 'invalid_list'
                    else:
                        #Check that current variable is not in evaluated_values from temp_cell 
                        if variable in evaluated_values:
                            result = 'invalid_list'
                        else:
                            result = 'valid_list'
    return result


def check_variables_for_uniqueness_between_column(variables_dict_table):
    """
    This function checks if all cells in variables_for_uniqueness_between column from
    Variables Dictionary file are valid values.

    Parameters
    ----------
    variables_dict_table : pandas dataframe
        The provided variables dictionary dataframe.

    Returns
    -------
    bool
        True: all values are valid.
        False: not all values are valid.

    """    
    #Init check_results values(bools)
    results = []
    
    for i,r in variables_dict_table.iterrows():
        #Get temp values
        temp_variable = r['variable']
        temp_cell = r['variables_for_uniqueness_between']
        
        #Define cell type
        temp_cell_type = define_variables_for_uniqueness_between_type(variables_dict_table, temp_cell, temp_variable)
        
        #Set value based on cell_type
        if temp_cell_type == 'none' or temp_cell_type =='valid_list':
            results.append(True)
        else:
            results.append(False)
    
    #Check if there are False values and return bool
    return all(results)


def check_variables_for_uniqueness_between_combs(variables_dict_table):
    """
    This function checks if all cells in 'variables_for_uniqueness_between' column to 
    'check_uniqueness_between_variables' combinations from Variables Dictionary file 
    are valid.

    Parameters
    ----------
    variables_dict_table : pandas dataframe
        The provided variables dictionary dataframe.

    Returns
    -------
    bool
        True: all values are valid.
        False: not all values are valid.

    """    
    #Init check_results values(bools)
    results = []
    
    for i,r in variables_dict_table.iterrows():
        #Get temp values
        temp_variable = r['variable']
        temp_cell = r['variables_for_uniqueness_between']
        temp_check_option = r['check_uniqueness_between_variables']
        
        #Define cell type
        temp_cell_type = define_variables_for_uniqueness_between_type(variables_dict_table, temp_cell, temp_variable)
        
        #Set reslut based on combination
        if temp_cell_type == 'none' and temp_check_option == 'no':
            results.append(True)
        elif temp_cell_type == 'valid_list' and temp_check_option == 'yes':
            results.append(True)
        else:
            results.append(False)
    
    #Check if there are False values and return bool
    return all(results)


def convert_str_values(value):
    """
    This functions converts the provided str value to a valid object that 
    we can work with. In the previous checks we have already ensure that only
    'any', 'none' or valid filled python list will be returned.

    Parameters
    ----------
    value : str
        Provided str to convert.

    Returns
    -------
    result : str or list
        'any', 'none' or valid filled python list.

    """
    if value == 'any':
        result = 'any'
    elif value == 'none':
        result = 'none'
    elif pd.isna(value):
        result = 'NA_value'
    else:
        #Try to convert
        ##Try to evaluate string to convert it to list
        try:
            evaluated_values = literal_eval(curly_2_straight_quotation(value))
        except:
            result = 'not_list'
        else:
            #Check if is list
            if not (type(evaluated_values) == list):
                result = 'not_list'
            else:
                #Check if list is not empty
                if len(evaluated_values) == 0:
                    result = 'empty_list'
                else:
                    result = evaluated_values
    return result


def check_class_variable(variable, class_type, metadata_table, color_treatment):
    """
    This function handles the class_type check for the provided variable.

    Parameters
    ----------
    variable : str
        Variable (column) to check in the metadata table.
    class_type : str
        The type of data expected in this variables[character or numeric].
    metadata_table : pandas dataframe
        Metadata table.
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
    #Check class_type
    if class_type == 'character':
        if metadata_table[variable].dtype == 'O' or metadata_table[variable].dtype == 'bool':
            print(rich_text_colored('Okey! Indicated as character and dtype is Object or Bool!','acceptable', color_treatment))
            return False
        else:
            print(rich_text_colored('Warning! Indicated as character but dtype is not Object or Bool!', 'program_warning2', color_treatment))
            return True
    elif class_type == 'numeric':
        if metadata_table[variable].dtype == 'int64' or metadata_table[variable].dtype == 'float64':
            print(rich_text_colored('Okey! Indicated as numeric and dtype is numeric!', 'acceptable', color_treatment))
            return False
        else:
            print(rich_text_colored('Warning! Indicated as numeric and dtype is not numeric!', 'program_warning2', color_treatment))
            return True
    else:
        print(rich_text_colored('Error! Permited options for class_type column are "character" or "numeric"!', 'dangerous', color_treatment))
        return True


def check_uniqueness_within_variable(variable, uniqueness_value, metadata_table, color_treatment):
    """
    This function handles the uniqueness within check for the provided variable. 

    Parameters
    ----------
    variable : str
        Variable (column) to check in the metadata table.
    uniqueness_value : str
        The type of uniqueness expected in this variable[unique or Non-unique].
    metadata_table : pandas dataframe
        Metadata table.
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
    #Get number of duplicates
    n_dup = len(metadata_table[metadata_table.duplicated(variable)])
    
    #Check uniqueness
    if uniqueness_value == 'unique':
        if n_dup == 0:
            print(rich_text_colored('Okey! Indicated as unique and no duplicates were found within variable!','acceptable', color_treatment))
            return False
        else:
            print(rich_text_colored('Warning! Indicated as unique but duplicates were found!','program_warning2', color_treatment))
            print(rich_text_colored('The following duplicates were found:','general_text', color_treatment))
            duplicates = list(set(metadata_table[metadata_table.duplicated(variable)][variable]))
            print_list_n_byline(duplicates, 5)
            return True
    elif uniqueness_value == 'nonunique':
        if n_dup == 0:
            print(rich_text_colored('Okey! Indicated as nonunique! No duplicates were found within variable!','acceptable', color_treatment))
        else:
            print(rich_text_colored('Okey! Indicated as nonunique! Duplicates were found within variable!','acceptable', color_treatment))
        return False
    else:
        print(rich_text_colored('Error! Permited options for "uniqueness_within_variable" column are "unique" or "nonunique"!','dangerous', color_treatment))
        return True


def check_uniqueness_by_variable(main_variable, check_variable, check_uniqueness_between, metadata_table, color_treatment):
    """
    This function checks if the are multiple matches for each value of the main_variable in the procided check_variable. 
    1) Duplicates will be searched in the main_variable
    2) If duplicates are found the function will find if some of these duplicates have multiple matches in the provided check_variable.

    Parameters
    ----------
    main_variable : str
        Current variable of the variables dictionary. It will be used as reference.
    check_variable : str
        Provided variable to check for multiple matches between variables.
    check_uniqueness_between : str
        Option provided in the "check_uniqueness_between_variables" column from variables dictionary [yes, or no].
    metadata_table : pandas dataframe
        Metadata table.
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
    if check_uniqueness_between == 'no':
        print(rich_text_colored(''.join(['Okey! Indicated as check no! Skipping check uniqueness between variables for ',main_variable, '!']),'acceptable', color_treatment))
        return False
    elif check_uniqueness_between == 'yes':
        if check_variable not in metadata_table.columns:
            #If not present skip checks and lunch warning
            print(rich_text_colored(''.join(['Warning! The provided variable ',check_variable, ' is not present in the Curated Metadata Table! Skipping check uniqueness between variable ',main_variable, ' and the provided variable ', check_variable,'!']), 'program_warning2', color_treatment))
            return True
        else:
            #get list of duplicates for main_variable and then for each dup see if check_variable has a unique match
            duplicates = list(pd.unique(metadata_table[metadata_table.duplicated(main_variable)][main_variable]))
            if len(duplicates)==0:
                print(rich_text_colored(''.join(['Okey! No duplicates were detected in variable ',main_variable,' thus each value has a unique match with the provided variable ',check_variable,'!']),'acceptable', color_treatment))
                return False
            else:
                #Set results
                duplicates_by_2ndvariable = []
                #Check each dup by check_variable
                for dup in duplicates:
                    #Get rows for this value
                    if pd.isna(dup):
                        temp_dup_rows = metadata_table[metadata_table[main_variable].isna()]
                    else:
                        temp_dup_rows = metadata_table[metadata_table[main_variable]==dup]
                    #Count number of unique matches
                    temp_matches = set(temp_dup_rows[check_variable])
                    #Add to results if there is more than one match
                    if len(temp_matches) > 1:
                        duplicates_by_2ndvariable.append(dup)
                #Show message result
                if len (duplicates_by_2ndvariable) == 0:
                    print(rich_text_colored(''.join(['Okey! Duplicates were detected in variable ',main_variable,' and each value has a unique match with the provided variable ',check_variable,'!']),'acceptable', color_treatment))
                    return False
                else:
                    print(rich_text_colored(''.join(['Warning! Duplicates were detected in variable ',main_variable,' but these have more than one match with the provided variable ',check_variable,'!']),'program_warning2', color_treatment))
                    print(rich_text_colored('The following duplicates with multiple matches were found:','general_text', color_treatment))
                    print_list_n_byline(duplicates_by_2ndvariable, 5)
                    return True 
    else:
        print(rich_text_colored('Error! Permited options for "check_uniqueness_between_variables" column are "yes" or "no"!','dangerous', color_treatment))
        return True
        
    
def check_allowed_values(variable, allowed_values_treatment, allowed_values_cell, NAs_allowed, metadata_table, color_treatment):
    """
    This function checks the allowed values for the provided variable in the metadata table
    based on the information provided in the columns "NAs_allowed", "allowed_values_treatment" 
    and "allowed_values" of the variables dictionary.

    Parameters
    ----------
    variable : str
        Variable (Column name) in metadata table that is going to be checked.
    allowed_values_treatment : str
        Type of treatment to check['any', 'range', 'subset', 'wholeset'].
        The values provided in "allowed_values_treatment" column of the variables dictionary.
    allowed_values_cell : str or list
        Allowed values provided for the variable.
        The allowed values provided in "allowed_values" column of the variables dictionary.
        After conversion with convert_str_values() function.
    NAs_allowed : str
        Type of NAs treatment['yes', 'no'].
        The values provided in "NAs_allowed" column of the variables dictionary.
    metadata_table : pandas dataframe
        Metadata table.
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
    if allowed_values_treatment == 'any':
        if NAs_allowed == 'yes':
            print(rich_text_colored('Okey! Allowed values treatment "any" and NAs are allowed! Skipping check!','acceptable', color_treatment))
            return False
        else:
            #Check if there is any NA
            check_NAs = metadata_table[variable].isna().any()
            if check_NAs == True:
                print(rich_text_colored('Warning! NAs were detected, but NAs were not allowed!','program_warning2', color_treatment))
                return True
            else:
                print(rich_text_colored('Okey! NAs were not allowed and no NAs were detected!','acceptable', color_treatment))
                return False
    elif allowed_values_treatment == 'range':
        #Get min and max values
        minim = allowed_values_cell[0]
        maxim = allowed_values_cell[1]
        #Treat NAs_allowed
        if NAs_allowed == 'yes':
            #Check if values are in range or NAs
            range_or_na = (((metadata_table[variable] >= minim) & (metadata_table[variable] <= maxim)) | metadata_table[variable].isna())
            #Treat result
            if range_or_na.all() == True:
                print(rich_text_colored('Okey! All values are in the allowed numeric range or are NAs!','acceptable', color_treatment))
                return False
            else:
                print(rich_text_colored(''.join(['Warning! Some elements in variable "',variable,'" are not in the provided allowed numeric range and are not NAs either!']),'program_warning2', color_treatment))
                print(rich_text_colored('The following conflicting elements were found:','general_text', color_treatment))
                conflicting_vals = list(pd.unique(metadata_table[range_or_na==False][variable]))
                print_list_n_byline(conflicting_vals, 5)
                return True
        else:
            #Just check if elements are in range
            range_bools = ((metadata_table[variable] >= minim) & (metadata_table[variable] <= maxim))
            #Treat result
            if range_bools.all() == True:
                print(rich_text_colored('Okey! All values are in the allowed numeric range!','acceptable', color_treatment))
                return False
            else:
                print(rich_text_colored(''.join(['Warning! Some elements in variable "',variable,'" are not in the provided allowed numeric range!']),'program_warning2', color_treatment))
                print(rich_text_colored('The following conflicting elements were found:','general_text', color_treatment))
                conflicting_vals = list(pd.unique(metadata_table[range_bools==False][variable]))
                print_list_n_byline(conflicting_vals, 5)
                return True
    elif allowed_values_treatment == 'wholeset':
        #Check that all allowed values are in variable column in metadata
        present_check = check_values_inv(metadata_table[variable], allowed_values_cell)
        #Get unique values and get difference between variable column and allowed values list
        variable_mt_col_unique = set(pd.unique(metadata_table[variable]))
        allowed_values_cell_unique = set(pd.unique(allowed_values_cell))
        diff_variable_mt_col_unique = list(variable_mt_col_unique.difference(allowed_values_cell_unique))
        diff_allowed_values_cell_unique = list(allowed_values_cell_unique.difference(variable_mt_col_unique))
        #Treat NA independend cases
        ##All detected and no extra values in metadata column
        if present_check==True and len(diff_variable_mt_col_unique)==0:
            print(rich_text_colored(''.join(['Okey! All values in the provided allowed values are present and there are no extra values in the "',variable,'" column in the metadata table!']),'acceptable', color_treatment))
            return False
        ##Not all detected and no extra values in metadata column
        elif present_check==False and len(diff_variable_mt_col_unique)==0:
            print(rich_text_colored('Warning! Not all values in the provided allowed values are present!','program_warning2', color_treatment))
            print(rich_text_colored(''.join(['There are no extra values in the "',variable,'" column in the metadata table!']), 'acceptable', color_treatment))
            print(rich_text_colored('The following provided values were absent in the metadata table:','general_text', color_treatment))
            print_list_n_byline(diff_allowed_values_cell_unique, 5)
            return True
        ##All detected and extra values are only NAs
        elif present_check==True and len(diff_variable_mt_col_unique)==1 and pd.isna(diff_variable_mt_col_unique[0]):
            #Treat NAs_allowed
            if NAs_allowed == 'yes':
                print(rich_text_colored(''.join(['Okey! All values in the provided allowed values are present!\nNAs were allowed and detected extra values are only NAs in the "',variable,'" column in the metadata table!\nNo more extra values were detected!']),'acceptable', color_treatment))
                return False
            else:                
                print(rich_text_colored('Warning!', 'program_warning2', color_treatment),rich_text_colored('All values in the provided allowed values are present!', 'acceptable', color_treatment))
                print(rich_text_colored(''.join(['However NAs were not allowed and detected extra values are only NAs in the "',variable,'" column in the metadata table!']),'program_warning2', color_treatment))
                print(rich_text_colored('No more extra values were detected!', 'acceptable', color_treatment))
                return True
        ##Not all deteceted and extra values are only NAs
        elif present_check==False and len(diff_variable_mt_col_unique)==1 and pd.isna(diff_variable_mt_col_unique[0]):
            #Treat NAs_allowed
            if NAs_allowed == 'yes':
                print(rich_text_colored('Warning! Not all values in the provided allowed values are present!', 'program_warning2', color_treatment))
                print(rich_text_colored(''.join(['NAs were allowed and detected extra values are only NAs in the "',variable,'" column in the metadata table!\nNo more extra values were detected!']),'acceptable', color_treatment))
                print(rich_text_colored('The following provided values were absent in the metadata table:','general_text', color_treatment))
                print_list_n_byline(diff_allowed_values_cell_unique, 5)
                return True
            else:
                print(rich_text_colored(''.join(['Warning! Not all values in the provided allowed values are present!\nAdditionally NAs were not allowed and detected extra values are only NAs in the "',variable,'" column in the metadata table!']),'program_warning2', color_treatment))
                print(rich_text_colored('No more extra values were detected!', 'acceptable', color_treatment))
                print(rich_text_colored('The following provided values were absent in the metadata table:','general_text', color_treatment))
                print_list_n_byline(diff_allowed_values_cell_unique, 5)
                return True
        #Every other case
        else:
            #If NAs are allowed the only thing that we do is filter the value from diff_variable_mt_col_unique
            if NAs_allowed == 'yes':
                #Filter nan from new list
                diff_variable_mt_col_unique_mod = [x for x in diff_variable_mt_col_unique if pd.isna(x)==False]
            else:
                diff_variable_mt_col_unique_mod = diff_variable_mt_col_unique
    
            #Show different messages
            if present_check==True:
                print(rich_text_colored('Warning!', 'program_warning2', color_treatment),rich_text_colored('All values in the provided allowed values are present!', 'acceptable', color_treatment))
                print(rich_text_colored(''.join(['However extra values were detected in the "',variable,'" column in the metadata table!']),'program_warning2', color_treatment))
                print(rich_text_colored('The following extra values were present in the metadata table:','general_text', color_treatment))
                print_list_n_byline(diff_variable_mt_col_unique_mod, 5)
                return True
            else:
                print(rich_text_colored(''.join(['Warning! Not all values in the provided allowed values are present!\nAdditionally extra values were detected in the "',variable,'" column in the metadata table!']),'program_warning2', color_treatment))
                print(rich_text_colored('The following provided values were absent in the metadata table:','general_text', color_treatment))
                print_list_n_byline(diff_allowed_values_cell_unique, 5)
                print(rich_text_colored('The following extra values were present in the metadata table:','general_text', color_treatment))
                print_list_n_byline(diff_variable_mt_col_unique_mod, 5)
                return True
    elif allowed_values_treatment == 'subset':
        #Treat NAs_allowed
        if NAs_allowed == 'yes':
            #Check if values in variable column are present in allowed values or NAs
            subset_or_na = (metadata_table[variable].isin(allowed_values_cell) | metadata_table[variable].isna())
            #Treat result
            if subset_or_na.all() == True:
                print(rich_text_colored('Okey! All values are allowed values or are NAs!','acceptable', color_treatment))
                return False
            else:
                print(rich_text_colored(''.join(['Warning! Some elements in variable "',variable,'" are not in the provided allowed values and are not NAs either!']),'program_warning2', color_treatment))
                print(rich_text_colored('The following conflicting elements were found:','general_text', color_treatment))
                conflicting_vals = list(pd.unique(metadata_table[subset_or_na==False][variable]))
                print_list_n_byline(conflicting_vals, 5)
                return True
        else:
            #Just check if elements are in allowed_values
            subset_bools = metadata_table[variable].isin(allowed_values_cell)
            #Treat result
            if subset_bools.all() == True:
                print(rich_text_colored('Okey! All values are allowed values!','acceptable', color_treatment))
                return False
            else:
                print(rich_text_colored(''.join(['Warning! Some elements in variable "',variable,'" are not in the provided allowed values!']),'program_warning2', color_treatment))
                print(rich_text_colored('The following conflicting elements were found:','general_text', color_treatment))
                conflicting_vals = list(pd.unique(metadata_table[subset_bools==False][variable]))
                print_list_n_byline(conflicting_vals, 5)
                return True
    else:
        print(rich_text_colored('Error! Permited options for "allowed_values_treatment" column are "any", "subset", "wholeset" or "range"!','dangerous', color_treatment))
        return True
                        

def check_requiredness_by_variable(variable, requiredness, metadata_table, color_treatment):
    """
    This function checks the presence of the provided variable in the metadata_table
    and gets the messages of the requiredness check by variable.
    
    Note:I do a previus initial comprobation with check_required_variables_in_metadata_table()

    Parameters
    ----------
    variable : str
        Variable (Column name) in metadata table that is going to be checked.
    requiredness : str
        Type of requiredness ['required', 'optional'].
        The values provided in "requiredness" column of the variables dictionary.
    metadata_table : pandas dataframe
        Metadata table.
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
    if requiredness == 'required':
        if variable not in metadata_table.columns:
            #If not present skip all remaining checks and lunch warning
            print(rich_text_colored('Warning! This required variable is not in the provided Curated Metadata Table! Skipping all remaining checks!', 'program_warning2', color_treatment))
            return True
        else:
            #If present continue with chesks
            print(rich_text_colored('Okey! This required variable is in the provided Curated Metadata Table!', 'acceptable', color_treatment))
            return False
    elif requiredness == 'optional':
        if variable not in metadata_table.columns:
            #If not present skip all remaining checks and lunch warning
            print(rich_text_colored('Warning! This optional variable is not in the provided Curated Metadata Table! Skipping all remaining checks!', 'program_warning2', color_treatment))
            return True
        else:
            #If present continue with chesks
            print(rich_text_colored('Okey! This optional variable is in the provided Curated Metadata Table!', 'acceptable', color_treatment))
            return False
    else:
        print(rich_text_colored('Error! Permited options for "requiredness" column are "required", or "optional"!','dangerous', color_treatment))
        return True

    
#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter metadata_table
    parser.add_argument(
            '-t','--metadata_table',
            action = 'store',
            required = True,
            help ='Curated Metadata Table [Expected sep=TABS]. Indicate the path to the Curated Metadata Table file.'
    )
    ##Parameter variables_dictionary
    parser.add_argument(
            '-d','--variables_dictionary',
            action = 'store',
            required = True,
            help ='Variables Dictionary [Expected sep=TABS]. Indicate path to the Variables Dictionary file. See Documentation for more information and format details.'
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
    metadata_path = args.metadata_table
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
        
        #1) Load files and do checks
        #Section header message
        print(rich_text_colored('\nLoading Files:', 'section_header', plain_text_bool))
        
        #1.1.Try to load Variables Dictionary File as pandas dataframe and do checks
        ##Show loading file message
        print(rich_text_colored('\nVariables Dictionary file:', 'general_text', plain_text_bool))
        print(variables_dictionary_path)
        
        ##A)Try to load the variables dict
        variables_dict_table = pd.read_csv(variables_dictionary_path, sep = '\t')
        
        ##B)Do initial checks on Variables Dictionary file for columns of interest
        
        ##Check headers
        check_headers_variablesfile = check_headers(VARIABLES_FILE_HEADERS, variables_dict_table)
        frase0_1 = 'Error! Some of the needed headers are not in the Variables Dictionary!\n Check your variables file!\n'
        frase1_1 = '\nThe headers needed for this mode are:'
        treat_headers_check(VARIABLES_FILE_HEADERS, check_headers_variablesfile, variables_dict_table, frase0_1, frase1_1, plain_text_bool)
        
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
        
        ##Check valid 'class_type' column values
        check_class_type_values = check_values(VALID_CLASS_TYPE, variables_dict_table['class_type'])
        frase0_4 = 'Error! Some of the values of the "class_type" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_4 = '\nValid values are:'
        treat_values_check(VALID_CLASS_TYPE, check_class_type_values, frase0_4, frase1_4, plain_text_bool)
        
        ##Check valid 'uniqueness_within_variable' column values
        check_uniqueness_values = check_values(VALID_UNIQUENESS, variables_dict_table['uniqueness_within_variable'])
        frase0_5 = 'Error! Some of the values of the "uniqueness_within_variable" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_5 = '\nValid values are:'
        treat_values_check(VALID_UNIQUENESS, check_uniqueness_values, frase0_5, frase1_5, plain_text_bool)
        
        ##Check valid 'check_uniqueness_between_variables' column values
        check_uniqueness_between_variables_options_values = check_values(VALID_CHECK_UNIQUENESS_BETWEEN_OPTIONS, variables_dict_table['check_uniqueness_between_variables'])
        frase0_6 = 'Error! Some of the values of the "check_uniqueness_between_variables" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_6 = '\nValid values are:'
        treat_values_check(VALID_CHECK_UNIQUENESS_BETWEEN_OPTIONS, check_uniqueness_between_variables_options_values, frase0_6, frase1_6, plain_text_bool)
        
        ##Check valid 'variables_for_uniqueness_between' column values
        check_variables_4_uniqueness_between_values = check_variables_for_uniqueness_between_column(variables_dict_table)        
        frase0_7 = 'Error! Some of the values of the "variables_for_uniqueness_between" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_7 = '\nValid values are:'
        treat_values_check(VALID_UNIQUENESS_BETWEEN_VARIABLES, check_variables_4_uniqueness_between_values, frase0_7, frase1_7, plain_text_bool)
        
        ##Check valid 'check_uniqueness_between_variables' to 'variables_for_uniqueness_between' column combinations
        check_vars_for_uniqueness_between_comb = check_variables_for_uniqueness_between_combs(variables_dict_table)
        frase0_8 = 'Error! Some of the values of the "allowed_values_treatment" to "allowed_values" combinations in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_8 = '\nValid combinations are:'
        treat_values_check(VALID_UNIQUENESS_BETWEEN_VARIABLES_COMB, check_vars_for_uniqueness_between_comb, frase0_8, frase1_8, plain_text_bool)
        
        ##Check valid 'NAs_allowed' column values
        check_NAs_allowed_values = check_values(VALID_NAS_ALLOWED, variables_dict_table['NAs_allowed'])
        frase0_9 = 'Error! Some of the values of the "NAs_allowed" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_9 = '\nValid values are:'
        treat_values_check(VALID_NAS_ALLOWED, check_NAs_allowed_values, frase0_9, frase1_9, plain_text_bool)
        
        ##Check valid 'allowed_values_treatment' column values
        check_allowed_values_treatment_values = check_values(VALID_ALLOWED_VALUES_TREATMENT, variables_dict_table['allowed_values_treatment'])
        frase0_10 = 'Error! Some of the values of the "allowed_values_treatment" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_10 = '\nValid values are:'
        treat_values_check(VALID_ALLOWED_VALUES_TREATMENT, check_allowed_values_treatment_values, frase0_10, frase1_10, plain_text_bool)
        
        ##Check that "class_type" and "allowed_values_treatment" column combinations are valid
        #Get "class_type" and "allowed_values_treatment" column combinations in variables dictionary
        type_sum_treatment = variables_dict_table['class_type'] + ' + ' + variables_dict_table['allowed_values_treatment']
        #Check that combinations are valid
        check_expected_type_treatment_sum_values = check_values(VALID_CLASSTYPE_TREATMENT_COMB, type_sum_treatment)
        frase0_11 = 'Error! Some of the "class_type" to "allowed_values_treatment" combinations in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_11 = '\nValid combinations are:'
        treat_values_check(VALID_CLASSTYPE_TREATMENT_COMB, check_expected_type_treatment_sum_values, frase0_11, frase1_11, plain_text_bool)
        
        ##Check that "requiredness" and "NAs_allowed" column combinations are valid
        #Get "requiredness" and "NAs_allowed" column combinations in variables dictionary
        requiredness_sum_NAs_allowed = variables_dict_table['requiredness'] + ' + NAs_allowed:' + variables_dict_table['NAs_allowed']
        #Check that combinations are valid
        check_expected_requiredness_NAsallowed_sum_values = check_values(VALID_REQUIREDNESS_NAS_COMB, requiredness_sum_NAs_allowed)
        frase0_12 = 'Error! Some of the "requiredness" to "NAs_allowed" combinations in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_12 = '\nValid combinations are:'
        treat_values_check(VALID_REQUIREDNESS_NAS_COMB, check_expected_requiredness_NAsallowed_sum_values, frase0_12, frase1_12, plain_text_bool)
        
        ##Check valid 'allowed_values' column values
        init_check_allowed_values = check_allowed_values_column(variables_dict_table)
        frase0_13 = 'Error! Some of the values of the "allowed_values" column in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_13 = '\nValid values are:'
        treat_values_check(VALID_ALLOWED_VALUES, init_check_allowed_values, frase0_13, frase1_13, plain_text_bool)
        
        ##Check valid 'allowed_values_treatment' on 'allowed_values' combinations
        check_treatment_allowed_vals_comb = check_allowed_values_treatment_comb(variables_dict_table)
        frase0_14 = 'Error! Some of the values of the "allowed_values_treatment" to "allowed_values" combinations in the Variables Dictionary are not valid!\n Check your variables file!\n'
        frase1_14 = '\nValid combinations are:'
        treat_values_check(VALID_TREATMENT_ALOWED_VALUES_COMB, check_treatment_allowed_vals_comb, frase0_14, frase1_14, plain_text_bool)
        
        #1.2.Try to load Curated Metadata Table File as pandas dataframe and do checks
        ##Show loading file message
        print(rich_text_colored('\nCurated Metadata Table file:', 'general_text', plain_text_bool))
        print(metadata_path)
        
        ##A)Load Curated Metadata File
        mt_table = pd.read_csv(metadata_path, sep = '\t')
        
        ##B) Do initial checks for provided metadata table
        
        #Initial preparation
        ##Get list of required columns
        required_cols = variables_dict_table[variables_dict_table['requiredness']=='required']['variable']
        ##Get list of all possible columns
        universe_cols = list(variables_dict_table['variable'])
        
        ##Check required variables
        check_required_variables_in_metadata_table(mt_table, required_cols, plain_text_bool)
        
        ##Check that all columns in metadata are in the variables dictionary
        check_metadata_table_vars_in_dict(mt_table, universe_cols, plain_text_bool)
        
        #2. Do main checks (selected analyses)
        #Section header message
        print(rich_text_colored('\nChecking Provided Variables:', 'section_header', plain_text_bool))
        
        #Init warnings detected list
        warnings_list = []
        
        #For each variable do analysis or selected main checks
        for i,r in variables_dict_table.iterrows():
            #Get temp check info
            temp_var = r['variable']
            temp_req = r['requiredness']
            temp_class = r['class_type']
            temp_uniqueness_within = r['uniqueness_within_variable']
            temp_uniqueness_between_check = r['check_uniqueness_between_variables']
            temp_uniqueness_between_vars = r['variables_for_uniqueness_between']
            temp_allowed_vals_NAs = r['NAs_allowed']
            temp_allowed_vals_treatment = r['allowed_values_treatment']
            temp_allowed_vals = r['allowed_values']
            
            ##Message - Show variable information
            print(rich_text_colored('\nVariable:', 'variable_color', plain_text_bool), temp_var)
            print('1)',rich_text_colored('Requiredness:', 'dict_column_color', plain_text_bool), temp_req)
            #Check if variable is required and if it is in metadata table
            temp_check_requiredness_by_variable = check_requiredness_by_variable(temp_var, temp_req, mt_table, plain_text_bool)
            #Add check result
            warnings_list.append(temp_check_requiredness_by_variable)
            #Check that temp variable is in metadata table(if there are no warnings continue with the rest of the checks)
            if temp_check_requiredness_by_variable == False:
                
                #Do checks for temp varible
                ##Check class_type
                print('2)',rich_text_colored('Class Type:', 'dict_column_color', plain_text_bool), temp_class)
                temp_check_class = check_class_variable(temp_var, temp_class, mt_table, plain_text_bool)
                warnings_list.append(temp_check_class)
                
                ##Check uniqueness_within_variable
                print('3)',rich_text_colored('Uniqueness Within:', 'dict_column_color', plain_text_bool), temp_uniqueness_within)
                temp_check_uniq_within = check_uniqueness_within_variable(temp_var, temp_uniqueness_within, mt_table, plain_text_bool)
                warnings_list.append(temp_check_uniq_within)
                
                ##Check uniqueness_between provided variables
                print('4)',rich_text_colored('Uniqueness Between:', 'dict_column_color', plain_text_bool))
                print(rich_text_colored('Check Uniqueness Between variables:', 'dict_column_color2', plain_text_bool), temp_uniqueness_between_check)
                #Convert cell values
                temp_vars_cell = convert_str_values(temp_uniqueness_between_vars)
                #Treat different values
                #Previously we have already check that all values are permited with other functions
                #Expected values are the string 'none' or a python list with colnames present in the metadata table
                if temp_vars_cell == 'none':
                    print('-',rich_text_colored('Variable to check:', 'dict_column_color2', plain_text_bool), temp_vars_cell)
                    temp_check_uniq_between = check_uniqueness_by_variable(temp_var, None, temp_uniqueness_between_check, mt_table, plain_text_bool)
                    warnings_list.append(temp_check_uniq_between)
                else:
                    for check_var in temp_vars_cell:
                        print('-',rich_text_colored('Variable to check:', 'dict_column_color2', plain_text_bool), check_var)
                        temp_check_uniq_between = check_uniqueness_by_variable(temp_var, check_var, temp_uniqueness_between_check, mt_table, plain_text_bool)
                        warnings_list.append(temp_check_uniq_between)
                        
                ##Check allowed values using "NAs_allowed"	"allowed_values_treatment"	"allowed_values"
                print('5)',rich_text_colored('Allowed Values:', 'dict_column_color', plain_text_bool))
                print('-',rich_text_colored('Allowed Values Treatment:', 'dict_column_color2', plain_text_bool), temp_allowed_vals_treatment)
                print('-',rich_text_colored('NAs Allowed:', 'dict_column_color2', plain_text_bool), temp_allowed_vals_NAs)
                #Convert cell values
                temp_allowed_values_cell = convert_str_values(temp_allowed_vals)
                #Show modified message for range treatment
                if temp_allowed_vals_treatment == 'range':
                    print('-',rich_text_colored('Values Allowed[min,max]:', 'dict_column_color2', plain_text_bool), temp_allowed_values_cell)
                else:
                    print('-',rich_text_colored('Values Allowed:', 'dict_column_color2', plain_text_bool), temp_allowed_values_cell)
                #Treat different values
                #Previously we have already check that all values are permited with other functions
                #Expected values are the string 'any' or a python list with allowed values or ranges
                temp_check_allowed_values = check_allowed_values(temp_var, temp_allowed_vals_treatment, temp_allowed_values_cell, temp_allowed_vals_NAs, mt_table, plain_text_bool)
                warnings_list.append(temp_check_allowed_values)
    
        #Add leyend and phrases since it is a check program
        if any(warnings_list) == True:
            #Show advise message
            print(rich_text_colored('\nThis could be due to:', 'due_to_header', plain_text_bool))
            print('- The absence of optional variables in the Curated Metadata Table', rich_text_colored('[Acceptable]', 'acceptable', plain_text_bool))
            print('- Mismatches between the observed dtypes for variables in the Curated Metadata Table and the provided class types', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- The presence of nonunique variables indicated as unique', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- The absence of variables in the Curated Metadata Table that were provided to perform a uniqueness between variables check', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- The presence of multiple matches when performing a uniqueness between variables check', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print('- The presence of incompatible values when performing allowed values check', rich_text_colored('[Warning]', 'legend_warning', plain_text_bool))
            print(rich_text_colored('\nYou should:','you_should_header', plain_text_bool))
            print('- Manually confirm which is your case')
            #Show leyend
            show_advise_legend(plain_text_bool)
    
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
    
    finally:
        #Print empty line for aesthetic purposes
        print('')    

if __name__ == '__main__':
    main()
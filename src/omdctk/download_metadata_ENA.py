#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: download_metadata_ENA.py
@author: Samuel Piquer-Esteban
@date: 7 Jul 2023 

"""

#Set program name
__program__ = 'download_metadata_ENA.py'

#Imports from MTD_CT common module
from omdctk import (DATE, VERSION, OMD_CTK_Exception, program_header, 
                    treat_output_directory_parameter_outfiles,
                    check_existence_directory_parameter, rich_text_colored)

#Import third-party modules
from argparse import ArgumentParser
from tabulate import tabulate
from urllib.request import (urlretrieve, Request, urlopen)
import pandas as pd
import mg_toolkit

#Program functions
def download_metadata_ENA_browser(study_accession, outputdir_path, color_treatment):
    """
    This function downloads the ENA Browser metadata file for the
    provided study_accession on the outputdir_path.

    Parameters
    ----------
    study_accession : str
        ENA project study accession.
    outputdir_path : str
        Provided output directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If the constructed URL is empty (response code 204) raises an exception, 
        since there will be no information to work with. The most likely 
        explanation is that the provided Study Accession is not in ENA. 

    Returns
    -------
    outputfile_tech : str
        Full path to the downloaded file.

    """      
    #Previous steps
    ##Construct url
    
    #url ENA Api Previous
    #url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + study_accession + '&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&format=tsv&download=true&limit=0'
    
    #url ENA Api: 2023-05-02 (they have remove some columns and modified some columns)
    #They have eliminated cram columns (cram_index_ftp,cram_index_aspera,cram_index_galaxy)
    #They have modified read_count(now uses USA annotation with commas for thousands) and base_count(no longer the int exac number, but a value in MB)
    url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + study_accession + '&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&format=tsv&download=true&limit=0'
    
    ##Get outfile name
    out_name_tech = study_accession + '_ENA_browser.tsv'
    ##Treat output_directory parameter / Get full output file path
    outputfile_tech = treat_output_directory_parameter_outfiles(out_name_tech, outputdir_path)
    
    #Show downloading message
    print(rich_text_colored('Downloading:', 'general_text', color_treatment))
    print(url)
    
    #Check if the provided url is empty (probably Accession not found in ENA) / Else continue
    response = urlopen(Request(url))
    if response.code == 204:
        raise OMD_CTK_Exception('Error! ULR has No Content. The most likely explanation is that the provided Study Accession is not in ENA!\n Check the provided --project parameter!')
    else:
        #Download file from ENA
        urlretrieve(url, outputfile_tech)
        #Show saved file message
        print(rich_text_colored('\nSaved file:', 'general_text', color_treatment))
        print(outputfile_tech)
        #Return full path to downloaded file
        return outputfile_tech


def download_metadata_mgtoolkit(study_accession, outputdir_path, color_treatment):
    """
    This function retrieves the mg-toolkit metadata for the
    provided study_accession and creates the mg-toolkit 
    metadata file on the provided outputdir_path.

    Parameters
    ----------
    study_accession : str
        ENA project study accession.
    outputdir_path : str
        Provided output directory.
    color_treatment: bool
        The color treatment option provided.
        True : Plain Text
        False : Colored Text

    Raises
    ------
    OMD_CTK_Exception
        If the retrieved pandas dataframe is empty raises an exception 
        since there will be no information to work with. The most likely 
        explanation is that the provided Study Accession is not in ENA.

    Returns
    -------
    meta_table_t : pandas dataframe
        Samples metadata retrieved from ENA with mg-toolkit.

    """   
    #Show downloading file message
    print(rich_text_colored('Downloading with:', 'general_text', color_treatment))
    #Get mg-toolkit version
    print('mg-toolkit')
    print(mg_toolkit.__version__)
    print('')
    
    #Previous steps
    ##Get outfile name
    out_name_meta = study_accession + '_mg-toolkit.tsv'
    ##Treat output_directory parameter / Get full output file path
    outputfile_meta = treat_output_directory_parameter_outfiles(out_name_meta, outputdir_path)
    
    #Fetch samples metadata with mg_toolkit
    mg_project = mg_toolkit.metadata.OriginalMetadata(study_accession)
    meta_dict = mg_project.fetch_metadata()
    meta_table = pd.DataFrame.from_dict(meta_dict)
    meta_table_t = meta_table.transpose()
    meta_table_t.reset_index(inplace = True)
    meta_table_t = meta_table_t.rename(columns = {'index':'Run'})
    
    #Check if the generated table is empty (Accession not found in ENA) / Else continue
    if meta_table_t.empty:
        raise OMD_CTK_Exception('Error! mg-toolkit returned an empty table. The most likely explanation is that the provided Study Accession is not in ENA!\n Check the provided --project parameter!')
    else:
        #Show saved file message
        print(rich_text_colored('Saving in file:', 'general_text', color_treatment))
        print(outputfile_meta)
        #Save outputfile_meta
        meta_table_t.to_csv(outputfile_meta, header = True, index = False, sep = '\t')
        #Return the pandas dataframe
        return meta_table_t


#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter project
    parser.add_argument(
            '-p','--project', 
            action = 'store',
            required = True,
            help = 'Project Accession. Indicate ENA Study Accession.'
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
    study_accession = args.project
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
        #0)Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        #1)Try to Download and Load ENA technical metadata (ENA Browser)
        
        #Section header message
        print(rich_text_colored('\nENA Technical Metadata (ENA Browser):\n', 'section_header', plain_text_bool))
        
        #Download ENA Browser metadata
        outputfile_tech = download_metadata_ENA_browser(study_accession, outputdir_path, plain_text_bool)
        
        #Load file as pandas dataframe
        ##Show loading file message
        print(rich_text_colored('\nLoading file:', 'general_text', plain_text_bool))
        print(outputfile_tech)
        ##Try to load the downloaded file
        tech_table = pd.read_csv(outputfile_tech, sep = '\t',thousands=',')
        
        #2)Try to Download ENA samples metadata (mg-toolkit)
        
        #Section header message
        print(rich_text_colored('\nENA Samples Metadata (mg-toolkit):\n', 'section_header', plain_text_bool))
        
        #Download mg-toolkit metadata
        mgtoolkit_table = download_metadata_mgtoolkit(study_accession, outputdir_path, plain_text_bool)
        
        #3)Try to merge both tables
        
        #Section header and messages
        print(rich_text_colored('\nCreating ENA Metadata Table:', 'section_header', plain_text_bool))
        print(rich_text_colored('\nCombining results:', 'general_text', plain_text_bool))
        print('Left Join by Run Accessions using ENA Technical Metadata Table as Reference')
        
        #Previous steps
        ##Get outfile name
        out_name_ena = study_accession + '_ENA_metadata.tsv'
        ##Treat output_directory parameter / Get full output file path
        outputfile_ena = treat_output_directory_parameter_outfiles(out_name_ena, outputdir_path)
        
        #Join both tables by run accession (left join)
        df_merge_ENA = pd.merge(tech_table, mgtoolkit_table, left_on = 'run_accession', right_on = 'Run', how = 'left')
        
        #Show saved file message
        print(rich_text_colored('\nSaving results in file:', 'general_text', plain_text_bool))
        print(outputfile_ena)
        #Save df_merge_ENA
        df_merge_ENA.to_csv(outputfile_ena, header = True, index = False, sep = '\t')
        
        
    except Exception as ex:
        print(rich_text_colored('\nThe system returned the following exception:\n', 'exception', plain_text_bool), ex)
    
    finally:
        #Print empty line for aesthetic purposes
        print('')
        
if __name__ == '__main__':
    main()
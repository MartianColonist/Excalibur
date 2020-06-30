#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 14:08:06 2020

@author: arnav
"""

from bs4 import BeautifulSoup
import requests
import pandas
import os
from hapi import moleculeName, isotopologueName

def HITEMP_table():
    url = 'https://hitran.org/hitemp/'

    web_content = requests.get(url).text
    
    soup = BeautifulSoup(web_content, "lxml")
    table = soup.find('table')

    n_rows = 0
    n_columns = 0
    column_names = []

    for row in table.find_all('tr'):
        td_tags = row.find_all('td')
        if len(td_tags) > 0:
            
            n_rows += 1
            
        if n_columns == 0: # Use the number of td tags in the first row to set the first column
            n_columns = len(td_tags)
            
            # Handle column names
            th_tags = row.find_all('th') 
            
        if len(th_tags) > 0 and len(column_names) == 0:
            for th in th_tags:
                column_names.append(th.get_text())
            
         
    hitemp = pandas.DataFrame(columns = column_names, index= range(0,n_rows))

    row_marker = 0
    for row in table.find_all('tr'):
        column_marker = 0
        columns = row.find_all('td')
        for column in columns:
            hitemp.iat[row_marker,column_marker] = column.get_text()
            column_marker += 1
        
        if len(columns) > 0:
            row_marker += 1
        

    hitemp = hitemp[:-1]
    hitemp.rename(columns = {'Iso counta':'Iso Count'}, inplace = True)
    hitemp.loc[len(hitemp)] = ['4', 'N2O', 'Nitrous Oxide', '5', '3626425', '0', '12899', '2019', '']
    hitemp.loc[:, 'ID'] = pandas.to_numeric(hitemp['ID'])
    hitemp.loc[:, 'Iso Count'] = pandas.to_numeric(hitemp['Iso Count'])


    hitemp.sort_values(by = 'ID', inplace = True)
    hitemp.reset_index(drop = True, inplace = True)


    counter = 0
    for tag in table.find_all('a'):
        hitemp.loc[counter, 'Download'] = 'https://hitran.org' + tag.get('href')
        counter += 1
        
    return hitemp


def create_directories(mol_ID, iso_ID):
    """
    Create new folders to store the relevant data

    Parameters
    ----------
    mol_ID : TYPE
        DESCRIPTION.
    iso_ID : TYPE
        DESCRIPTION.

    Returns
    -------
    line_list_folder : TYPE
        DESCRIPTION.

    """
    
    input_folder = '../input'
    molecule_folder = input_folder + '/' + moleculeName(mol_ID) + '-' + isotopologueName(mol_ID, iso_ID)
    line_list_folder = molecule_folder + '/HITEMP'
    
    if os.path.exists(input_folder) == False:
        os.mkdir(input_folder)
    
    if os.path.exists(molecule_folder) == False:
        os.mkdir(molecule_folder)

    if os.path.exists(line_list_folder) == False:
        os.mkdir(line_list_folder)
        
    return line_list_folder


def download_line_list(mol_ID, iso_ID):
    table = HITEMP_table()
    row = table.loc[table['ID'] == mol_ID]
    
    download_link = row.loc[row.index.values[0], 'Download']
    print(download_link)


def summon_HITEMP(molecule, isotopologue):
    output_folder = create_directories(molecule, isotopologue)
    

table = HITEMP_table()
row = table[table['ID'] == 4]
print(row.index.values[0])
#print(row.loc[:, 'Iso Count'])
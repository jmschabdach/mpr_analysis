"""
This script collates freesurfer results in to a single SQLite database.
There are tables for each file defined
"""

import sqlite3
from glob import glob
import pandas as pd
import argparse
import os

##
# Open a connection to a sqlite3 database file
# @param filename A str specifying the path to the sqlite3 .db file
# @returns CUR A cursor object for the sqlite3 database connection
def connectToDatabase(filename):
    # start database
    CONNECTION = sqlite3.connect(filename)
    CUR = CONNECTION.cursor()

    return CUR, CONNECTION


##
# Close the connection to a sqlite3 database
# @param CUR A cursor object from a sqlite3 database connection
def closeDatabase(CONNECTION):
    CONNECTION.close()


##
# Get a list of image ids (single image selected per subject)
# @param CUR A cursor object for the sqlite3 database connection
# @returns 
def getSingleImageIdPerSubject(CUR):
    images = []
   
    # For each unique subject in the sessionData table
    for subject in [subj[0] for subj in CUR.execute("SELECT DISTINCT subj FROM sessionData").fetchall()]:
       print(subject)

       # For each scan session belonging to the subject
       for session in [sess[0] for sess in CUR.execute("SELECT DISTINCT session FROM sessionData WHERE subj = \""+str(subject)+"\"").fetchall()]:
           print(session)

           # Get the directories containing the FreeSurfer output for each image ("run"s)
           runs = sorted([run[0] for run in CUR.execute("SELECT run FROM sessionData WHERE subj =\""+str(subject)+"\" AND session = \""+str(session)+"\"")])
           print("Number of runs:", len(runs))

           # Select the first run in the sorted list as the image for that subject
           if len(runs) >= 1:
               images.append(runs[0])

    # Get the id of the session (the ID common across all tables in the database) using the list of images
    ids = [CUR.execute("SELECT session_id FROM sessionData WHERE run = \""+str(img)+"\"").fetchall()[0] for img in images]

    # Convert the list of ids into a string that can be used to query other tables
    idsString = "("+",".join("%s" % tup[0] for tup in ids)+")"

    return idsString


def getMeasureStatsAsDf(CUR, idsString):
    # Get the measure labels, values, and image ids from the Measures table if the image ids are in the idsString sql list
    measures = CUR.execute("SELECT Measure, value, session_id FROM Measures WHERE session_id IN "+idsString).fetchall()

    # Convert the list of tuples into a dataframe with the structure (attribute, value, entity) 
    df = pd.DataFrame(measures, columns=['measure', 'value', 'image_id'])
    # Reshape the datafram from an EAV structure to a table where each column is a measure, each row is the id of an image, and the cells are the values
    df = df.pivot(index='image_id', columns='measure', values='value')

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--dbfile', help='Path to database file')

    args = parser.parse_args()
    dbfn = args.dbfile

    # Open the database
    cursor, connection = connectToDatabase(dbfn)

    # Get the image ids
    idsStringForQuery = getSingleImageIdPerSubject(cursor)

    # Get the measures for the image ids in a dataframe
    measuresDf = getMeasureStatsAsDf(cursor, idsStringForQuery)

    # Close the database (important! do before ending the file and after all "execute" commands)
    closeDatabase(connection)

    # Print the measures dataframe
    print(measuresDf)


if __name__ == "__main__":
    main()
    

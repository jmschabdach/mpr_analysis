"""
This script collates freesurfer results in to a single SQLite database.
There are tables for each file defined
"""

import sqlite3
from glob import glob
import pandas as pd
import argparse
import os

class Subj:
    """
    Class function for handling the freesurfer data for an individual
    """
    def __init__(self, subjpath):
        # define subject path
        self.subjpath = subjpath

        # define subject, session and run
        self.subj = subjpath.split("/")[-3]
        self.session = subjpath.split("/")[-2]
        self.run = subjpath.split("/")[-1]

        # blank session ID
        self.sess_id = None

    def write_session_data(self, sess_id, cursor, conn):
        """ Writes data for each session """
        sess_dict = {"subj": self.subj,
                     "session": self.session,
                     "run": self.run,
                     "session_id": str(sess_id)}

        self.sess_id = sess_id

        # write session information
        cursor.execute("INSERT INTO sessionData VALUES('" +\
                       "', '".join([str(sess_dict[header]) for header in\
                                   ["subj", "session", "run", "session_id"]])\
                       + "')")
        conn.commit()

    def write_to_results(self, file_name, side_lr, column_names, conn):
        """ Write data to results file """
        if side_lr:
            file_name_path = os.path.join(self.subjpath, "stats",
                                       '.'.join([side_lr, file_name, "stats"]))
        else:
            file_name_path = os.path.join(self.subjpath, "stats",\
                                   '.'.join([file_name, "stats"]))

        # import table using pandas
        data_frame = pd.read_csv(file_name_path, sep='\s+',
                                 comment="#",
                                 header=None)
        data_frame.columns = column_names

        # add some columns
        data_frame['side'] = side_lr
        data_frame['session_id'] = self.sess_id

        # write table to sql database
        data_frame.to_sql(file_name, conn, if_exists="append")
        conn.commit()

    def write_measures(self, conn):
        '''
        Writes specific meaures from aseg.stats, including estimated total intracranial volume
        '''
        aseg_path = os.path.join(self.subjpath, "stats/aseg.stats")

        aseg_open = open(aseg_path, "r")
        measures = {line.split(",")[0].split()[2]: [line.split(",")[0].split()[2],
                                                    line.split(",")[3].strip(' ')]
                    for line in aseg_open.readlines() if 'Measure' in line}
        aseg_open.close()

        data_frame = pd.DataFrame.from_dict(measures,
                                            orient="index",
                                            columns=["Measure", "Value"])

        # add session id
        data_frame['SESSION_ID'] = self.sess_id

        # write table to sql database
        data_frame.to_sql("Measures", conn, if_exists="append")
        conn.commit()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', help='Path to FreeSurfer output directory')

    args = parser.parse_args()

    bidspath = args.path
    print(bidspath)

    # set headings of freesurfer tables
    COL_NAMES = ["StructName",
                 "NumVert",
                 "SurfArea",
                 "GrayVol",
                 "ThickAvg",
                 "ThickStd",
                 "MeanCurv",
                 "GausCurv",
                 "FoldInd",
                 "CurvInd"]
    
    # freesurfer filenames
    FNAMES = ["aparc.a2009s",
              "aparc.DKTatlas",
              "aparc"]
    
    # set up brainstem values
    COL_NAMES_BS = ["index_fs", "number", "blank", "volume", "StrucName"]
    
    # set up brainnetome brainstem column names
    COL_NAMES_ASEG = ["Index_fs",
                      "SegId",
                      "NVoxels",
                      "Volume_mm3",
                      "StructName",
                      "normMean",
                      "normStdDev",
                      "normMin",
                      "normMax",
                      "normRange"]
    
    # set up SQL database
    # start database
    CONNECTION = sqlite3.connect("fsResults.db")
    CUR = CONNECTION.cursor()
    
    CUR.execute("CREATE TABLE IF NOT EXISTS sessionData (subj, session, run, session_id)")
    CONNECTION.commit()
    
    # get a list of existing subjects/sessions
    CUR.execute("SELECT subj, session FROM sessionData")
    EXISTING = [subsess[0] + subsess[1] for subsess in CUR.fetchall()]
    
    SESSION_ID = len(EXISTING) + 1
    # iterate through all subjects
    fsoutdirs = os.path.join(bidspath, "sub-*/ses-*/sub-*_ses-*")
    for subj_dir in glob(fsoutdirs, recursive=True):
        print(subj_dir)
        fs_session = Subj(subj_dir)
        if not fs_session.subj + fs_session.session in EXISTING:
            print(subj_dir)
            fs_session.write_session_data(SESSION_ID, CUR, CONNECTION)
    
            # iterate through filenames
            for fname in FNAMES:
                for side in ["lh", "rh"]:
                    try:
                        fs_session.write_to_results(fname, side, COL_NAMES, CONNECTION)
                    except OSError as err_os:
                        print(err_os)
    
            try:
                fs_session.write_to_results("aseg", None, COL_NAMES_ASEG, CONNECTION)
                fs_session.write_measures(CONNECTION)
            except OSError as err_os:
                print(err_os)
    
            SESSION_ID += 1
    
    CONNECTION.close()


if __name__ == "__main__":
    main()
    

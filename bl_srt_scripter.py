#!/usr/local/anaconda3/bin/python
'''

Creat a quick script that chooses stars based on LST
for the Sardinia Radio Telescope

Author: Howard 
Edits: Vishal
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import re
from argparse import ArgumentParser
import optparse


if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.add_option("--time", action='store', dest='utime', type=str,help="Input UTC time in format inside quote \"YYYY-MM-DD HH:MM:SS\"")

    options,args = parser.parse_args()
    utime =  options.utime

    print(utime)
    print(type(utime))

    if utime:
        time = Time(utime)
    else:
        time = Time.now()

    SRT_obs = EarthLocation(lat=39.4930*u.deg, lon=-9.2451*u.deg, height=600*u.m)
    #time = Time('2021-05-01 08:00:00') #- utcoffset # Time of observation
    #time = Time('2024-03-16 15:30:00') # Time in UTC 
    timestart = time

    buffer = 5 # in minutes between observations
    session_length = 9.*u.hour # 6 hours in minutes
    #times = time * session_length

    # Read in the relevant starlists
    #FILE_IN1 = "SRT_GBT_list.csv"
    #FILE_IN1 = "SRT_south_list_edited.csv"
    #FILE_IN1 = "Simultaneous_SETI_TESS_targets.csv"
    FILE_IN1 = "All_TESS_north.csv"
    #OUTFILE = 'outfile.csv'

    df1 = pd.read_csv(FILE_IN1,index_col=False)
    df1.drop_duplicates('TIC ID',inplace=True)

    # sort the stars in each by nearest to LST.
    #test = df1['RA (deg)']* u.deg
    #df1['RA_hr'] = df1['RA (deg)'] / 15.*u.hour
    #df1['Dec (deg)'] * u.deg
    #targets = SkyCoord(df1['RA (deg)'].iloc[:]*u.deg,df1['Dec (deg)'].iloc[:]*u.deg,frame='icrs')

    ra_quantities = [ra * u.hourangle for ra in df1['RA (hr)'].tolist()]
    dec_quantities = [dec * u.deg for dec in df1['Dec (deg)'].tolist()]
    targets = SkyCoord(ra_quantities, dec_quantities, frame='icrs')
    
    t_altaz = targets.transform_to(AltAz(obstime=time,location=SRT_obs))
    print(t_altaz)

    # At beginning of session, identify hightest source:

    df_observe = pd.DataFrame(columns=df1.columns)

    n_targets = np.arange(session_length.value *2)

    for item in n_targets:
        #targets = SkyCoord(df1['RA (deg)'].iloc[:]*u.deg,df1['Dec (deg)'].iloc[:]*u.deg,frame='icrs')
        #targets = SkyCoord(df1['RA (deg)'].iloc[:]*u.deg,df1['Dec (deg)'].iloc[:]*u.deg,frame='icrs')
        t_altaz = targets.transform_to(AltAz(obstime=time,location=SRT_obs))
        highest = np.argmax(t_altaz.alt)
        print(highest)
        df_observe = df_observe._append(df1.iloc[highest])
        #print(time.value,t_altaz[highest].alt.value,df1['RA_hr'].iloc[highest],df1['Dec (deg)'].iloc[highest],df1['TIC ID'].iloc[highest])
        df1.drop(df1.index[highest],inplace=True)
        ra_quantities = [ra * u.hour for ra in df1['RA (hr)'].tolist()]
        dec_quantities = [dec * u.deg for dec in df1['Dec (deg)'].tolist()]
        targets = SkyCoord(ra_quantities, dec_quantities, frame='icrs')
        time = time+0.5*u.hour


    outtxtfile = []
    outlistfile = []
    
    new_ra_quantities = [ra * u.hourangle for ra in df_observe['RA (hr)'].tolist()]
    new_dec_quantities = [dec * u.deg for dec in df_observe['Dec (deg)'].tolist()]
    newtargets = SkyCoord(new_ra_quantities,new_dec_quantities,frame='icrs')
    df_observe['RA (deg)'] = newtargets.ra.degree

    for ind,d in df_observe.iterrows():   
        for i in range(3):
            #print "TIC" + str(int(d['TIC ID'])) + "_ON  AlwaysOn TOTALPOWER EQ " + str(d['RA (deg)']) + "d " + str(d['Dec (deg)']) + "d"
            outtxtfile.append(str(d['TIC ID']) + "_ON  AlwaysOn TOTALPOWER EQ " + str(d['RA (deg)']) + "d " + str(d['Dec (deg)']) + "d")

            #print "TIC" + str(int(d['TIC ID'])) + "_OFF  AlwaysOff TOTALPOWER EQ " + str(d['RA (deg)']) + "d " + str(d['Dec (deg)']) + "d"
            if(d['Dec (deg)'] + 5.0 < 90.0):
                outtxtfile.append(str(d['TIC ID']) + "_OFF  AlwaysOn TOTALPOWER EQ " + str(d['RA (deg)']) + "d " + str(d['Dec (deg)']+5.0) + "d")
            else:    
                outtxtfile.append(str(d['TIC ID']) + "_OFF  AlwaysOn TOTALPOWER EQ " + str(d['RA (deg)']) + "d " + str(d['Dec (deg)']-5.0) + "d")


    outfilename="targets_"+"_".join(re.split('-|:| ',timestart.value))+".txt"

    with open(outfilename, 'w') as f:
        for item in outtxtfile:
            print(item, file=f)

    #with open(outfilename, 'w') as f:
    #    for item in outfile:
    #        print >> f, item
    #df_observe.to_csv('outfile.csv')


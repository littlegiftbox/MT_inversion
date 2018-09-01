#!/usr/bin/env python3

import sys
import csv
import os
import subprocess
from obspy.clients.iris import Client
from obspy import UTCDateTime

###############################################
# Step 0. Configuration
net_name='HV'  #network name
sta_list=['HSSD','ERZ1','HOVE','HUAD'] #station name
chn_list=['HHE','HHN','HHZ']


if len(sys.argv)!=2:
    print ('Usage: python mtinversion.py process_stage')
    sys.exit()

process_stage=int(sys.argv[1])

# 1: make Green's function
# 2: read cataloig and download data
# 3: add location + change oz files
# 4: preprocess sac file and make .dat file for inversion
# 5: do moment tensor inversion & plot

# b2s.par, model_hvo.d, dfile, mt_inv.in, run_cps, 
# run_cps_2_helm, run_filtsyniso and catalog should be in the folder


# AFTER ALL THE PROCEDURE
# NEED TO ADJUST ZCOR MANUALLY

# 6: based on the zcor/quality on each station, grid search for quality = 4
#    make sac file, do cross-correlation to find the best zcor for quality < 4
#    especially zcor < 0

################################################

#Step 1. Make Green's Function 

os.system('mkdir GREEN_FUNCTIONS')
os.chdir('GREEN_FUNCTIONS')

if process_stage <= 1:
    os.system('cp ../b2s.par .')
    os.system('../run_cps')
    os.system('../run_cps_2_helm')
    disp_list=os.popen('ls *.disp').readlines()
    print(disp_list)
    for disp in disp_list:
        disp=disp.strip()
        new_disp=disp[:-5]
        print(new_disp)
        os.system('../run_filtsyniso {} {}'.format(disp,new_disp))

os.chdir('..')

# Step 2. Read csv catalog downloaded from USGS
csvFile = open ("Hawaii4.5_June.csv","r")
dict_reader = csv.DictReader(csvFile)
catalog=dict_reader
os.putenv("SAC_DISPLAY_COPYRIGHT", "0")



for item in catalog:
    #print(item)
    if item["type"]=="earthquake":
        continue
    if float(item["mag"]) < 5.0:
        continue

    date_str=(item["time"]).split('T')
    ymd_list=(date_str[0]).split('-')
    # make folder for each day 
    # assume there is only 1 Mw 5.0 Earthquake per day
    os.system('mkdir '+date_str[0]) 
    os.chdir(date_str[0])

    print('Now doing '+date_str[0]+', magnitude '+item["mag"])
    
    os.system('mkdir MT_INVERSION')
    os.system('cp ../GREEN_FUNCTIONS/{}* MT_INVERSION'.format(net_name.lower()))
    os.system('cp ../mt_inv.in MT_INVERSION')

    os.system('mkdir RAW')
    os.chdir('RAW')

    #download sac file & instrument responce
    if process_stage <= 2: 
        for sta_name in sta_list:
            for chn_name in chn_list:
                #loop over sta/channels
                sac_name='SAC.'+ymd_list[0]+ymd_list[1]+ymd_list[2]+'.'+net_name+'.'+sta_name+'.'+chn_name;
                pz_name=net_name+'.'+sta_name+'.pz'
                t0=UTCDateTime(item["time"])
                client = Client()
                client.timeseries(net_name, sta_name, "", chn_name, t0-60, t0+5*60, filename=sac_name,output='sacbb')
                client.sacpz(net_name, sta_name, "", chn_name, t0, filename=pz_name)

    if process_stage <= 3:

        #get station location in pz file 
        #meanwhile delete 1 zeros, overwrite the file
        for sta_name in sta_list:

            pz_name=net_name+'.'+sta_name+'.pz'
            pz_file=open(pz_name,'r')
            pz_lines=pz_file.readlines()
            pz_file.close()
            pz_file=open(pz_name,'w')
            flag=0
            for line in pz_lines:

                #skip empty lines
                if len(line) == 1:
                    continue

                new_line=line.split()
                if line.startswith('*'):
                    if new_line[1]=='LATITUDE':
                        sta_lat=float(new_line[3])
                    if new_line[1]=='LONGITUDE':
                        sta_lon=float(new_line[3])
                    if new_line[1]=='ELEVATION':
                        sta_ele=float(new_line[3])
                    continue
                

                if new_line[0] == 'ZEROS':
                    pz_file.write("ZEROS\t"+str(int(new_line[1])-1)+"\n")
                    flag=1
                    continue

                if flag==1:
                    flag=0
                    continue

                pz_file.write(line)

            pz_file.close()

            #add event location and station location
            s = "rh SAC*{}* \n".format(sta_name)
            s += "ch evlo {} evla {} evdp {} mag {} \n".format(float(item["longitude"]), 
                                                           float(item["latitude"]), 
                                                           float(item["depth"]), 
                                                           float(item["mag"]))
            s += "ch stlo {} stla {} stel {} \n".format(sta_lon,sta_lat,sta_ele) 
            s += "wh \n"
            s += "q \n"
            subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s.encode())


    #preprocess the data & make dat files
    if process_stage <= 4:

        for sta_name in sta_list:

            pz_name=net_name+'.'+sta_name+'.pz'
            
            #change header
            s = "rh SAC*{}*{} \n".format(sta_name,chn_list[0])
            s += "ch cmpaz 90 cmpinc 90 \n"
            s += "wh \n"
            s += "rh SAC*{}*{} \n".format(sta_name,chn_list[1])
            s += "ch cmpaz 0 cmpinc 90 \n"
            s += "wh \n"
            #rotate
            s += "cut 0 300 \n"
            s += "r SAC*{}*{} SAC*{}*{} \n".format(sta_name,chn_list[0],sta_name,chn_list[1]) #read E,N component
            s += "rot to gcp \n"
            s += "r more SAC*{}*{} \n".format(sta_name,chn_list[2]) #z component
            s += "rmean \n"
            s += "transfer from polezero s {} \n".format(pz_name)
            #int, mul, filter and etc
            s += "int \n"
            s += "mul 100 \n"
            s += "bp co 0.02 0.05 n 6 p 1 \n"
            s += "interpolate delta 0.1 \n"
            s += "w tmp2 tmp1 tmp3 \n"
            s += "q \n"
            subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s.encode())
            os.system('sac2helm out={}'.format(sta_name+'.dat'))
            os.system('cp *.dat ../MT_INVERSION')
            #rename tmp2, tmp1, tmp3
            os.system('mv tmp2 {}.{}.HHR.SAC'.format(net_name,sta_name))
            os.system('mv tmp1 {}.{}.HHT.SAC'.format(net_name,sta_name))
            os.system('mv tmp3 {}.{}.HHZ.SAC'.format(net_name,sta_name))

    #do moment tensor inversion and 
    if process_stage <= 5:
        os.chdir('../MT_INVERSION')
        os.system('tdmt_invc_iso mt_inv.in 2> mt_inv.pyout')
        os.system('tdmt_plot_gmt5.perl mtinvout')


    #go back to parent directory
    os.chdir('../..')

csvFile.close()




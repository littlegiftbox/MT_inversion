import os
import obspy.signal
from obspy.core import read
from obspy.signal.cross_correlation import correlate, xcorr_max
import matplotlib.pyplot as plt
import numpy as np

###############################################
# Step 0. Configuration
net_name='HV'  #network name
sta_list=['HSSD','ERZ1','HOVE','HUAD'] #station name 
chn_list=['HHE','HHN','HHZ']

# IMPORTANT: the sequence of sta_list should be the same as mt_inv.in 
# For really bad results, sometime we still cannot find a reasonable answer
# to reset autozcor, run "python mtinversion.py 5" in the terminal 

################################################
# Functions

def get_tdmt_info(i):
    auto_zcor=(os.popen('grep {} mt_inv.out'.format(obsrv_name)).readlines())[-1].split('=')[-1].strip()
    theo_zcor=(float(dist)/8+60)*10; #delta=0.1, cut 60s before original time
    VR=(os.popen('grep Station mt_inv.pyout').readlines())[-5+i].split('=')[1].split()[0]
    VR_all=(os.popen('grep VarRed mt_inv.out').readlines())[-1].split('=')[-1].replace('+','').strip()
    return auto_zcor,theo_zcor, VR, VR_all


def change_zcor(obsrv_name,dist,azimuth,new_zcor,pts,i):
    #Change zcor and overwrite mt_inv.in
    f=open('mt_inv.in','r')
    mt_in=f.readlines()
    f.close()

    f=open('mt_inv.in','w')
    for j in range(0,len(mt_in)):
        if j!=i:
            f.write(mt_in[j])
        else:
            f.write('{} {} {} {} {}\n'.format(obsrv_name,dist,azimuth,int(new_zcor),pts))
    f.close()    

def find_best_zcor(old_zcor,step,step_len,VR,obsrv_name,dist,azimuth,pts,i):
    #search range old_zcor+-step*step_len
    best_vr=VR
    best_zcor=old_zcor
    #np.linspace: if not specify, default step=50
    new_zcors=np.linspace(float(old_zcor)-step*step_len,float(old_zcor)+step*step_len,step_len*2+1)
    for new_zcor in new_zcors:
        change_zcor(obsrv_name,dist,azimuth,new_zcor,pts,i)
        os.system('tdmt_invc_iso mt_inv.in 2> mt_inv.pyout')
        tmp,tmp,new_vr,tmp=get_tdmt_info(i)
        if float(new_vr) > float(best_vr):
            best_vr=new_vr;
            best_zcor=new_zcor;
    print('Step=',step,'new_zcor=',best_zcor,'new_vr=',best_vr) 
    return best_zcor,best_vr   


##################################################

#Loop over each event folder
events=os.listdir(os.getcwd())
for event in events:
    if event[4]!='-' or event[7]!='-':
        continue
    #go into the event folder
    os.chdir(event)
    os.chdir('MT_INVERSION')
    os.system('cp ../../b2s.par .')

    #in the event folder, read mt_inv.in 
    mt_file=open('mt_inv.in','r')
    mt_in=mt_file.readlines()
    mt_file.close()

    #grep info from mt_inv.out
    strike=(os.popen('grep Strike mt_inv.out').readlines())[-1].split('=')[1].split()[0]
    rake=(os.popen('grep Rake mt_inv.out').readlines())[-1].split('=')[1].split()[0]
    dip=(os.popen('grep Dip mt_inv.out').readlines())[-1].split('=')[1].split()[0]
    moment=(os.popen('grep Mo mt_inv.out').readlines())[-1].split('=')[1].split()[0].replace('+','')
    moment=moment[:4]+moment[-3:]
    print('\n')
    print('Working on event ',event)
    print('Strike, rake, dip, moment: ',strike,rake,dip,moment)


    sta_num=len(sta_list)


    #read info and fix auto-generated zcor
    for i in range(1,sta_num+1):

        obsrv_name,dist,azimuth,zcor,pts=(mt_in[i]).split()
        model_name=((mt_in[sta_num+i]).split())[0]
        auto_zcor,theo_zcor,VR,VR_all=get_tdmt_info(i)
        change_zcor(obsrv_name,dist,azimuth,auto_zcor,pts,i)



    for i in range(1,sta_num+1):

        obsrv_name,dist,azimuth,zcor,pts=(mt_in[i]).split()
        model_name=((mt_in[sta_num+i]).split())[0]
        print('Working on station ',model_name)

        auto_zcor,theo_zcor,VR,VR_all=get_tdmt_info(i)
        print('auto_zcor=', auto_zcor,' theo_zcor=',theo_zcor,' VR=',VR,' Average VR=',VR_all)

        #search for half-wavelength, bp filter 0.02~0.05 Hz, 20~50s period
        #now I only test for [-50,50] zcor deviated from theotical zcor

        #define searching step-width based on VR quality
        #another thing that is worth exploring: 
        #the search bound corresponding to each zcor
        if float(VR)>90:
            #good quality, minor shift
            #step=1 step_len=5
            best_zcor,best_vr=find_best_zcor(auto_zcor,1,5,VR,obsrv_name,dist,azimuth,pts,i)

        if float(VR)>70: 
            best_zcor,best_vr=find_best_zcor(auto_zcor,5,10,VR,obsrv_name,dist,azimuth,pts,i)
            best_zcor,best_vr=find_best_zcor(best_zcor,1,5,best_vr,obsrv_name,dist,azimuth,pts,i)

        elif float(VR)>50:
            best_zcor,best_vr=find_best_zcor(auto_zcor,5,20,VR,obsrv_name,dist,azimuth,pts,i)
            best_zcor,best_vr=find_best_zcor(best_zcor,1,5,best_vr,obsrv_name,dist,azimuth,pts,i)

        else:
            #bad quality, use theoretical zcor instead
            best_zcor,best_vr=find_best_zcor(theo_zcor,5,50,VR,obsrv_name,dist,azimuth,pts,i)
            best_zcor,best_vr=find_best_zcor(best_zcor,1,5,best_vr,obsrv_name,dist,azimuth,pts,i)




    #run and plot new moment tensor results
    os.system('tdmt_invc_iso mt_inv.in 2> mt_inv.pyout')
    os.system('tdmt_plot_gmt5.perl mtinvout')
    os.chdir('../..')





"""
#obsolete cross-correlation code

os.chdir('2018-05-17/MT_INVERSION/')

st_model = read('hvo_30.0d0.6.ver')
st_obsrv = read('../RAW/HV.HSSD.HHZ.SAC')
# type(st) = obspy.stream.Stream
data_model = st_model[0].data
data_obsrv = st_obsrv[0].data
max_shift=1000 #pts
cc = correlate(data_model,data_obsrv,max_shift) 
cc_time = np.linspace(-max_shift,max_shift,2*max_shift+1)
shift, value = xcorr_max(cc)
print('cc: ',cc)
print(shift,value)

plt.figure(1)
plt.plot(data_model/max(data_model))
plt.plot(data_obsrv/max(data_obsrv))
plt.xlabel('data points')
plt.ylabel('disp?') 
plt.savefig('waveforms.png')


plt.figure(2)
plt.plot(cc_time,cc)
plt.xlabel('data points')
plt.ylabel('xcorr coefficient')
plt.savefig('xcorr.png')

plt.figure(3)
plt.plot(np.roll(data_model,-shift*2)/max(data_model))
plt.plot(data_obsrv/max(data_obsrv))
plt.xlabel('data points')
plt.ylabel('disp?') 
plt.savefig('waveforms_shift.png')


"""


# Telecommunication-Radar-system
 Radar is a detection system that uses radio waves to determine the range,  angle, or velocity of objects. It can be used to decry ships, missiles, motor vehicles, weather formation, spacecraft, aircraft and so on.
import scipy.io as spio
import os
import matplotlib.pyplot as plt
import numpy as np
# import datetime
# DataPath='D:/Python/Uni_HWI/scripts20/RAW/'
DataPath='D:/Hochschule Wismar/stud.ip/Sem 
2/Lecture_Nachrichtentechnisches_Projekt_-_2020 (3)/RAW data/'
Files=os.listdir(DataPath)
currentfile=str(DataPath)+str(Files[12])
# importing MATLAB mat file (containing radar raw data)
mat = spio.loadmat(currentfile, squeeze_me=True)
datenums=mat['datenums']
ranges=mat['ranges']
data=mat['data']
# datenums ~ days since year 0
# here only the time is important for us -> hours, minutes, seconds
# => fraction / remainder of the integer
t=(datenums-np.floor(np.min(datenums)))*24
# number of range gates, data points, receivers
noRG=np.size(data,0)
noDP=np.size(data,1)
noRx=np.size(data,2)
RXsel=0;
# time series 
y=data [1,:,RXsel]
# I/Q diagram
plt.figure()
plt.plot(np.real(y),np.imag(y),'.')
plt.xlabel('real')
plt.ylabel('imag')
plt.title('I/Q diagram for one receiver and one range with signal')
# power of the complex valued voltages (measured by the receiver)
PWR=20*np.log10(np.abs(data[:,:,RXsel]))
#PWR=20*np.log10(np.abs(y))
PWR.shape
Type (PWR)
plt.figure()
plt.plot(PWR[:,0],ranges)
plt.ylabel('range /km')
plt.xlabel('power /dB')
plt.title('Power plot-Line profile for 1st receiver and range')
plt.figure()
plt.pcolor(t,ranges,PWR)
plt.xlabel('time / HH')
plt.ylabel('range /km')
plt.title('power /dB')
plt.clim(20,70)
plt.colorbar()
# combine the data of all three receivers - complex sum!
PWR=np.zeros([noRG,noRx])
for rx in range(noRx):
 PWR[:,rx]=20*np.log10(np.abs(np.mean(data[:,:,rx],1)))
plt.figure()
plt.plot(PWR,ranges)
plt.legend(['RX1','RX2','RX3'])
plt.xlabel('power /dB')
plt.ylabel('ranges /km')
datacomb=np.sum(data,2)
# power for the combined receivers
PWRcomb=20*np.log10(np.abs(datacomb))
plt.figure()
plt.pcolor(t,ranges,PWRcomb)
plt.xlabel('time / HH')
plt.ylabel('range /km')
plt.title('Power for combined receivers for all ranges /dB')
plt.clim(20,70)
plt.colorbar()
def make_fft(t,y):
 dt = t[1]-t[0] # dt -> temporal resolution ~ sample rate
 f = np.fft.fftfreq(t.size, dt) # frequency axis
 Y = np.fft.fft(y) # FFT
 f=np.fft.fftshift(f)
 Y= np.fft.fftshift(Y)/(len(y))
 return f,Y
tsec=t*60*60
f,spec=make_fft(tsec,data[-25,:,RXsel])
plt.figure()
plt.plot(f,10*np.log10(abs(spec)))
plt.grid ('on')
plt.xlabel('f / Hz')
plt.ylabel('amplitude /dB')
# Spectra for all ranges and all receivers
Spectr=np.zeros([noRG,noDP,noRx])+1j*np.zeros([noRG,noDP,noRx])
for rx in range(noRx):
 for rg in range(noRG):
 f,Spectr[rg,:,rx]=make_fft(tsec,data[rg,:,rx])
plt.figure()
for rx in range(noRx):
 plt.subplot(1,3,rx+1)
 plt.pcolor(f,ranges,10*np.log10(abs(Spectr[:,:,rx])))
 plt.clim([-15, 15])
 plt.xlabel('f /Hz')
 plt.ylabel('range /km')
 plt.xlim([-1,1])
 plt.title('Spectra Amplitude')
 # plt.colorbar()
phases=[]
plt.figure() 
for rx in range(noRx): # rx=0
 plt.subplot(1,3,rx+1)
 phases=np.angle(Spectr[:,:,rx])/np.pi*180
 plt.pcolor(f,ranges,phases,cmap='jet')
 plt.clim([-180, 180])
 plt.colorbar()
 plt.title('Spectra Phase') 
 
# perform coherent integrations
def make_ci(t,y, ci):
 nptsn=int(np.floor(len(y)/ci))
 yn=np.empty(nptsn)+1j*np.empty(nptsn)
 tn=np.empty(nptsn)
 for i in range(0,nptsn):
 yn[i]=np.mean(y[i*ci:i*ci+ci-1])
 tn[i]=np.mean(t[i*ci:(i+1)*ci])
 return tn,yn
# make FFT spectrum, frequency axis
def make_fft(t,y):
 dt = t[1]-t[0] # dt -> temporal resolution ~ sample rate
 f = np.fft.fftfreq(t.size, dt) # frequency axis
 Y = np.fft.fft(y) # FFT
 f=np.fft.fftshift(f)
 Y= np.fft.fftshift(Y)/(len(y))
 return f,Y
# number of coherent integrations of I/Q raw data (time series)
ci=2 #low noise as comapred to below Ci values
#ci=12 Moderate Noise
#ci=30 After performing CI, we found that Noise is very large
y=data[-25,:,0]
tn,yn=make_ci(t,y,ci)
plt.figure()
plt.subplot(1,2,1)
plt.plot(np.real(y),np.imag(y),'-*')
plt.xlim([-100,100])
plt.ylim([-100,100])
plt.subplot(1,2,2)
plt.plot(np.real(yn),np.imag(yn),'-*')
plt.xlim([-100,100])
plt.ylim([-100,100])
# lentgh of the "new" integrated time series
noDPn=int(np.floor(noDP/ci))
# predefine matrix for integrated raw data
datan=np.zeros([noRG,noDPn,noRx])+1j*np.zeros([noRG,noDPn,noRx])
for rx in range(noRx):
 for rg in range(noRG):
 tn,datan[rg,:,rx]=make_ci(t,data[rg,:,rx],ci)
#Powwer-combined recievers , time series and all ranges
datacomb=np.sum(datan,2)
PWRcomb=20*np.log10(np.abs(datacomb))
plt.figure()
plt.pcolor(tn,ranges,PWRcomb)
plt.xlabel('time / HH')
plt.ylabel('range /km')
plt.title('combined power /dB')
plt.clim(20,70)
plt.colorbar()
# Cross-Spectra for all ranges and all receivers
XSpectr=np.zeros([noRG,noDP,noRx])+1j*np.zeros([noRG,noDP,noRx])
XSpectr[:,:,0]=Spectr[:,:,0]*np.conj(Spectr[:,:,1])
XSpectr[:,:,1]=Spectr[:,:,0]*np.conj(Spectr[:,:,2])
XSpectr[:,:,2]=Spectr[:,:,1]*np.conj(Spectr[:,:,2])
plt.figure()
plt.pcolor(f,ranges,10*np.log10(abs(XSpectr[:,:,1]))/2)
plt.pcolor(f,ranges,np.angle(XSpectr[:,:,1]))
plt.colorbar()
plt.clim([-15, 15])
plt.figure()
for rx in range(noRx):
 plt.subplot(1,3,rx+1)
 plt.pcolor(f,ranges,10*np.log10(abs(XSpectr[:,:,rx])/2),cmap='jet')
 plt.clim([-15, 15])
 plt.title('XSp ampl')
phases=[]
plt.figure() 
for rx in range(noRx): # rx=0
 plt.subplot(1,3,rx+1)
 phases=np.angle(XSpectr[:,:,rx])/np.pi*180
 plt.pcolor(f,ranges,phases,cmap='jet')
 plt.clim([-180, 180])
 # plt.colorbar()
 plt.title('XSp phase')
 plt.figure()
for rx in range(noRx):
 plt.subplot(1,3,rx+1)
 ampl=10*np.log10(abs(XSpectr[:,:,rx])/2)
 SNRsel=ampl<-11
 ampl[SNRsel]="nan"
 plt.pcolor(f,ranges,ampl,cmap='jet')
 plt.clim([-10, 20])
 plt.title('XSp ampl')
plt.figure() 
for rx in range(noRx): # rx=0
 plt.subplot(1,3,rx+1)
 phases=np.angle(XSpectr[:,:,rx])/np.pi*180
 SNRsel=10*np.log10(abs(XSpectr[:,:,rx])/2)<-11 
 # phases[SNRsel]=float("nan")
 phases[SNRsel]="nan"
 plt.pcolor(f,ranges,phases,cmap='jet')
 plt.clim([-180, 180])
 # plt.colorbar()
 plt.title('XSp phase')
 plt.xlim([-.75,.75])
# cross-correlation for one range and two receivers -> testing reasons
RG=17
x1=data[RG,:,0]
x2=data[RG,:,1]
# t=np.linspace(0,15,16)
# x1=[1,0,0,1,1,1,0,0,0,1,1,0,0,1,0,0,]
# x2=np.round(np.random.rand(16))
plt.figure()
plt.plot(t,np.abs(x1),t,np.abs(x2))
# xcor=signal.correlate(x1s,x2s,method='direct',mode='full')
xcor=signal.correlate(x1,x2,method='direct',mode='full')
# xcor2=np.correlate(a=x1,v=x2)
# xcor3=np.cov(x1,x2)
#
xcor4=signal.correlate(np.real(x1),np.real(x2))+1j*signal.correlate(np.imag(x1),
np.imag(x2))
# xcor=xcor/max(xcor)
tc=(t-min(t))*60*60
t2=-1*tc[::-1]
t2=np.append(t2, tc[1:])
plt.figure()
plt.subplot(2,1,1)
plt.plot(t2,abs(xcor))
plt.xlabel('tau')
plt.ylabel('abs(xcor)')
plt.subplot(2,1,2)
plt.plot(t2,np.angle(xcor)/np.pi*180)
plt.xlabel('tau')
plt.ylabel('phase /°')
# cross-correlation for all ranges and all receivers
Xcor=np.zeros([noRG,noDP*2-1,noRx])+1j*np.zeros([noRG,noDP*2-1,noRx])
for rg in range(noRG):
 Xcor[rg,:,0]=signal.correlate(data[rg,:,0],data[rg,:,1])
 Xcor[rg,:,1]=signal.correlate(data[rg,:,0],data[rg,:,2])
 Xcor[rg,:,2]=signal.correlate(data[rg,:,1],data[rg,:,2])
plt.figure()
for rx in range(noRx):
 plt.subplot(2,3,rx+1)
 plt.pcolor(t2,ranges,10*np.log10(np.abs(Xcor[:,:,rx])),cmap='jet')
 plt.xlabel('abs(xcor) /dB')
 plt.ylabel('range /km')
 plt.clim([20,90])
 plt.colorbar()
for rx in range(noRx):
 plt.subplot(2,3,rx+4)
 plt.pcolor(t2,ranges,np.angle(Xcor[:,:,rx])/np.pi*180,cmap='jet')
 plt.clim([-180, 180])
 plt.colorbar()
 plt.ylabel('range /km')
 plt.xlabel('phase /°')
# line plots for the mean of xcor (along the data points - time )
plt.figure()
plt.subplot(1,3,1)
plt.plot(PWR,ranges)
plt.legend(['Rx1','Rx2','Rx3'])
plt.xlabel('power /dB')
plt.grid('on')
plt.subplot(1,3,2)
plt.plot(10*np.log10(Xcor[:,noDP,:]),ranges)
plt.legend(['1-2','1-3','2-3'])
plt.xlabel('xcor ampl. /dB')
plt.grid('on')
plt.subplot(1,3,3)
plt.plot(np.angle(Xcor[:,noDP,:])/np.pi*180,ranges)
plt.legend(['1-2','1-3','2-3'])
plt.xlabel('phase / °')
plt.grid('on')
# coherence for all ranges and all receivers
npts=256
cohtemp=np.zeros(256)
Coh=np.zeros([noRG,npts,noRx])
for rg in range(noRG):
 f,cohtemp=signal.coherence(data[rg,:,0],data[rg,:,1],Fs)
 Coh[rg,:,0]=np.fft.fftshift(cohtemp)
 f,cohtemp=signal.coherence(data[rg,:,0],data[rg,:,2],Fs)
 Coh[rg,:,1]=np.fft.fftshift(cohtemp)
 f,cohtemp=signal.coherence(data[rg,:,1],data[rg,:,2],Fs)
 Coh[rg,:,2]=np.fft.fftshift(cohtemp)
f=np.fft.fftshift(f)
plt.figure()
for rx in range(noRx):
 plt.subplot(1,3,rx+1)
 plt.pcolor(f,ranges,(np.abs(Coh[:,:,rx])))
 plt.xlabel('f / Hz')
 plt.ylabel('range /km')
 plt.clim([0.2,1])
 plt.xlim([-1,1])
 plt.colorbar()
plt.subplot(1,3,1)
plt.title('coherence RX1-2')
plt.subplot(1,3,2)
plt.title('coherence RX1-3')
plt.subplot(1,3,3)
plt.title('coherence RX2-3')

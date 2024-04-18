
+"""
DSP HomeWork 6
Due April 17th, 2024
By: Alexander Olson
"""
import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import time
from matplotlib.backends.backend_pdf import PdfPages
#%% Functions
# Ideal Low pass Filter
def IdealFilD (n,Oc):
    X = np.zeros(len(n))
    for k in range(0,len(X)):
        if k == len(X)/2:
            X[k] = 1*(Oc/np.pi)
        else:
            X[k] = (Oc/(np.pi))*((np.sin(Oc*n[k])/(Oc*n[k]))) # Funciton for ideal filter in discert time
    return X
def IdealFilF (f,Oc ):
    """
    Summary: Is the Ideal function for a Low Pass filter in the Frequncey Domain
    INPUTS: f = current frequncy, Oc = Cutoff or Frequncy range
    OUPUTS: Magnitude response
    """
    if(abs(f)<= Oc): 
        return 1
    else: return 1*np.exp(-100);
                             
def fouierTransform(Np, O, x):
    """
    Summary: Summation Defintion of the FouierTransform
    INPUTS: Np = Number of Samples, O = range of Omega Values, x = discert signal
    OUTPUT: X = sigle in frequncy domain
    """
    X= 0.0
    for k in range(0, int(Np), 1): # perform summation for definition
      X = X + x[k]*np.exp(-1j*O*k)
    return X
# Windowing funciton
def win(n,c):
    # n = number of samples, c = m+1 where m is the number of points from the center of the window
    if( abs(n)> c):
        return 0
    return 0.54+0.46*np.cos((n*np.pi)/c)
def convolve(x,h):
    # Function to convolve two functions, x is an array of values, h is another array of values that are from a window function,
    # M is the number of values on each side of the origin
    M = len(x)
    N = len(h)
    # Length of the result
    L = M + N - 1
    # Initialize result with zeros
    result = [0] * L
    
    # Compute convolution sum
    for n in range(L):
        for k in range(max(0, n-N+1), min(M, n+1)):
            result[n] += x[k] * h[n - k]
    return result      
#%% Given/detrimined Conditions/Data Sets
fs = 20000 # Sampling Frequncy, Hz
ftones = np.array([697, 770, 852, 941, 1209, 1336, 1477, 1633]) # Given Tones in Hz
OmegaC = ftones[0]+((ftones[1] - ftones[0])/2) # Cutoff frequncy, Hz
Omega = np.linspace(-np.pi, np.pi, fs+1) # Range of Omega/Frequncy Values in radian
tf = 5.0 # duration of sound in seconds
Np = fs*tf # Number of Samples
nVecx = np.arange(Np) # Vector of discrete-times for input audio signal x[n]
inputSig = np.zeros(int(tf * fs)) # Initialize vector for audio signal x[n]
#Input Signal
for m in range(0, len(ftones), 1): # Loop over tones to build x[n]
  inputSig = inputSig + np.sin(2 * np.pi * nVecx * ftones[m] / fs)
inSigFT = fouierTransform(Np , Omega, inputSig) # Fouier Transfrom of the input signal
# Ideal Low pass filterin frequncy domain
X = np.zeros(len(Omega))
for i in range(0,len(X)):
    X[i] = IdealFilF(Omega[i],(OmegaC*2*np.pi)/fs)
filFT = inSigFT *X # Filtered Signal in the Frqency Domain
hIdealN = IdealFilD(nVecx-(len(nVecx)/2),OmegaC*2*np.pi/fs) #impluse function of the ideal LPF in discert time
#%% Generating the Window Function/Paramters
m = 300 # Window Size
winN = np.zeros(len(nVecx))
for j in (nVecx):
    winN[int(j)] = win(j-(len(nVecx)/2), m +1 ) # Finding the window funcution in discert time
# Fouier Transfrom of the Window funciton
WinMod = winN/sum(winN) # Normilizeing the Window
winF = fouierTransform(Np, Omega, WinMod)
#%% Creating the Practical Filter and Appyling to the Ideal Filter
#Pracitcal Filter
hPracN = hIdealN*winN
hPracF = fouierTransform(Np , Omega, hPracN)
# Applying Filter to create output signal
yN = np.zeros(len(inputSig))
yNtest = np.convolve(inputSig,hPracN ) # Using inbuilt function
nVecC = (np.arange(len(yNtest)))
# Using self made convolution
hPracNMod = []
for n in range(len(hPracN)):
    if hPracN[n] != 0:
        hPracNMod.append(hPracN[n])
yNReal = convolve(inputSig, hPracNMod)
nVecMod = np.arange(len(yNReal))
# Generating Refernce Signal
refSig = np.sin(2 * np.pi * nVecx * 697 / fs)
# All of there fourier transforms
yNRealF = fouierTransform(Np, Omega,yNReal )
refSigF = fouierTransform(Np, Omega,refSig )
#%% Part A Plotting
# Formatting
figA = plt.figure(figsize=(12, 16))
Ga = plt.GridSpec(3, 2)
# Plotting Ideal LPF in discert time

IdealD = plt.subplot(Ga[0,:])
plt.stem(nVecx-(len(nVecx)/2),hIdealN )
plt.xlim([-500,500])
plt.title('Ideal Low Pass Filter (Discert Time)')
plt.xlabel('n (Number of Samples)')
plt.ylabel('h[n]')
plt.grid()
# Plotting Ideal LPF in Frequncy Domain
#Magnitude Response
IdealF = plt.subplot(Ga[1,:])
plt.plot(Omega, X)
RadOmegaC = "%.3f"%(OmegaC*2*np.pi/fs)
plt.scatter((OmegaC*2*np.pi)/fs, 0, c = 'k', label = f'Cutoff Frequncy {OmegaC} Hz | {RadOmegaC} Radians')
plt.scatter(-((OmegaC*2*np.pi)/fs), 0, c = 'k')
plt.xlim([-np.pi-.1, np.pi+.1])
plt.title('Ideal Low Pass Filter Magnitude Response')
plt.xlabel('Ω (Radians)')
plt.ylabel('|H(Ω)|')
plt.legend()
plt.grid()
# Plotting Ideal LPF in Frequncy Domain
#Magnitude Response
IdealFdB = plt.subplot(Ga[2,:])
plt.plot(Omega, 20*np.log10(X))
plt.scatter((OmegaC*2*np.pi)/fs, 0, c = 'k', label = f'Cutoff Frequncy {OmegaC} Hz | {RadOmegaC} Radians')
plt.scatter(-((OmegaC*2*np.pi)/fs), 0, c = 'k')
plt.xlim([-np.pi-.1, np.pi+.1])
plt.title('Ideal Low Pass Filter Magnitude Response in dB')
plt.legend()
plt.xlabel('Ω (Radians)')
plt.ylabel('|H(Ω)|')
plt.grid()
#Frequncy Response

#Part B Plotting 

# Formatting/ Plotting
figB= plt.figure(figsize=(12, 16))
GB = plt.GridSpec(3, 2)
WinD = plt.subplot(Ga[0,:])
plt.stem(nVecx-(len(nVecx)/2), winN, label = f'0.54+0.46cos((nπ)/c): Where c is equal to {m+1}')
plt.xlim([-310,310])
plt.title('Hamming Window Funciton (Discert Time)')
plt.xlabel('n (Number of Samples)')
plt.ylabel('w[n]')
plt.grid()
plt.legend()
# Plotting Window Funciton in Frequncy Domain
#Magnitude Response
Winf = plt.subplot(Ga[1,:])
plt.plot(Omega, abs(winF))
plt.xlim([-.2, .2])
plt.title(' Hamming Window Magnitude Response')
plt.xlabel('Ω (Radians)')
plt.ylabel('|W(Ω)|')
plt.grid()
# Plotting Window Funciton in Frequncy Domain
#Magnitude Response in dB
WinFdB = plt.subplot(Ga[2,:])
plt.plot(Omega, 20*np.log10(abs(winF)))
plt.xlim([-.5, .5])
plt.title('Hamming Window Magnitude Response in dB')
plt.xlabel('Ω (Radians)')
plt.ylabel('|W(Ω)|')
plt.grid()

# Part C Plotting
figC= plt.figure(figsize=(12, 16))
HPracN = plt.subplot(Ga[0,:])
plt.stem(nVecx-(len(nVecx)/2), hPracN)
plt.title('Pratical Low Pass Filter (Discert Time)')
plt.xlabel('n (Number of Samples)')
plt.ylabel('h[n]')
plt.xlim([-500, 500])
plt.grid()
hPraCF = plt.subplot(Ga[1,:])
plt.plot(Omega, abs(hPracF))
plt.scatter((OmegaC*2*np.pi)/fs, 0, c = 'k', label = f'Cutoff Frequncy {OmegaC} Hz | {RadOmegaC} Radians')
plt.scatter(-((OmegaC*2*np.pi)/fs), 0, c = 'k')
plt.scatter(-((770*2*np.pi)/fs), 0, label = ' 770 hz')
plt.scatter(-((852*2*np.pi)/fs), 0,label = ' 852 hz')
plt.scatter(-((941*2*np.pi)/fs), 0,label = ' 941 hz')
plt.xlim([-.5, .5])
plt.title('Pratical Low Pass Filter (Frequncy Domain)')
plt.xlabel('Ω (Radians)')
plt.ylabel('|W(Ω)|')
plt.legend()
plt.grid()
hPraCF = plt.subplot(Ga[2,:])
plt.plot(Omega, 20*np.log10(abs(hPracF)))
plt.scatter((OmegaC*2*np.pi)/fs, 0, c = 'k', label = f'Cutoff Frequncy {OmegaC} Hz | {RadOmegaC} Radians')
plt.scatter(-((OmegaC*2*np.pi)/fs), 0, c = 'k')
plt.title('Pratical Low Pass Filter in dB (Frequncy Domain)')
plt.xlabel('Ω (Radians)')
plt.ylabel('|W(Ω)|')
plt.grid()

# Part D Plotting
figC= plt.figure(figsize=(12, 16))
gd = plt.GridSpec(4, 1)
yNplot = plt.subplot(Ga[0,:])
plt.plot(nVecx-(len(nVecx)/2),refSig,label = "Refernce Signal")
plt.plot(nVecMod-(len(nVecMod)/2),yNReal,label = "Filter Applyed Using conolve (Self Made)")
plt.plot(nVecC-(len(nVecC)/2),yNtest, label = "Filter Applyed Using Np.Convolve")
plt.xlim([-50,50])
plt.title('Filtered Signal (Discert Time and when in Steady State) ')
plt.xlabel('n (Number of Samples)')
plt.ylabel('y[n]')
plt.legend()
plt.grid()
yNplot = plt.subplot(Ga[1,:])
plt.plot(nVecx-(len(nVecx)/2),refSig,label = "Refernce Signal")
plt.plot(nVecC-(len(nVecC)/2),yNtest, label = "Filter Applyed Using Np.Convolve")
plt.plot(nVecMod-(len(nVecMod)/2),yNReal,linestyle = '--',label = "Filter Applyed Using conolve (Self Made)")
plt.xlim([-50250,-49750])
plt.title('Filtered Signal (Discert Time and showing edge cases) ')
plt.xlabel('n (Number of Samples)')
plt.ylabel('y[n]')
plt.legend()
plt.grid()

yNplot = plt.subplot(Ga[2,:])
plt.plot(Omega,abs(refSigF), label = 'Reference Signal',c= 'k', linewidth = 3)
plt.plot(Omega,abs(yNRealF), label = 'convolve (self made)')
plt.plot(Omega,abs(filFT) ,linestyle = '--',label = 'Np.convolve')
plt.title('Filtered Signal (Frequncy Domain)')
plt.xlabel('Ω (Radians)')
plt.ylabel('|W(Ω)|')
plt.xlim([0,.3])
plt.legend()
plt.grid()
#%% E Playing Audio
# Play the input audio sigal through speakers
'''
print('Playing audio signal ...')
sd.play(refSig, fs) # start playing audio signal
time.sleep(tf) # stall script while audio is playing
sd.stop() # stop playing audio signal
sd.play(yNReal, fs) # start playing audio signal
time.sleep(tf) # stall script while audio is playing
sd.stop() # stop playing audio signal
sd.play(yNtest, fs) # start playing audio signal
time.sleep(tf) # stall script while audio is playing
sd.stop() # stop playing audio signal
'''
#%% Printing Figures
p = PdfPages("DSPHw6Graphs.pdf")
# get_fignums Return list of existing  
# figure numbers 
fig_nums = plt.get_fignums()   
figs = [plt.figure(n) for n in fig_nums]     
# iterating over the numbers in list 
for fig in figs:       
   # and saving the files 
   fig.savefig(p, format='pdf')  
p.close()
"""Toutes les fonctions utiles pour les TPs de R205 !"""

from sounddevice import *
from pylab import *
import scipy.signal as sc       # pour utiliser les fonctions liées aux filtres 

def plotFreq(sig,fe,c,titre):
    """Trace le spectre d'amplitude du signal sig, qui est échantillonné à fe Hz. c indique la couleur souhaité (ex : 'b' pour bleu). titre indique le titre."""
    longueur=len(sig)
    specOrig=abs((1/longueur)*fftshift(fft(sig)))
    freq = arange(-fe/2,fe/2,fe/longueur) 
    style = c+'-'
    #display(style)
    plot(freq, specOrig,style)
    grid()
    title(titre, fontsize=16)
    xlim([-fe/2, fe/2])
    xlabel('fréquence (Hz)', fontsize=16)
    ylabel('Amplitude', fontsize=16)
    
def plotTimeFreq(sig,fe):
    """Effectue une analyse temps-fréquence du signal sig, échantillonné à fe Hz."""
    specgram(sig,Fs=fe,NFFT = 4096)
#     title('Analyse temps fréquence', fontsize=16)
    xlabel('Time (s)', fontsize=16)
    ylabel('Frequency (Hz)', fontsize=16)

def filtrage(s, Fs, Fc, typ):
    """applique un filtre passe bas ou passe-haut de fréquence de coupure Fc sur le signal s
    - s : signal que l'on souhaite filtrer
    - Fs : fréquence d'échantillonnage du signal
    - Fc : fréquence de coupure du filtre
    - typ : 'low' ou 'high' , respectivement passe-bas ou passe-haut 
    - Exemple : sfiltre = filtrage(son, 44000, 2000,'low')
    - Cette instruction va filtrer le signal "son" (intialement échantillonné à 44000Hz) 
    à l'aide d'un passe-bas de fréquence de coupure 2000Hz. 
    Le signal filtré sera stocké dans la variable sfiltre.
    """
    ordre=8
    b, a = sc.butter(ordre, Fc/(Fs/2), typ)
    return sc.lfilter(b, a, s)


def puissancedBm(signal):
    """Fonction permet de calculer la puissance en dBm de signal.
    ===========================
    Parametre : 
    - signal : signal a transmettre
    ===========================
    """
    Pref = 0.001
    puissanceW = linalg.norm(signal,2)**2/size(signal)
    puissance = 10*math.log(puissanceW/Pref,10)
    return round(puissance,4)

import math
def channel(signal, Type, varargin_1,varargin_2,varargin_3=None) :
    """La fonction signalRecu = channel(signal, Type, varargin) permet de simuler l'attenuation de differents types de canaux transmission.
    Le nombre et la nature des parametres changent en fonction du type du canal choisi.
    ===========================
    Parametres d'entree :
    - signal : signal a transmettre
    - Type : type de canal : 'espacelibre', 'coaxial' et 'fibre'
    - varargin : 
    si Type = 'espacelibre'
     varargin_1 : longueur du canal (en km)
     varargin_2 : frequence porteuse (en MHz)
    si Type = 'coaxial'
     varargin_1 : longueur du canal (en km)
     varargin_2 : pertes d'insertion (en dB/100m)
    si Type = 'fibre'
     varargin_1 : longueur du canal (en km)
     varargin_2 : affaiblissement (en dB/km)
     varargin_3 : longueur d'onde pour la transmission (en nm)
    Parametres de sortie :
    - signalRecu : signal en sortie du canal
    """
    if Type == 'espacelibre' :
        alpha = 2
        pertes_dBtot = 32.45 + 20*math.log(varargin_2,10) + 10*alpha*math.log(varargin_1,10)
        
    elif Type == 'coaxial' : 
        distance_km = varargin_1
        pertes_dBper100metre = varargin_2
        pertes_dBtot = pertes_dBper100metre*10*distance_km
        
    elif Type == 'fibre' :
        distance_km = varargin_1
        pertes_dBperkm = varargin_2
        L0_nm = varargin_3
        c = 3*10**8
        Fp_Hz = c/L0_nm
        pertes_dBtot = pertes_dBperkm*distance_km
        
    gainG  = 10**(-pertes_dBtot/20)
    signalRecu = []
    for i in signal :
        signalRecu.append(i * gainG)
    return signalRecu
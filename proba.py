# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 23:19:07 2023

@author: Nikola
"""
import numpy as np
from math import *
from argparse import ArgumentParser

class Transformacje:
    
    # TU MI COS NIE SIEDZI 
    def __init__(self,model):
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Elipsoida Krasowskiego":
            self.a = 6370245.0
            self.b = 6356863.01877
            
        self.sp = (self.a - self.b) / self.a
        self.e2 = (2 * self.sp - self.sp ** 2) 
        
        
    def Np(self,fi): #promień krzywizny w 1 wektykale
        wynik =[]   
        for f in zip(fi):    
            N=self.a/np.sqrt(1-self.e2*np.sin(f)**2)
            wynik.append(N)
        return wynik
    
    def przetworzenie(self,plik, metoda):
        dane = np.genfromtext(plik)
        if metoda == "Np":
            wynik = self.Np(dane[:0])
            np.savetxt("wyniki",delimiter= ";")
        return wynik
            
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-plik" , type = str, help = "sciezka do pliku")
    parser.add_argument("-trans", type = str, help = "wybrana transformacja")
    args = parser.parse_args()
    
    model = {"wgs84":"wgs84"} # tu jakies parametry trzeba wstawic
    trans = {"Npu":"Npu"}
    
    try:
        args.plik = None
        args.trans = None
        args.plik = "dane.txt"
        args.trans = "Np"
        obiekt = Transformacje(model)
        dane = obiekt.przetworzenie(args.p,trans[args.trans])
    finally:
        print("Plik wynikowy zapisany.")
    print("proba")

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
        

    
    def hirvonen(self,X,Y,Z): #XYZ zamieniamy na flh
        p = np.sqrt(X**2+Y**2)
        fi = np.arctan(Z/(p*(1-self.e2)))
        while True: #pętla
            N=self.a/np.sqrt(1-self.e2*np.sin(fi)**2)
            h=p/np.cos(fi)-N
            fip=fi
            fi=np.arctan(Z/(p*(1-self.e2*N/(N+h))))
            if abs(fip-fi)<(0.000001/206265):
                break
        l=np.arctan2(Y,X) #lambda
        return(fi,l,h)
    
    def flh2XYZ(self,fi,l,h):
        while True:
            N=self.a/np.sqrt(1-self.e2*np.sin(fi)**2)
            X=(N+h)*np.cos(fi)*np.cos(l)
            Xp=X
            Y=(N+h)*np.cos(fi)*np.sin(l)
            Z=(N*(1-self.e2)+h)*np.sin(fi)
            if abs(Xp-X)<(0.000001/206265):
                break
        return(X,Y,Z)

    def NEU(self,fi,l,v):
        """
        Funckja obliczająca wektor w układzie neu
    
        Parameters:
        -----------
        R : R : [array of float64] : macierz obrotu
        v : [array of float64] : wektor w układzie XYZ
        
        Returns:
        -------
        NEU : [array of float64] : współrzedne topocentryczne (North , East (E), Up (U))
    
        """
        N=[(-np.sin(fi) * np.cos(l)), (-np.sin(fi) * np.sin(l)), (np.cos(fi))]
        E=[(-np.sin(l)), (np.cos(l)),  (0)]
        U=[( np.cos(fi) * np.cos(l)), ( np.cos(fi) * np.sin(l)), (np.sin(fi))]
        R=np.transpose(np.array([N,E,U]))
        NEU=np.zeros(v.shape)
        for a in range(v.shape[0]):
            for b in range(3):
                for c in range(3):
                    NEU[a,c]+=v[a,b]*R[c,b]
        return (NEU)

    
if __name__ == "__main__":
    # def Mp(fi):
    #     M=self.a*(1-self.e2)/np.sqrt((1-self.e2*(np.sin(fi))**2)**3) 
    #     return(M)
    
    parser = ArgumentParser()
    parser.add_argument("-plik" , type = str, help = "sciezka do pliku")
    parser.add_argument("-trans", type = str, help = "wybrana transformacja")
    parser.add_argument("-model" , type = str, help = "wybrany model")
    args = parser.parse_args()
    
    # model = {"wgs84":"wgs84"} # tu jakies parametry trzeba wstawic
    model = {"wgs84":"wgs84", "grs80":"grs80", "Elipsoida Krasowskiego":"Elipsoida Krasowskiego"}
    trans = {"hirvonen": "hirvonen", "flh2XYZ": "flh2XYZ"}
    
    try:
        dane = np.genfromtxt(args.plik,delimiter=",")
        obiekt = Transformacje(model[args.model])
        a = obiekt.a
        e2 = obiekt.e2
        result = []
        print(dane)
        for xyz in dane:    
            if trans[args.trans]=="hirvonen":
                line = obiekt.hirvonen(xyz[0],xyz[1],xyz[2])
                result.append(line)
            elif trans[args.trans]=="flh2XYZ":
                print(xyz)
                line = obiekt.flh2XYZ(xyz[0],xyz[1],xyz[2])
                result.append(line)
        # elif trans[args.trans] == "Mp":
        #     # print(obiekt.Np(dane[:,0]))
        #     wynik = obiekt.Mp(dane[:,0])
        print(result)
        np.savetxt("wyniki.txt",result,delimiter=",")
            
        
    finally:
        print("Plik wynikowy zapisany.")

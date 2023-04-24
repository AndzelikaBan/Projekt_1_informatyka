# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 01:08:31 2023

@author: andze
"""
import numpy as np
from math import *
class Calculations:
   

    def Np(self,fi,a,e2): #promień krzywizny w 1 wektykale w punkcie P
        N=a/np.sqrt(1-e2*np.sin(fi)**2)
        return(N)

    def Mp(self,fi,a,e2):
        M=a*(1-e2)/np.sqrt((1-e2*(np.sin(fi))**2)**3) 
        return(M)

    def Hirvonen(self,X,Y,Z,a,e2): #XYZ zamieniamy na flh
        p = np.sqrt(X**2+Y**2)
        fi = np.arctan(Z/(p*(1-e2)))
        while True: #pętla
            N=self.Np(fi,a,e2) 
            h=p/np.cos(fi)-N
            fip=fi
            fi=np.arctan(Z/(p*(1-e2*N/(N+h))))
            if abs(fip-fi)<(0.000001/206265):
                break
        l=np.arctan2(Y,X) #lambda
        return(fi,l,h)
    
    def dms(self,x):
        znak=' '
        if x<0:
            znak='-'
            x=abs(x) #wartosć dodatnia
        x=x*180/np.pi
        d=int(x) #stopnie
        m=int((x-d)*60) #minuty
        s=(x-d-m/60)*3600 #sekundy
        print(znak,"%3d⁰%2d'%7.5f''"%(d,m,s)) #d tutaj oznacza liczbę całkowitą, f liczbę rzeczywistą
        #7.5 oznacza 7 cyfr, 5 po przecinku
        
    def flh2XYZ(self,fi,l,h,a,e2):
        while True:
            N=self.Np(fi,a,e2)
            X=(N+h)*np.cos(fi)*np.cos(l)
            Xp=X
            Y=(N+h)*np.cos(fi)*np.sin(l)
            Z=(N*(1-e2)+h)*np.sin(fi)
            if abs(Xp-X)<(0.000001/206265):
                break
        return(X,Y,Z)
    def Transformacje(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Elipsoida Krasowskiego":
            self.a = 6370245.0
            self.b = 6356863.01877
        else:
            raise NotImplementedError(f"{model} Ta operacja jest niemożliwa!")
        self.sp = (self.a - self.b) / self.a
        self.e2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
        
        
        

    

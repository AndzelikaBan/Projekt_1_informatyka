# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 01:08:31 2023

@author: andze
"""
import numpy as np
from math import *
import argparse 
class Calculations:
   
    def __init__(self, model: str = "wgs84"):
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
        self.e2 = (2 * self.sp - self.sp ** 2) # eccentricity**2
        
    def Np(self,fi): #promień krzywizny w 1 wektykale w punkcie P
        N=self.a/np.sqrt(1-self.e2*np.sin(fi)**2)
        return(N)

    def Mp(self,fi):
        M=self.a*(1-self.e2)/np.sqrt((1-self.e2*(np.sin(fi))**2)**3) 
        return(M)

    def Hirvonen(self,X,Y,Z): #XYZ zamieniamy na flh
        p = np.sqrt(X**2+Y**2)
        fi = np.arctan(Z/(p*(1-self.e2)))
        while True: #pętla
            N=self.Np(fi) 
            h=p/np.cos(fi)-N
            fip=fi
            fi=np.arctan(Z/(p*(1-self.e2*N/(N+h))))
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
        
    def flh2XYZ(self,fi,l,h):
        while True:
            N=self.Np(fi)
            X=(N+h)*np.cos(fi)*np.cos(l)
            Xp=X
            Y=(N+h)*np.cos(fi)*np.sin(l)
            Z=(N*(1-self.e2)+h)*np.sin(fi)
            if abs(Xp-X)<(0.000001/206265):
                break
        return(X,Y,Z)


    # def saz2neu1(self,s,alfa,z):
    #     dneu=np.array([s*np.sin(z)*np.cos(alfa),
    #                     s*np.sin(z)*np.sin(alfa),
    #                     s*np.cos(z)])
    #     return(dneu)



    # def XYZ2neu(self,dneu,fi,l):
    #     R=np.array([[-np.sin(fi)*np.cos(l),-np.sin(l),np.cos(fi)*np.cos(l)],
    #                 [-np.sin(fi)*np.sin(l),np.cos(l),np.cos(fi)*np.sin(l)],
    #                 [ np.cos(fi),            0.     ,np.sin(fi)]])
        
    #     return(R.T @ dneu)
    def Rneu(self, fi, l):
        """
        Funkcja, która, przyjmujac współrzedne krzywoliniowe utworzy macierz obrotu 
        potrzebną do przeliczenia współrzędnych do układu współrzędnych neu
    
        INPUT:
        ----------
        phi : [float] : wspołrzędna fi punktu początkowego układu lokalnego
        lam : [float] :wspołrzędna l punktu początkowego układu lokalnego
    
        OUTPUT:
        -------
        R : [array of float64] : macierz obrotu
    
        """
        N=[(-np.sin(fi) * np.cos(l)), (-np.sin(fi) * np.sin(l)), (np.cos(fi))]
        E=[(-np.sin(l)), (np.cos(l)),  (0)]
        U=[( np.cos(fi) * np.cos(l)), ( np.cos(fi) * np.sin(l)), (np.sin(fi))]
        R=np.transpose(np.array([N,E,U]))
        return (R, N, E, U)
    
    def NEU(self, R,v):
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
        NEU=np.zeros(v.shape)
        for a in range(v.shape[0]):
            for b in range(3):
                for c in range(3):
                    NEU[a,c]+=v[a,b]*R[c,b]
        return (NEU)
    def deg2rad(self, kat):
        '''
        Funkcja przelicza kąt w stopniach na radiany.
        
        Argumenty:
        ----------
        kat - kąt podany w stopniach | typ: float lub int
        
        Wyniki:
        -------
        katr - kąt przeliczony z argumentu 'kat' na radiany | typ: float
        
        Dodatkowy opis:
        ---------------
        '''
        katr = kat * pi / 180
        return(katr)
    
    def pl1992(self,fi,l,m=0.9993):
        """
        Przeniesienie wspolrzednych krzywoliniowych geodezyjnych punktu A
        do ukladu PL-1992

        Parameters
        ----------
        fi : float
            Szerokosc geodezyjna punktu A.
            Jednostka -- RAD
        l : float
            Dlugosc geodezyjna punktu A.
            Jednostka -- RAD
        a : float, optional
            Polos wielka. The default is 6378137.
            Jednostka -- METR
        e2 : float, optional
            I mimosrod elipsoidy. The default is 0.00669438002290.
            Jednostka -- brak
        m : float, optional
            Wspolczynnik zmiany skali. The default is 0.9993.
            Jednostka -- brak

        Returns
        -------
        x92 : float
            Wspolrzedna X punktu A w ukladzie PL-1992.
            Jednostka -- METR
        y92 : TYPE
            Wspolrzedna Y punktu A w ukladzie PL-1992.
            Jednostka -- METR

        """
        l0 = np.deg2rad(19)
        # 1 parametry elipsoidy     
        b2 = self.a**2*(1-self.e2)
        ep2 = (self.a**2-b2)/b2
        # 2. Wielkosci pomocnicze     
        dell = l - l0
        t = np.tan(fi)
        ni2 = ep2*(np.cos(fi)**2)
        N = self.Np(fi)
        
        # 3. Długosc luku poludnika 
        A0 = 1- (self.e2/4)-(3*self.e2**2/64)-(5*self.e2**3/256)
        A2 = (3/8)*(self.e2+(self.e2**2/4)+(15*self.e2**3/128))
        A4 = (15/256)*(self.e2**2+((3*self.e2**3)/4))
        A6 = (35*self.e2**3)/3072
        
        sigma = self.a *(A0*fi-A2*np.sin(2*fi)+A4*np.sin(4*fi)-A6*np.sin(6*fi))
        
        # wsolrzedne prostokatne lokalne na plaszczyznie gaussa-krugera
        
        xgk =  sigma    +    ( ((dell**2/2)*N*np.sin(fi)*np.cos(fi))    *    (1   +   ((dell**2/12)*(np.cos(fi)**2)*(5 - t**2 + 9*ni2 + 4*ni2**2))      +         ((dell**4/360)*(np.cos(fi)**4)*(61 - 58*t**2 + t**4 + 270*ni2 - 330*ni2*t**2))))
        
        ygk =  (dell* N * np.cos(fi))  *   ( 1 +  ((dell**2/6)   *   (np.cos(fi)**2)   *  (1 - t**2 + ni2))     +     (((dell**4/120)*(np.cos(fi)**4)) * (5 - (18*t**2) + t**4 + (14 * ni2) - (58*ni2*t**2))))
        
        x92 = xgk * m - 5300000
        y92 = ygk*m + 500000
        return  x92, y92,xgk,ygk
    def pl2000(self,fi,l,m=0.999923):
        """
        Przeniesienie wspolrzednych geodezyjnych krzywoliniowych fi lam punktu A
        do ukladu wspolrzednych PL-2000

        Parameters
        ----------
        fi : float
            Szerokosc geodezujna punktu A.
            Jednostka -- RAD
        l : float
            Dlugosc geodeztjna punktu A.
            Jednostka -- RAD
        a : float, optional
            Wielka połos. The default is 6378137.
            Jednostka -- METR
        e2 : float, optional
            I mimosrod elipsoidy. The default is 0.00669438002290.
            Jednostka -- brak
        m : float, optional
            Wspolczynnik zmiany skali . The default is 0.999923.
            Jednostka -- bral

        Returns
        -------
        x2000 : float
            Wspolrzedna X punktu A w ukladzie wspolrzednych PL-2000.
            Jednostka -- METR
        y2000 : float
            Wspolrzedna X punktu A w ukladzie wspolrzednych PL-2000.
            Jednostka -- METR

        """
        l0=0 
        strefa = 0
        if l >np.deg2rad(13.5) and l < np.deg2rad(16.5):
            strefa = 5
            l0 = np.deg2rad(15)
        elif l >np.deg2rad(16.5) and l < np.deg2rad(19.5):
            strefa = 6
            l = np.deg2rad(18)
        elif l >np.deg2rad(19.5) and l < np.deg2rad(22.5):
            strefa =7
            l0 = np.deg2rad(21)
        elif l >np.deg2rad(22.5) and l < np.deg2rad(25.5):
            strefa = 8
            l0 = np.deg2rad(24)
        else:
            print("Punkt poza strefami odwzorowawczymi układu PL-2000")        
        
        # 1 parametry elipsoidy     
        b2 = self.a**2*(1-self.e2)
        ep2 = (self.a**2-b2)/b2
        # 2. Wielkosci pomocnicze     
        dell = l - l0
        t = np.tan(fi)
        ni2 = ep2*(np.cos(fi)**2)
        N = self.Np(fi)
        
        # 3. Długosc luku poludnika 
        A0 = 1- (self.e2/4)-(3*self.e2**2/64)-(5*self.e2**3/256)
        A2 = (3/8)*(self.e2+(self.e2**2/4)+(15*self.e2**3/128))
        A4 = (15/256)*(self.e2**2+((3*self.e2**3)/4))
        A6 = (35*self.e2**3)/3072
        
        sigma = self.a *(A0*fi-A2*np.sin(2*fi)+A4*np.sin(4*fi)-A6*np.sin(6*fi))
        
        # wsolrzedne prostokatne lokalne na plaszczyznie gaussa-krugera
        
        xgk =  sigma    +    ( ((dell**2/2)*N*np.sin(fi)*np.cos(fi))    *    (1   +   ((dell**2/12)*(np.cos(fi)**2)*(5 - t**2 + 9*ni2 + 4*ni2**2))      +         ((dell**4/360)*(np.cos(fi)**4)*(61 - 58*t**2 + t**4 + 270*ni2 - 330*ni2*t**2))))
        
        ygk =  (dell* N * np.cos(fi))  *   ( 1 +  ((dell**2/6)   *   (np.cos(fi)**2)   *  (1 - t**2 + ni2))     +     (((dell**4/120)*(np.cos(fi)**4)) * (5 - (18*t**2) + t**4 + (14 * ni2) - (58*ni2*t**2))))
        
        x2000 = xgk * m 
        y2000 = ygk*m + (strefa *1000000) +500000
        return  x2000, y2000,xgk,ygk
    
    
    
if __name__ == "__main__":
#utworzenie obiektu 
    parser=argparse.ArgumentParser()
    parser.add_argument("file_name", type = str,help="Nazwa pliku wejsciowego")
    args=parser.parse_args()
    

    # X = 23
    # Y= 89
    # Z = 67

    calc = Calculations() 
    tablica = np.genfromtxt(args.file_name, delimiter=',', skip_header = 4)
    # with open("wyniki.txt","a") as plik_w:
    #     plik_w.write(str(calc.Hirvonen(X, Y, Z)))
    print(tablica)    

    wybor = ("wybierz co chcesz zrobić:\n 1 Hirvonen \n Neu")
        
        
        

    

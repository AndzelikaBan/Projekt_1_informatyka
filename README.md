# PROJEKT 1 - TRANSFORMACJE

OPIS DZIAŁANIA PROGRAMU: 
  Program służy do przeliczania współrzędnych między różnymi układami na różnych elipsoidach odniesienia 
  (WGS84, GRS80, elipsoida Krasowskiego). Mamy 4 opcje przeliczania współrzędnych:
   * XYZ -> FLH
   * FLH -> XYZ
   * FL -> NEU
   * FL -> PL2000
   * FL -> PL1992
    
WYMAGANIA:
  * python 3.9
  * biblioteki: numpy, math, argparse
  * obsługuje system Windows 10 i Windows 11

JAK KORZYSTAĆ?:
W celu poprawnego korzystania z programu konieczne będzie utowrzenie pliku ze współrzędnymi (.txt).
Współrzędne muszą być w kolejności: X Y Z - najlepiej wyrażone w metrach. Jeden wiersz odpowiada współrzędnym
jednego punktu i jego dane powinny być oddzielone "," (przecinkiem). Części dziesiętne muszą znajdować się po
"." (kropce). Poniżej przykładowe dane:

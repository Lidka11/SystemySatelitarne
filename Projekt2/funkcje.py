from readrnx_studenci import*
import numpy as np
import math as mat

def Rneu(phi_rad, lam_rad):
    macierz = np.array([[-np.sin(phi_rad) * np.cos(lam_rad), -np.sin(lam_rad), np.cos(phi_rad)*np.cos(lam_rad)],
                        [-np.sin(phi_rad) * np.sin(lam_rad), np.cos(lam_rad), np.cos(phi_rad) * np.sin(lam_rad)],
                        [np.cos(phi_rad), 0, np.sin(phi_rad)]]) 
    return macierz

def xyz2blh(X, Y, Z):
    a = 6378137
    e2 = 0.00669438002290
    p = np.sqrt(X**2 + Y**2)
    phi = np.arctan2(Z, p * (1 - e2))
    phi_prev = 2 * np.pi
    while abs(phi - phi_prev) > 1e-10:
        N = a / (np.sqrt(1 - e2 * np.sin(phi) * np.sin(phi)))
        h = p / np.cos(phi) - N
        phi_prev = phi
        phi = np.arctan2(Z, p * (1 - e2 * N / (N + h)))
    lam = np.arctan2(Y, X)
    h = p / np.cos(phi) - N
    phi = np.rad2deg(phi)
    lam = np.rad2deg(lam)
    return phi, lam, h

def Np(B):
    a = 6378137
    e2 = 0.00669438002290
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N

def hirvonen(X,Y,Z):
    e2 = 0.00669438002290
    r = (X**2 + Y**2)**0.5
    B = mat.atan(Z/(r*(1-e2)))
    
    while 1:
        N = Np(B)
        H = r/np.cos(B) - N
        Bst = B
        B = mat.atan(Z/(r*(1-(e2*(N/(N+H))))))    
        if abs(Bst-B)<(0.00001/206265):
            break
    L = mat.atan2(Y,X)
    N = Np(B)
    H = r/np.cos(B) - N
    return B, L, H 


def obl_wspolrzednych_satelity(sow, week, sat, inav, nav):
    index_satelity = inav == sat
    nav_satelity = nav[index_satelity, :]
    toe = nav_satelity[:, 17]
    roznica = np.abs(sow - toe)
    index_najmniejszej_roznicy = np.argmin(roznica)
    wiersz_nav = nav_satelity[index_najmniejszej_roznicy, :]
    e = wiersz_nav[14] #ekscentr (mimośród) orbity
    pierwiastek_a = wiersz_nav[16] #pierwiastek z dużej półosi orbity
    Omega_0 = wiersz_nav[19] #rektascenzja (długość geograficzna) węzła wstępującego na początek tygodnia GPS
    perigee = wiersz_nav[23] #argument perygeum
    M0 = wiersz_nav[12] #anomalia średnia na epokę odniesienia
    Toe = wiersz_nav[17] #SOW epoka wyznaczenia efemerydy, dana w sekundach tygodnia GPS (epokaodniesienia efemerydy)
    i0 = wiersz_nav[21] #kąt inklinacji na epokę odniesienia
    dt = (wiersz_nav[9]/1000)*np.pi/180
    af0 = wiersz_nav [6] #s współczynnik wielomianu do poprawki zegara satelity (opóźnienie)
    af1 = wiersz_nav[7] #s/s współczynnik wielomianu do poprawki zegara (dryft)
    af2 = wiersz_nav[8] #s/s2 współczynnik wielomianu do poprawki zegara (częstotliwość dryftowania)
    gps_week = wiersz_nav[27] #numer tygodnia systemu GPS
    t_razem_z_tygodiami = week * 7 * 86400 + sow
    toe_razem_z_tygodniami = gps_week * 7 * 86400 + Toe
    tk = t_razem_z_tygodiami - toe_razem_z_tygodniami
    dn = wiersz_nav[11]
    Omega_kropka = wiersz_nav[24]
    IDOT = wiersz_nav[25]
    Cuc = wiersz_nav[13]
    Cus = wiersz_nav[15]
    Cic = wiersz_nav[18]
    Cis = wiersz_nav[20]
    Crc = wiersz_nav[22]
    Crs = wiersz_nav[10]
    ni = 3.986005 * 10**14 
    omegaE = 7.2921151467 * 10**-5 
    a = pierwiastek_a**2
    n0 = np.sqrt(ni/(a**3))
    n = n0 + dn
    Mk = M0 + n * tk
    Ek = Mk
    Ei=0 
    while abs(Ek-Ei)>10**-12:
        Ei = Ek
        Ek = Mk + e * math.sin(Ei)
    vk = np.arctan2(np.sqrt(1-e**2)*np.sin(Ek), np.cos(Ek) - e)
    PHI_K = vk + perigee
    duk = Cus* np.sin(2 * PHI_K) + Cuc* np.cos(2 * PHI_K)
    drk = Crs* np.sin(2 * PHI_K) + Crc* np.cos(2 * PHI_K)
    dik = Cis* np.sin(2 * PHI_K) + Cic* np.cos(2 * PHI_K) 
    uk = PHI_K + duk
    rk = a*(1 - e * np.cos(Ek)) + drk
    ik = i0 + IDOT * tk + dik
    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)
    Omega_K = Omega_0 + (Omega_kropka - omegaE) * tk - omegaE * Toe
    Xk = xk * np.cos(Omega_K) - yk*np.cos(ik) * np.sin(Omega_K)
    Yk = xk * np.sin(Omega_K) + yk*np.cos(ik) * np.cos(Omega_K)
    Zk = yk * np.sin(ik)
    c = 299792458.0 #[m/s] – prędkość światła
    dts = af0 + af1 * tk + af2 * (tk**2)
    dt_rel = (-2 * np.sqrt(ni))/(c**2) * e * pierwiastek_a * np.sin(Ek)
    dt_s_rel = dts + dt_rel
    XYZ = [Xk, Yk, Zk]
    return Xk, Yk, Zk, dts, dt_s_rel

def klobuchar(phi, lamb, el, az, tgps):
    c=299792458.0
    alfa = [2.4214E-08, 0.0000E+00, -1.1921E-07, 5.9605E-08]
    beta = [1.2902E+05, 0.0000E+00, -1.9661E+05, -6.5536E+04]
    phis = phi/180
    lambs = lamb/180
    els = el/180
    azs = az/180
    # 1. kąt geocentryczny
    psi = 0.0137 / (els + 0.11) - 0.022
    # 2. szerokość geograficzna IPP
    phi_ipp = phis  + psi * np.cos(np.deg2rad(az))
    if phi_ipp > 0.416:
        phi_ipp = 0.416
    elif phi_ipp < -0.416:
        phi_ipp = -0.416
    # 3. długość geograficzna IPP
    lamb_ipp = lambs + (psi * np.sin(np.deg2rad(az)) / np.cos(phi_ipp*np.pi))
    # 4. szerokość geomagnetyczna IPP
    phim = phi_ipp + 0.064 * np.cos((lamb_ipp - 1.617) * np.pi)
    #5 Wyznaczenie czasu lokalnego
    t=43200*lamb_ipp+tgps
    t=t%86400
    # 6Wyznaczenie amplitudy opóźnienia jonosferycznego
    AION=alfa[0]/180+alfa[1]/180*phim+alfa[2]/180*phim**2+alfa[3]/180*phim**3
    if AION<0:
        AION=0
    # 7. Wyznaczenie okresu opóźnienia jonosferycznego
    PION=beta[0]+beta[1]*phim+beta[2]*phim**2+beta[3]*phim**3
    #8 Wyznaczenie fazy opóźnienia jonosferycznego
    phi_ION=2*np.pi*(t-50400)/PION
    #9 Funkcja mapująca
    mf=1.0+16.0*((0.53-els)**3)
    #10 Opóźnienie jonosferyczne w kierunku satelity dla częstotliwości L1 GPS, w metrach
    if np.abs(phi_ION)<=np.pi/2:
        deltaIL1=c*mf*(5*(10**(-9))+AION*(1-(phi_ION**2)/2+(phi_ION**4)/24)) #dzien
    if np.abs(phi_ION)>np.pi/2:
        deltaIL1=c*mf*5*(10**(-9)) #noc
    return deltaIL1

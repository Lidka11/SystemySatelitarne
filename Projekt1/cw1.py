import tkinter as tk
from tkinter import ttk
from tkcalendar import Calendar, DateEntry
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from matplotlib.pyplot import rc, rcParams, grid 
import math
from pylab import *
from matplotlib.widgets import Slider
def julday(y,m,d,h=0):
    '''
    Simplified Julian Date generator, valid only between
    1 March 1900 to 28 February 2100
    '''
    if m <= 2:
        y = y - 1
        m = m + 12
    # A = np.trunc(y/100)
    # B = 2-A+np.trunc(A/4)
    # C = np.trunc(365.25*y)
    # D = np.trunc(30.6001 * (m+1))
    # jd = B + C + D + d + 1720994.5
    jd = math.floor(365.25*(y+4716))+math.floor(30.6001*(m+1))+d+h/24-1537.5
    return jd

def read_yuma(almanac_file):
    ''' 
    Reading and parsing YUMA asci format
    INPUT:
        Almanac: YUMA format 
    OUTPUT:
        almanac_data -  type list of list [strings value], number of lists is equal to number of satellite
                        one list contain satellites according to the order:         
                        ['SV ID', 'Health', 'Eccentricity', 'Time of Applicability(s)', 'Inclination(rad)', 
                        'Rate of Right Ascen(r/s)', 'SQRT(A)  (m 1/2)', 'Right Ascen at Week(rad)', 
                        'Argument of Perigee(rad)', 'Mean Anom(rad)', 'Af0(s)', 'Af1(s/s)', 'Week no']
        
    '''
    
    if almanac_file:
        alm = open(almanac_file)
        
        alm_lines = alm.readlines()
        all_sat = []
        for idx, value in enumerate(alm_lines):
            # print(idx, value)
            
            if value[0:3]=='ID:':
                one_sat_block = alm_lines[idx:idx+13]
                one_sat = []
                for line in one_sat_block:
                    data = line.split(':')
                    one_sat.append(float(data[1].strip()))
                all_sat.append(one_sat)
        alm.close()
        all_sat = np.array(all_sat)
        
        return (all_sat)


def read_alm(file):
    '''
    Parameters
    ----------
    file : .alm file
    Returns
    -------
    nav_data : 
    nav_data[0] - svprn
    nav_data[1] - health
    nav_data[2] - eccentricity
    nav_data[3] - SQRT_A (square root of major axis a [m**(1/2)])
    nav_data[4] - Omega (Longitude of ascending node at the beginning of week [deg])
    nav_data[5] - omega (Argument of perigee [deg])
    nav_data[6] - M0 (Mean anomally [deg])
    nav_data[7] - ToA (Time of Almanac [second of week])
    nav_data[8] - delta_i (offset to nominal inclination angle; i = 54 + delta_i [deg])
    nav_data[9] - Omega_dot (Rate of Right Ascension [deg/s * 1000])
    nav_data[10]- Satellite clock offset [ns]
    nav_data[11]- Satellite clock drift [ns/s]
    nav_data[12]- GPS week
   
    '''
    m = 0
    with open(file, "r") as f:
        block = []
        nav_data = []
        for s in f:
            # print(s)
            
            if m<13:
                m+=1
                block.append(s)
            else:
                block_array = np.genfromtxt(block,delimiter=10).T
                nav_data.extend(block_array)
                
                m = 0
                block = []
            
    nav_data = np.array(nav_data)        
    return nav_data

def get_prn_number(nav_data):
    prns = []
    for nav in nav_data:
        nsat = nav[0]
        if 0<nsat<=37:
            prn = int(nsat)
            prns.append(prn)
        elif 38<=nsat<=64:
            prn = 100 + int(nsat-37)
            prns.append(prn)
        elif 111<=nsat<=118:
            prn = 400 + int(nsat-110)
            prns.append(prn)
        elif 201<=nsat<=263:
            prn = 200 + int(nsat-200)
            prns.append(prn)    
        elif 264<=nsat<=310:
            prn = 300 + int(nsat-263)
            prns.append(prn)
        elif 311<=nsat:
            prn = 300 + int(nsat-328)
            prns.append(prn)           
        else: 
            prn = 500 + int(nsat)
            prns.append(prn)
    return prns

def get_alm_data(file):
    nav_data = read_alm(file)
    prns = get_prn_number(nav_data)
    nav_data[:,0] = prns
    return nav_data

def create_prn_alm(nav_data):
    prns = []
    for nav in nav_data:
        nsat = nav[0]
        if 0<nsat<=37:
            prn = 'G'+str(int(nsat)).zfill(2)
            prns.append(prn)
        elif 38<=nsat<=64:
            prn = 'R'+str(int(nsat-37)).zfill(2)
            prns.append(prn)
        elif 111<=nsat<=118:
            prn = 'Q'+str(int(nsat-110)).zfill(2)
            prns.append(prn)
        elif 201<=nsat<=263:
            prn = 'E'+str(int(nsat-200)).zfill(2)
            prns.append(prn)    
        elif 264<=nsat<=310:
            prn = 'C'+str(int(nsat-263)).zfill(2)
            prns.append(prn)
        elif 311<=nsat:
            prn = 'C'+str(int(nsat-328)).zfill(2)
            prns.append(prn)           
        else: 
            prn = 'S'+str(int(nsat)).zfill(2)
            prns.append(prn)
    return prns

def get_alm_data_str(alm_file):
    alm_data = read_alm(alm_file)
    nav_data= read_alm(alm_file)
    prn= get_prn_number(alm_data) # prn- numer kodu satelity
    nav_data[:,0]=prn
    prns = create_prn_alm(alm_data) 
    return alm_data, prns    

def get_gps_time(y, m, d, h = 0, mnt = 0, s = 0):
    days = julday(y, m, d) - julday(1980, 1, 6)
    week = days // 7
    day = days%7
    sow = day * 86400 + h * 3600 + mnt * 60 + s
    return int(week), sow

def obliczenie_wspolrzednych_satelity(wiersz_nav, y, m, d, h=0, mnt=0, s=0):
    mi = 3.986005 * 10**14
    OMGE = 7.2921151467 * 10**-5
    svprn = wiersz_nav[0]
    health = wiersz_nav[1]
    e = wiersz_nav[2]
    toa = wiersz_nav[7]
    i = (54 + wiersz_nav[8]) * np.pi/180
    Omega_dot = (wiersz_nav[9]/1000) * np.pi/180
    sqrtA = wiersz_nav[3]
    Omega = wiersz_nav[4] * np.pi/180
    omega = wiersz_nav[5] * np.pi/180
    M0 = wiersz_nav[6] * np.pi/180
    af0 = wiersz_nav[10]
    af1 = wiersz_nav[11]
    gps_week = wiersz_nav[12]
    week,t=get_gps_time(y,m,d,h,mnt,s)
    t_tygodnie = week * 7 * 86400 +t
    toa_tygodnie = gps_week * 7 * 86400 + toa
    tk = t_tygodnie - toa_tygodnie
    a = sqrtA**2
    n = math.sqrt(mi/a**3)
    M_k = M0 + n * tk
    E_k = M_k
    E_i = 0
    while abs(E_k - E_i) > 10**-12:
        E_i = E_k 
        E_k = M_k + e * np.sin(E_i)
    vk = np.arctan2(np.sqrt(1 - e**2) * np.sin(E_k), np.cos(E_k) - e)
    uk = vk + omega
    rk = a * (1 - e * math.cos(E_k))
    Omega_k = Omega + (Omega_dot - OMGE) * tk - OMGE* toa
    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)
    X = xk * np.cos(Omega_k) - yk * np.cos(i) * np.sin(Omega_k)
    Y = xk * np.sin(Omega_k) + yk * np.cos(i) * np.cos(Omega_k)
    Z = yk * np.sin(i)
    return X, Y, Z

def blh2xyz(phi, lam, h):
    phi_rad = np.deg2rad(phi)
    lam_rad = np.deg2rad(lam)
    a = 6378137
    e2 = 0.00669438002290
    N = a / (np.sqrt(1 - e2 * np.sin(phi_rad) * np.sin(phi_rad)))

    X = (N + h) * np.cos(phi_rad) * np.cos(lam_rad)
    Y = (N + h) * np.cos(phi_rad) * np.sin(lam_rad)
    Z = (N * (1 - e2) + h) * np.sin(phi_rad)
    return np.array([X, Y, Z]) #w radianach

def Rneu(phi, lam):
    phi_rad = np.deg2rad(phi)
    lam_rad = np.deg2rad(lam)
    macierz = np.array([[-np.sin(phi_rad) * np.cos(lam_rad), -np.sin(lam_rad), np.cos(phi_rad) * np.cos(lam_rad)],
                        [-np.sin(phi_rad) * np.sin(lam_rad), np.cos(lam_rad), np.cos(phi_rad) * np.sin(lam_rad)],
                        [np.cos(phi_rad), 0, np.sin(phi_rad)]])
    return macierz
def rysuj_wykres_elewacji(elevations, maska, prns):
    elevations_all_satellites = []
    azimuths_all_satellites = []

    for sat in range(len(elevations[0])):
        sat_elevations = []
        sat_azimuths = []
        for min in elevations:
            sat_elevations.append(min[sat][0])
            sat_azimuths.append(min[sat][1])
        elevations_all_satellites.append(sat_elevations)
        azimuths_all_satellites.append(sat_azimuths)

    fig, ax = plt.subplots() 

    for idx, elevations_satellite in enumerate(elevations_all_satellites):
        filtered_elevations = [el if el > maska else None for el in elevations_satellite]
        sat_index= idx
        prn= prns[sat_index]
        ax.plot(range(len(filtered_elevations)), filtered_elevations, label=f" {prn}")

    ax.set_xlabel('Czas w godzinach od 00:00', fontsize=12)
    ax.set_ylabel('Elewacja [stopnie]', fontsize=12)
    ax.set_title('Wykres elewacji dla wszystkich satelitów w zależności od czasu', fontsize=16)
    godziny = [f"{hour}" for hour in range(24)]
    plt.xticks(np.arange(0, 24 * 60, 60), labels=godziny)

    legend = ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=10)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    legend.set_draggable(True)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True)

    return fig

def rysuj_wykres_widocznych_satelitow(elevations, maska):
    widoczne_satelity = [np.sum(elevations[i][:,0] > maska) for i in range(len(elevations))]
    czas_w_godzinach = [i / 60 for i in range(len(widoczne_satelity))]
    fig, ax = plt.subplots()
    ax.plot(czas_w_godzinach, widoczne_satelity, color= "black", alpha=0.8)
    ax.fill_between(czas_w_godzinach, 0, widoczne_satelity, color="green", alpha= 0.2)
    ax.set_xlabel('Czas w godzinach od 00:00', fontsize=12)
    ax.set_ylabel('Liczba widocznych satelitów', fontsize=12)
    ax.set_title('Wykres widoczności satelitów', fontsize=16)

    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=12)
    return fig


def plot_single_satellite_skyplot(elevations, azimuths, indices, elewacje,godzina, minut, maska=None):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(polar=True)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90 + 10, 10))
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    ax.set_rlim(0, 90)

    for idx in indices:
        sat_elevations = elevations[idx]
        sat_azimuths = azimuths[idx]

        if maska:
            valid_indices = [i for i, el in enumerate(sat_elevations) if el > maska]
            sat_elevations = [sat_elevations[i] for i in valid_indices]
            sat_azimuths = [sat_azimuths[i] for i in valid_indices]

        azimuths_rad = np.radians(sat_azimuths)
        radii = 90 - np.array(sat_elevations)
        ax.plot(azimuths_rad, radii)
        
    for el, az, prn in elewacje:
        ax.plot(np.radians(az), 90 -el, marker='o', markersize=8, linestyle='', label=f'{prn}')
    formatted_time = f"{godzina}:{str(minut).zfill(2)}"
    fig.suptitle('Skyplot dla godziny ' + formatted_time, fontsize=16)
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', labelsize=12)
    plt.legend(fontsize=12)

    return fig  # Dodaj zwrócenie figury

def rysuj_wykres_elewacji_1(elevations, maska, prns):
    import numpy as np
    import matplotlib.pyplot as plt
    
    widoczne_satelit = []
    for i in range(len(elevations)):
        visible = [(j, elevations[i][j, 0] > maska) for j in range(len(elevations[i]))]
        widoczne_satelit.append(visible)
    widoczne_satelit = np.array(widoczne_satelit)
    nazwy_satelitow = [f'{prns[idx]}' for idx in range(len(widoczne_satelit[0]))]
    czasy = list(range(24 * 60))
    fig = plt.figure(figsize=(10, 6))

    for satelita in range(len(widoczne_satelit[0])):
        widocznosci = [nazwy_satelitow[satelita] if widoczne_satelit[minuta][satelita][1] else '' for minuta in range(len(widoczne_satelit))]
        plt.scatter(czasy, widocznosci, marker='o')
    plt.xlabel('Czas [godziny]', fontsize=12)
    plt.ylabel('Nazwa satelity', fontsize=12)
    plt.title('Widoczność satelitów w zależności od czasu', fontsize=16)
    plt.yticks(nazwy_satelitow)
    plt.grid(True, axis='y', linestyle='--')
    plt.tick_params(axis='y', labelsize=12)
    plt.tick_params(axis='x', labelsize=12)
    plt.xticks(np.arange(0, 24 * 60, 60), labels=range(24))
    plt.yticks(fontsize=12)
    return fig



def ground_track(indeksy, godzina, minut, elevations, prns, xyz_pos, maska):
    '''
    Ploting satellite ground track:
        Satellite groundtrack longitude and latitude are calculated 
        from XYZ(ECEF) satellites bythe use of function(above): groundtrack_latlon(xs, ys, zs)
    
    INPUT:
        lon_lat_ground ; list of list [lon, lat]
        Coastline.txt - file with coastline coordinate
    '''
    # GENERAL setup
    # rc('text', usetex=True)
    # rc('font', family='serif');
    # rc('font',family='helvetica');
    rc('grid', color='gray', linewidth=0.1, linestyle='--')
    fontsize = 20
    rc('xtick', labelsize = fontsize)
    rc('ytick', labelsize = fontsize)
    rc('font', size = fontsize)
    params = {'legend.fontsize': 8,  'legend.handlelength': 2}
    rcParams.update(params)
    coastline_data= np.loadtxt('Coastline.txt',skiprows=1)
    w, h = plt.figaspect(0.5)
    fig = plt.figure(figsize=(w,h))
    ax = fig.gca()
    formatted_time = f"{godzina}:{str(minut).zfill(2)}"
    fig.suptitle('GPS Ground Track dla godziny ' + formatted_time, fontsize=16)
    plt.plot(coastline_data[:,0],coastline_data[:,1],'g', linewidth=0.5);
    ax.set_xlabel('Longitude $[\mathrm{^\circ}]$',fontsize=14)
    ax.set_ylabel('Latitude $[\mathrm{^\circ}]$',fontsize=14)
    plt.xlim(-180,180);
    plt.ylim(-90,90);
    plt.yticks([-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90]);
    plt.xticks([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180]);
    for i in range(len(elevations)):
        visible_satellites = np.where(elevations[i][:, 0] > maska)[0]
        for j in indeksy:
            if j in visible_satellites:
                lon = np.rad2deg(np.arctan2(xyz_pos[i][j][1], xyz_pos[i][j][0]))
                lat = np.rad2deg(np.arcsin(xyz_pos[i][j][2] / np.linalg.norm(xyz_pos[i][j])))
                plt.plot(lon, lat, color=(1, 0, 0, 0.5), marker='o', linestyle='', markersize=0.5)            
            
    
    indeks_godziny = godzina * 60 + minut
    
    for j in indeksy:
        if elevations[indeks_godziny][j][0] > maska:
            lon = np.rad2deg(np.arctan2(xyz_pos[indeks_godziny][j][1], xyz_pos[indeks_godziny][j][0]))
            lat = np.rad2deg(np.arcsin(xyz_pos[indeks_godziny][j][2] / np.linalg.norm(xyz_pos[indeks_godziny][j])))
            plt.plot(lon, lat, 'bo', markersize=6)
            plt.text(lon, lat, str(prns[j]), fontsize=10, color='black')
    ax.grid(True)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    # Zwracanie obiektu Figure
    return fig

global godzina, minut, maska, phi, lam, h, y, m, d


def najw_oblicz(phi, lam, h, y, m, d, godzina, minut, maska):
    nav= get_alm_data('Almanac2024053.alm')

    #nav2, prns= get_alm_data_str('Almanac2024053.alm')
    prns= create_prn_alm(nav)
    satelity= nav[:,0]<100
    nav= nav[satelity,:]

    wiersz_nav = nav[0,:]
    elevations=[]
    xyz_pos=[]


    A=[]
    for t in range(0,24):
        for min in range(0,60):
            xyz_sat=[]
            el_ar=[]
            A_rows = []
            for wiersz_nav in nav:
                xyz_satelity=np.array(obliczenie_wspolrzednych_satelity(wiersz_nav, y,m,d,t, min,0))
                if len(xyz_satelity) == 3:
                    xyz_sat.append(xyz_satelity)
                a = 6378137
                e2 = 0.00669438002290
                xyz_odbiornika = np.array(blh2xyz(phi, lam, h)) 
                N= a/(np.sqrt(1-e2*(np.sin(np.deg2rad(phi))**2)))
                wektor_s_o = xyz_satelity - xyz_odbiornika
                phi_rad = np.deg2rad(phi)
                lam_rad = np.deg2rad(lam)
                R = np.array([[-np.sin(phi_rad) * np.cos(lam_rad), -np.sin(lam_rad), np.cos(phi_rad)*np.cos(lam_rad)],
                                        [-np.sin(phi_rad) * np.sin(lam_rad), np.cos(lam_rad), np.cos(phi_rad) * np.sin(lam_rad)],
                                        [np.cos(phi_rad), 0, np.sin(phi_rad)]])
                xrneu= R.T@wektor_s_o
                neu= xrneu
                az= np.rad2deg(np.arctan2(neu[1], neu[0]))
                if az<0:
                    az=az+360
                #print("az",az)
                el= np.rad2deg(np.arcsin(neu[2]/np.sqrt(neu[0]**2+neu[1]**2+neu[2]**2)))
                #print("el",el)
                el_ar.append(np.array([el,az]))
                psr= np.sqrt((xyz_satelity[0]-xyz_odbiornika[0])**2+(xyz_satelity[1]-xyz_odbiornika[1])**2+(xyz_satelity[2]-xyz_odbiornika[2])**2)
                if el>maska:
                    wiersz_macierzy_A= np.array([(xyz_odbiornika[0]-xyz_satelity[0])/psr,(xyz_odbiornika[1]-xyz_satelity[1])/psr,(xyz_odbiornika[2]-xyz_satelity[2])/psr, 1]) 
                    A_rows.append(wiersz_macierzy_A)
            A.append(np.array(A_rows))
            elevations.append(np.array(el_ar))
            xyz_pos.append(np.array(xyz_sat))
    elevations= np.array(elevations)
    xyz_pos= np.array(xyz_pos)
    def oblicz_dop(Q):
        Q_zmniejszona= Q[:3,:3]
        Qneu= R.T@Q_zmniejszona@R
        GDOP= np.sqrt(Q[0,0]+Q[1,1]+Q[2,2]+Q[3,3])
        PDOP= np.sqrt(Q[0,0]+Q[1,1]+Q[2,2])
        TDOP= np.sqrt(Q[3,3])
        HDOP= np.sqrt(Qneu[0,0]+Qneu[1,1])
        VDOP= np.sqrt(Qneu[2,2])
        PDOPneu= np.sqrt(Qneu[0,0]+Qneu[1,1]+Qneu[2,2])
        return TDOP, PDOP, GDOP, HDOP, VDOP, PDOPneu

    DOPS = []
    for a in A:
        try:
            q = np.linalg.inv(a.T.dot(a))
        except Exception as e:
            continue
        dops_to_append = oblicz_dop(q)
        dops_to_append = np.array(dops_to_append)
        DOPS.append(dops_to_append)

    A = np.array(
        [a for a in A for a in a]
    )
    Q = np.linalg.inv(A.T@A)
    DOPS = np.array(DOPS)

    elevations_all_satellites = []
    azimuths_all_satellites = []

    for sat in range(len(elevations[0])):
        sat_elevations = []
        sat_azimuths = []
        for min in elevations:
            sat_elevations.append(min[sat][0])
            sat_azimuths.append(min[sat][1])
        elevations_all_satellites.append(sat_elevations)
        azimuths_all_satellites.append(sat_azimuths)
    indeksy=[]
    elewacje = []
    for idx, (el, az) in enumerate(elevations[godzina*60+minut]):
        if el > maska:
            sat_index = idx 
            indeksy.append(sat_index)
            prn= prns[sat_index]
            elewacje.append((el, az, prn))
    nazwy_dops = ["GDOP", "PDOP", "TDOP", "HDOP", "VDOP", "PDOPneu"]
    kolory = ["r", "g", "b", "y", "c", "m"]
    return  elevations, elevations_all_satellites, azimuths_all_satellites, indeksy, elewacje, DOPS, xyz_pos, nazwy_dops, kolory, Q, prns

def on_closing():
    root.quit()


def draw_plot_in_frame2(elevations, maska, prns):
    fig = rysuj_wykres_elewacji(elevations, maska, prns)
    canvas = FigureCanvasTkAgg(fig, master=frame2)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def draw_plot_in_frame3(DOPS, xyz_pos, nazwy_dops, kolory):
    fig = create_dops_plot(DOPS, xyz_pos, nazwy_dops, kolory)
    canvas = FigureCanvasTkAgg(fig, master=frame3)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def draw_plot_in_frame4(elevations, maska):
    fig = rysuj_wykres_widocznych_satelitow(elevations, maska)
    canvas = FigureCanvasTkAgg(fig, master=frame4)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def draw_plot_in_frame5(elevations, maska, prns):
    fig = rysuj_wykres_elewacji_1(elevations, maska, prns)
    canvas = FigureCanvasTkAgg(fig, master=frame5)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def draw_plot_in_frame6(elevations_all_satellites, azimuths_all_satellites, indeksy,elewacje,godzina, minut,  maska):
    fig = plot_single_satellite_skyplot(elevations_all_satellites, azimuths_all_satellites, indeksy,elewacje,godzina, minut, maska=0)
    canvas = FigureCanvasTkAgg(fig, master=frame6)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def draw_plot_in_frame7(indeksy, godzina, minut, elevations, prns, xyz_pos, maska):
    fig = ground_track(indeksy, godzina, minut, elevations, prns, xyz_pos, maska)
    canvas = FigureCanvasTkAgg(fig, master=frame7)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def create_dops_plot(DOPS, xyz_pos, nazwy_dops, kolory):
    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(DOPS.shape[1]):
        ax.plot(range(len(xyz_pos)), DOPS[:, i], label=nazwy_dops[i], color=kolory[i])
        
    ax.set_xlabel('Czas w godzinach od 00:00', fontsize=12)
    ax.set_ylabel('Wartość DOPS', fontsize=12)
    ax.set_title('Wykres DOPS w zależności od czasu', fontsize=16)
    godziny = [f"{hour}" for hour in range(24)]  # Formatowanie godzin
    plt.xticks(np.arange(0, 24 * 60, 60), labels=godziny)
    ax.legend(fontsize=12)
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=12)
    return fig

def get_date():
    selected_date = cal.selection_get()
    e_date.delete(0, tk.END)
    e_date.insert(0, selected_date)
def show_calendar():
    cal.place(x=300, y=200)
    b_date.config(command=hide_calendar)

def hide_calendar(event=None):
    selected_date_str = cal.get_date()
    selected_date = datetime.datetime.strptime(selected_date_str, '%m/%d/%y')
    formatted_date = selected_date.strftime('%Y-%m-%d')
    e_date.delete(0, tk.END)
    e_date.insert(0, formatted_date)
    cal.place_forget()
    b_date.config(command=show_calendar)


def clear_frame(frame):
    for widget in frame.winfo_children():
        widget.destroy()

def reset_font_config():
    rc('grid', color='gray', linewidth=0.1, linestyle='--')
    fontsize = 10
    rc('xtick', labelsize=fontsize)
    rc('ytick', labelsize=fontsize)
    rc('font', size=fontsize)
    params = {'legend.fontsize': 5, 'legend.handlelength': 1}
    rcParams.update(params)

def pobierz_dane():
    global frame2, frame3, frame4, frame5, frame6, frame7
    clear_frame(frame2)
    clear_frame(frame3)
    clear_frame(frame4)
    clear_frame(frame5)
    clear_frame(frame6)
    clear_frame(frame7)
    #reset_font_config()
    lam = int(e1.get())
    phi = int(e2.get())
    godzinaa= e3.get()
    godzina= int(godzinaa)
    h=int(e12.get())
    minutt=e4.get()
    minut = int(minutt)
    maska = int(e7.get())
    if not (0 <= phi <= 90):
        print("Błąd: Phi musi być z przedziału [0, 90].")
    if not (0 <= godzina <= 23):
        print("Błąd: Godzina musi być z przedziału [0, 23].")
    if not (0 <= minut <= 59):
        print("Błąd: Minuta musi być z przedziału [0, 59].")
    if not (0 <= maska <= 180):
        print("Błąd: Maska musi być z przedziału [0, 180].")
    selected_date = e_date.get()
    if selected_date:
        try:
            global y, m, d
            y, m, d = map(int, selected_date.split('-'))
        except ValueError:
            print("Błąd: Nieprawidłowy format daty. Użyj formatu RRRR-MM-DD.")
            return
    print("Wartości zapisane:", lam, phi, h, minut, maska, y, m, d)
    elevations, elevations_all_satellites, azimuths_all_satellites, indeksy, elewacje, DOPS, xyz_pos, nazwy_dops, kolory, Q, prns = najw_oblicz(phi, lam, h, y, m, d, godzina, minut, maska)
    draw_plot_in_frame2(elevations, maska, prns)
    draw_plot_in_frame3(DOPS, xyz_pos, nazwy_dops, kolory)
    draw_plot_in_frame4(elevations, maska)
    draw_plot_in_frame5(elevations, maska, prns)
    draw_plot_in_frame6(elevations_all_satellites, azimuths_all_satellites, indeksy, elewacje,godzina, minut, maska)
    draw_plot_in_frame7(indeksy, godzina, minut, elevations, prns, xyz_pos, maska)

# Tworzenie głównego okna
root = tk.Tk()
root.title("GNSS - Interfejs Użytkownika")
root.geometry("1200x700")

notebook = ttk.Notebook(root)
notebook.pack(fill='both', expand=True)

frame1 = ttk.Frame(notebook)
frame2 = ttk.Frame(notebook)
frame3 = ttk.Frame(notebook)
frame4 = ttk.Frame(notebook)
frame5 = ttk.Frame(notebook)
frame6= ttk.Frame(notebook)
frame7= ttk.Frame(notebook)

notebook.add(frame1, text='Dane')
notebook.add(frame2, text='Elewacje od czasu')
notebook.add(frame3, text='DOPS')
notebook.add(frame4, text='Widoczność satelitów')
notebook.add(frame5, text='Widoczność satelitów od czasu')
notebook.add(frame6, text='Skyplot')
notebook.add(frame7, text='Groundtrack')


t0= tk.Label(frame1, text='Podaj swoje współrzędne', font= ('Arial', 10, 'bold'))
t0.place(x=50,y=50)
t1= tk.Label(frame1, text='lambda')
t1.place(x=50,y=80)
e1= tk.Spinbox(frame1, from_=0, to=180, width=5)
e1.place(x=100,y=80)
e1.delete(0, tk.END) 
e1.insert(0, '21')
t2= tk.Label(frame1, text='phi')
t2.place(x=50,y=105)
e2= tk.Spinbox(frame1, from_=0, to=90, width=5)
e2.place(x=100,y=105)
e2.delete(0, tk.END)  
e2.insert(0, '52')
t12= tk.Label(frame1, text='h')
t12.place(x=50,y=130)
e12= tk.Spinbox(frame1, from_=0, to=500, width=5)
e12.place(x=100,y=130)
e12.delete(0, tk.END)  
e12.insert(0, '100')
t3= tk.Label(frame1, text='godz')
t3.place(x=100,y=170)
e3= tk.Spinbox(frame1, from_=0, to=23, width=5)
e3.place(x=80,y=200)
e3.delete(0, tk.END)  
e3.insert(0, '12')
t4= tk.Label(frame1, text='min')
t4.place(x=150,y=170)
e4= tk.Spinbox(frame1, from_=0, to=59, width=5)
e4.place(x=150,y=200)
t5= tk.Label(frame1, text=':')
t5.place(x=138,y=198)
t6= tk.Label(frame1, text='Podaj godzinę dla skyplot i groundtrack', font= ('Arial', 10, 'bold'))
t6.place(x=50,y=150)
t7= tk.Label(frame1, text='Podaj maskę', font= ('Arial', 10, 'bold'))
t7.place(x=50,y=290)
available_values = ['5', '10', '15']
e7 = ttk.Combobox(frame1, values=available_values, width=5)
e7.set('10') 
e7.place(x=80,y=320)
e7.delete(0, tk.END)  
e7.insert(0, '10')
t8= tk.Label(frame1, text='Podaj datę', font= ('Arial', 10, 'bold'))
t8.place(x=50,y=230)
e_date= tk.Entry(frame1)
e_date.place(x=50,y=260)
cal= Calendar(frame1, selectmode='day', year=2024, month=2, day=29)
cal.bind("<<CalendarSelected>>", hide_calendar)
b_date= tk.Button(frame1, text='Wybierz datę', command=show_calendar)
b_date.place(x=200,y=255)
but1= tk.Button(frame1, text='Rysuj wykresy', command=pobierz_dane)
but1.place(x=100,y=360)
root.protocol("WM_DELETE_WINDOW", on_closing)

root.mainloop()

from funkcje import*
import matplotlib.pyplot as plt

# scieżka do pliku nawigacyjnego
nav_file = "C:\sem 4\Systemy nawigacji satelitarnej\projekt2\BRDC00WRD_R_20240650000_01D_GN.rnx"
# scieżka do pliku obserwacyjnego
obs_file = "C:\sem 4\Systemy nawigacji satelitarnej\projekt2\JOZ200POL_R_20240650000_01D_30S_MO.rnx"

time_start= [2024, 3, 5, 0, 0, 0]  
time_end =    [2024, 3, 5, 23, 59, 59] 
obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
nav, inav = readrnxnav(nav_file)
zdrowe = nav[:, 30] == 0
nav = nav[zdrowe, :]
inav = inav[zdrowe]

el_mask = 10 # maska elewacji
week, tow = date2tow(time_start)[0:2]
week_end, tow_end=date2tow(time_end)[0:2]
t=tow


sats = np.array([25, 31, 32, 29, 28, 24, 20, 11, 12,  6])
dt=30
index_t=iobs[:,2]==t
Pobs= obs[index_t,0]
sats=iobs[index_t,0]

#wspolrzedne odbiornika
phi_odbr = 52.100390639667985 
lam_odbr = 20.932535365103107 
h_odbr = -12258.136260229163
tgps = 43200
t0=291.15
p0=1013.25
Rh0=0.5
c1=77.64
c2=-12.96
c3=3.718*10**5
h=152

xr0 = [3660000.,  1400000.,  5000000.]
dtr=0
tau=0.07
xyz_r_ref = [3664880.9100,  1409190.3850,  5009618.2850] # wspolrzedne referencyjne 
we=7.2921151467*10**-5
c = 299792458.0 
visible_sats = []
dx_val=[]
dy_val=[]
dz_val=[]

consider_tropo = input("Uwzględnić poprawkę troposferyczną? (tak/nie): ").strip().lower() == 'tak'
consider_iono = input("Uwzględnić poprawność jonosferyczna? (tak/nie): ").strip().lower() == 'tak'

for t in range(tow, tow_end+1, 30): 
    index_t = iobs[:, 2]==t
    pobs = obs[index_t, 0]
    sats = iobs[index_t, 0]
    tau_list = [0.07] * len(sats)
    for i in range(5):
        A = []
        Y = []
        X = []
        tau_new_list = []
        visible_sat = []
        dx_values=[]
        dy_values=[]
        dz_values=[]
        for id, sat in enumerate(sats):
            ts = t - tau_list[id] + dtr
            xs0, ys0, zs0, dts, dt_s_real = obl_wspolrzednych_satelity(ts, week, sat, inav, nav)
            xyz0_s = np.array([xs0, ys0, zs0]) 
            alfa=we*tau_list[id]
            Rz=np.array([[np.cos(alfa), np.sin(alfa), 0], [-np.sin(alfa), np.cos(alfa), 0], [0, 0, 1]])
            xs,ys,zs = Rz@[xs0, ys0, zs0]
            xyz_s = np.array([xs, ys, zs]) 
            rho = np.sqrt((xs - xr0[0])**2 + (ys - xr0[1])**2 + (zs - xr0[2])**2)
            tau_new = rho/c
            tau_new_list.append(tau_new)
            wektor_s_o = xyz_s - xr0 
            phi_r, lam_r, h_r = hirvonen(xr0[0], xr0[1], xr0[2])
            R_h = Rneu(phi_r, lam_r)
            neu = R_h.T@wektor_s_o
            az = np.rad2deg(np.arctan2(neu[1], neu[0]))
            if az < 0:
                az = az + 360
            el = np.rad2deg(np.arcsin(neu[2] / np.sqrt(neu[0]**2 + neu[1]**2 + neu[2]**2)))
            if el > el_mask:
                if i == 0:
                    dT = 0
                    dJ = 0
                    tropo=0
                    deltaIL1=0
                else:
                    hort=h-31
                    p=p0*(1-0.0000226*hort)**5.225
                    temp=t0-0.0065*hort
                    Rh=Rh0*mat.exp(-0.0006396*hort)
                    Nd0=c1*(p/temp)
                    e=6.11*Rh*10**((7.5*(temp-273.15))/(temp-35.85))
                    Nw0 = c2 * (e/temp) + c3 * (e/(temp**2))
                    N0=Nd0+Nw0
                    hd= 40136+148.72*(temp-273.15)
                    dTdo=0.002277*p
                    dTwo=0.002277*((1255/temp)+0.05)*e
                    tropo=(1/np.sin(np.deg2rad(el)))*(dTdo+dTwo)
                    deltaIL1=klobuchar(xr0[0], xr0[1], el, az, t)
                    visible_sat.append(sat)
                pcalc = rho + c*dtr - c*dt_s_real
                if consider_tropo:
                    pcalc+=tropo
                if consider_iono:
                    pcalc+=deltaIL1
                y = pobs[id] - pcalc
                Y.append(y)
                A_wiersz = [-(xyz_s[0] - xr0[0])/rho,  -(xyz_s[1] - xr0[1] )/rho, -(xyz_s[2] - xr0[2])/rho , 1]
                A.append(A_wiersz)

        A_ar = np.array(A)
        Y_ar = np.array(Y)
        X =(np.linalg.inv(A_ar.T@A_ar))@A_ar.T@Y_ar
        xr0 += X[:3]
        dx=xr0[0] - xyz_r_ref[0]
        dy=xr0[1] - xyz_r_ref[1]
        dz=xr0[2] - xyz_r_ref[2]
        dtr = dtr + X[3]/c
        for ix in range(len(tau_new_list)):
            tau_list[ix] = tau_new_list[ix]
    visible_sats.append(visible_sat)
    dx_val.append(dx)
    dy_val.append(dy)
    dz_val.append(dz)

time = np.arange(0, len(dx_val) * 30, 30) / 3600  # 30 sekundowe odstępy, przeliczane na godziny
sigma_dx = np.std(dx_val)
rms_dx = np.sqrt(np.mean(np.square(dx_val)))

sigma_dy = np.std(dy_val)
rms_dy = np.sqrt(np.mean(np.square(dy_val)))

sigma_dz = np.std(dz_val)
rms_dz = np.sqrt(np.mean(np.square(dz_val)))

plt.figure(figsize=(12, 8))
plt.subplot(311)
plt.plot(time, dx_val, label='dx', color='blue')
plt.xlabel('Czas [godziny]')
plt.ylabel('dx [m]')
plt.grid()
plt.text(0.05, 0.95, f'$\\sigma = {sigma_dx:.2f}$\nRMS = {rms_dx:.2f}m', 
         verticalalignment='top', horizontalalignment='left', 
         transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

plt.subplot(312)
plt.plot(time, dy_val, label='dy', color='dodgerblue')
plt.xlabel('Czas [godziny]')
plt.ylabel('dy [m]')
plt.grid()
plt.text(0.05, 0.95, f'$\\sigma = {sigma_dy:.2f}$\nRMS = {rms_dy:.2f}m', 
         verticalalignment='top', horizontalalignment='left', 
         transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

plt.subplot(313)
plt.plot(time, dz_val, label='dz', color='royalblue')
plt.xlabel('Czas [godziny]')
plt.ylabel('dz [m]')
plt.grid()
plt.text(0.05, 0.95, f'$\\sigma = {sigma_dz:.2f}$\nRMS = {rms_dz:.2f}m', 
         verticalalignment='top', horizontalalignment='left', 
         transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

plt.tight_layout()
plt.show()

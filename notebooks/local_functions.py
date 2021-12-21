import matplotlib.pyplot as plt
import numpy as np
import datetime


# During feeding we have N_ext concentration and it changes due to growth for some time
def N_feeding_original(x,t,Nintmax,Nintmin,Vmax,Ks,KN,dNextoutdt,dNextindt,miu,dmoutdt,Nintcrit):
    """ See powerpoint with the formulae 
    During feeding we supply an initially high Next and let the algae to consume it 
    """
    Next = x[0] # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin) # units: [-]
    uN = Vmax * Next / (Ks + Next) # units: [umol N/g DW/h]
    #fN = (Nint - Nintmin)/(KN + (Nint - Nintmin)) # units: [-]
    #fN = ((Nint - Nintmin)/Nint) / ((Nintmax - Nintmin)/Nintmax) # units: [-]
    if Nint >= Nintcrit:
        fN = 1
    else:
        fN = ((Nint - Nintmin)/Nint) / ((Nintcrit - Nintmin)/Nintcrit)  # units: [ ]
    # fD =  Kd**2/(Kd**2+m**2) #test - density effect
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    dNextdt = -Neff * uN * m - dNextoutdt + dNextindt # units: [umol N/l/h]
    dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * fN  #units: [%g N/g DW/h]
    dmdt = miu * fN * m - dmoutdt #units: [g DW/l/h]
    
    return [dNextdt,dNintdt,dmdt]

# added density function
def N_feeding_density(x,t,Nintmax,Nintmin,Vmax,Ks,KN,dNextoutdt,dNextindt,miu,dmoutdt,Kd):
    """ See powerpoint with the formulae 
    During feeding we supply an initially high Next and let the algae to consume it 
    """
    Next = x[0] # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin) # units: [-]
    uN = Vmax * Next / (Ks + Next) # units: [umol N/g DW/h]
    fN = (Nint - Nintmin)/(KN + (Nint - Nintmin)) # units: [-]
    # fN = ((Nint - Nintmin)/Nint) / ((Nintmax - Nintmin)/Nintmax) # units: [-]
    fD =  Kd**2/(Kd**2+m**2) #test - density effect
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    dNextdt = -Neff * uN * m - dNextoutdt + dNextindt # units: [umol N/l/h]
    dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * fN  #units: [%g N/g DW/h]
    dmdt = miu * fN * fD * m - dmoutdt #units: [g DW/l/h]
    
    return [dNextdt,dNintdt,dmdt]


def N_feeding_constant_Nint(x,t,Nintmax,Nintmin,Vmax,Ks,KN,dNextoutdt,dNextindt,miu,dmoutdt):
    """ See powerpoint with the formulae 
    During feeding we supply an initially high Next and let the algae to consume it 
    """
    Next = x[0] # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin) # units: [-]
    uN = Vmax * Next / (Ks + Next) # units: [umol N/g DW/h]
    fN = (Nint - Nintmin)/(KN + (Nint - Nintmin)) # units: [-]
    # fN = ((Nint - Nintmin)/Nint) / ((Nintmax - Nintmin)/Nintmax) # units: [-]
    # fD =  Kd**2/(Kd**2+m**2) #test - density effect
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    dNextdt = -Neff * uN * m - dNextoutdt + dNextindt # units: [umol N/l/h]
    dNintdt = 0 #umol_to_percent_DW * Neff * uN - Nint * miu * fN  #units: [%g N/g DW/h]
    dmdt = miu * fN * m - dmoutdt #units: [g DW/l/h]
    
    return [dNextdt,dNintdt,dmdt]

def N_growing(x,t,Nintmax,Nintmin,Vmax,Ks,KN,dNextoutdt,dNextindt,miu,dmoutdt):
    """ See powerpoint with the formulae 
    ! we assume during growing that there is Next = 0 
    """

    Next = 0 # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin) # units: [-]
    uN = Vmax * Next / (KN + Next) # units: [umol N/g DW/h]
    fN = (Nint - Nintmin)/(Ks + (Nint - Nintmin)) # units: [-]
    # fN = ((Nint - Nintmin)/Nint) / ((Nintmax - Nintmin)/Nintmax) # units: [-]
    # fD =  Kd**2/(Kd**2+m**2) #test - density effect
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    dNextdt = 0 # units: [umol N/l/h]
    dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * fN  #units: [%g N/g DW/h]
    dmdt = miu * fN * m - dmoutdt #units: [g DW/l/h]
    
    return [dNextdt,dNintdt,dmdt]

def N_growing_constant_Nint(x,t,Nintmax,Nintmin,Vmax,Ks,KN,dNextoutdt,dNextindt,miu,dmoutdt):
    """ See powerpoint with the formulae 
    ! we assume during growing that there is Next = 0 
    """

    Next = 0 # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin) # units: [-]
    uN = Vmax * Next / (KN + Next) # units: [umol N/g DW/h]
    fN = (Nint - Nintmin)/(Ks + (Nint - Nintmin)) # units: [-]
    # fN = ((Nint - Nintmin)/Nint) / ((Nintmax - Nintmin)/Nintmax) # units: [-]
    # fD =  Kd**2/(Kd**2+m**2) #test - density effect
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    dNextdt = 0 # units: [umol N/l/h]
    dNintdt = 0 #umol_to_percent_DW * Neff * uN - Nint * miu * fN  #units: [%g N/g DW/h]
    dmdt = miu * fN * m - dmoutdt #units: [g DW/l/h]
    
    return [dNextdt,dNintdt,dmdt]

def plot_result_extra(t_model,Next_model,Nint_model,m_model,Nint=None, yerrNint=None, tNint=None,m=None, yerrm=None, tm=None, Next=None, yerrNext=None, tNext=None,Next_val = None,tNext_val=None,Nint_val = None,tNint_val=None,m_val = None,tm_val=None,Next_spor = None,tNext_spor=None,Nint_spor = None,tNint_spor=None,
                      m_spor = None,tm_spor=None):
    """ Plot time series of the results """
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(12,10))
    ax[0].plot(t_model,Next_model,'.--',color = 'dodgerblue')
    ax[0].errorbar(t_model[::],Next_model[::],yerr=yerrNext,fmt='.',color='black',ecolor='lightblue',
                       elinewidth=1, capsize=0)
    if Next is not None:
        ax[0].errorbar(tNext[::],Next[::],yerr=None,fmt='.',color='black',ecolor='lightgray',
                       elinewidth=3, capsize=0)
        ax[0].errorbar(tNext_val[::],Next_val[::],yerr=None,fmt='*',color='orange',ecolor='lightgray',
                       elinewidth=3, capsize=0)
        ax[0].errorbar(tNext_spor[::],Next_spor[::],yerr=None,fmt='.',color='red',ecolor='lightgray',
               elinewidth=3, capsize=0)
        ax[0].plot()
        # ax[0].plot(tNext[::],Next[::],'ko')
    ax[0].set_ylabel('Next \n [µmol $L^{-1}$]', fontsize=16, weight='bold')
    #ax[0].set_ylim([0,Next_model[0]])
    
    # plt.axis([0, 180, 0, 430])
    ax[1].plot(t_model,Nint_model,'x--', color = 'lightgreen')
    ax[1].set_ylabel('Nint \n [% g N $g DW^{-1}$]\n', fontsize=16, weight='bold')
    ax[1].set_ylim([0.5,4.5])
    ax[1].errorbar(t_model[::],Nint_model[::],yerr=yerrNint,fmt='.',color='black',ecolor='lightgreen',
                   elinewidth=1, capsize=0)
    if Nint is not None:
        ax[1].errorbar(tNint[::],Nint[::],yerr=None,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        ax[1].errorbar(tNint_val[::],Nint_val[::],yerr=None,fmt='*',color='orange',ecolor='lightgray',
        elinewidth=3, capsize=0)
        ax[1].errorbar(tNint_spor[::],Nint_spor[::],yerr=None,fmt='.',color='red',ecolor='lightgray',
        elinewidth=3, capsize=0)
        # Add Continuous Errors? https://jakevdp.github.io/PythonDataScienceHandbook/04.03-errorbars.html
        ax[1].plot()
        # ax[1].plot(tNint[::],Nint[::],'ko')
        
    # plt.axis([0, 180, 2.0, 4.0])
    ax[2].plot(t_model,m_model,'+--',color = 'green')
    ax[2].set_xlabel('time [hours]')
    ax[2].errorbar(t_model[::],m_model[::],yerr=yerrm,fmt='.',color='black',ecolor='lightgreen',elinewidth=1, capsize=0)
    if m is not None:        
        ax[2].errorbar(tm[::],m[::],yerr=None,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        ax[2].errorbar(tm_val[::],m_val[::],yerr=None,fmt='*',color='orange',ecolor='lightgray',
        elinewidth=3, capsize=0)
        ax[2].errorbar(tm_spor[::],m_spor[::],yerr=None,fmt='.',color='red',ecolor='lightgray',
        elinewidth=3, capsize=0)
        ax[2].plot()
        # ax[2].plot(tm[::],m[::],'ko')
    ax[2].set_ylabel('m \n [g DW $L^{-1}$]', fontsize=16, weight='bold')
    ax[2].set_ylim([0,0.6])
    # plt.axis([0, 180, 0.03, 0.3])
    # plot_result(np.hstack(T),np.hstack(NEXT),np.hstack(NINT),np.hstack(M))
    # ax[0].plot(df.iloc[:][1::2],df.iloc[:][2::2],'--s')

    
def plot_result_starvation(t_model,Next_model,Nint_model,m_model):
    """ Plot time series of the results """
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(2,1,sharex=True,figsize=(12,6))
    ax[0].plot(t_model,Nint_model,'x--', color = 'lightgreen')
    ax[0].set_ylabel('Nint \n [% g N $g DW^{-1}$]', fontsize=12, weight='bold')
    ax[0].set_ylim([0.5,4.5])
    ax[1].plot(t_model,m_model,'+--',color = 'green')
    ax[1].set_xlabel('time [hours]')
    ax[1].set_ylabel('m \n [g DW $L^{-1}$]', fontsize=12, weight='bold')
    ax[1].set_ylim([0,1])
    
def plot_result(t_model,Next_model,Nint_model,m_model,Nint=None, yerrNint=None, tNint=None,m=None, yerrm=None, tm=None, Next=None, yerrNext=None, tNext=None):
    """ Plot time series of the results """
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(12,10))
    ax[0].plot(t_model,Next_model,'r.--')
    if Next is not None:
        ax[0].errorbar(tNext[::],Next[::],yerr=yerrNext,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        ax[0].plot()
        # ax[0].plot(tNext[::],Next[::],'ko')
    ax[0].set_ylabel('Next \n [umol / l]', fontsize=12, weight='bold')
    # plt.axis([0, 180, 0, 430])
    ax[1].plot(t_model,Nint_model,'gx--')
    ax[1].set_ylabel('Nint \n [% g N / g DW]', fontsize=12, weight='bold')
    if Nint is not None:
        ax[1].errorbar(tNint[::],Nint[::],yerr=yerrNint,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        # Add Continuous Errors? https://jakevdp.github.io/PythonDataScienceHandbook/04.03-errorbars.html
        ax[1].plot()
        # ax[1].plot(tNint[::],Nint[::],'ko')
        
    # plt.axis([0, 180, 2.0, 4.0])
    ax[2].plot(t_model,m_model,'b+--')
    ax[2].set_xlabel('time [hours]')
    if m is not None:
        ax[2].errorbar(tm[::],m[::],yerr=yerrm,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        ax[2].plot()
        # ax[2].plot(tm[::],m[::],'ko')
    ax[2].set_ylabel('m \n [g DW / l]', fontsize=12, weight='bold')
    # plt.axis([0, 180, 0.03, 0.3])
    # plot_result(np.hstack(T),np.hstack(NEXT),np.hstack(NINT),np.hstack(M))
    # ax[0].plot(df.iloc[:][1::2],df.iloc[:][2::2],'--s')
    
    
def remove_sporulation_from_m_list(m_all, thrashold):
    return [value for value in m_all if value <= thrashold]

def Time_to_Hours(T1,T2,h=None):
    T1_arr = T1.split(',')
    T1_arr = [int(i) for i in T1_arr]
    T2_arr = T2.split(',')
    T2_arr = [int(i) for i in T2_arr] 
    new_T1 = datetime.datetime(T1_arr[0],T1_arr[1],T1_arr[2],T1_arr[3],T1_arr[4])
    new_T2 = datetime.datetime(T2_arr[0],T2_arr[1],T2_arr[2],T2_arr[3],T2_arr[4])
    h = int((new_T2-new_T1).total_seconds()/60/60) # hours
    return h

def controlled_N_new(x,t,Nintmax,Nintmin,Vmax,Ks,dNextoutdt,dNextindt,miu,dmoutdt,Nintcrit,S,Z,KI,K0,Ka,Topt,
                Tmin,Tmax,losses20,teta,Sopt,Smin,Smax,n,umol_to_percent_DW,Temp,I0):


    Next = x[0] # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin)  # units: [ ]
    if Next <= 0:
        uN = 0
    else: 
        uN = Vmax * Next / (Ks + Next)  # units: [umol N/g DW/h]
    if Nint >= Nintcrit:
        fN = 1
    else:
        fN = ((Nint - Nintmin)/Nint) / ((Nintcrit - Nintmin)/Nintcrit)  # units: [ ]
    fP = 1  #(N:P < 12)


    # density - light penetration effects:
    SD = m * 5 / (0.2 * 0.178) # Stocking Density. units: [g DW/ m^2]
    I_average = (I0(t) / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
    fI = I_average / (I_average + KI) # units: [-]
    # print(f"t={t} I0={I0(t)} fI={fI}")
    
        # Temperature effects:         
    if Temp <= Topt:
        Tx = Tmin
    else:
        Tx = Tmax

    fT = np.exp(-2.3 * ((Temp - Topt) / (Tx - Topt))**n) # Temp from temperature data

    # S (salinity) effects
    if S < Sopt:
        Sx = Smin
        b = 2.5
        if S < 5:
            fS = ((S - Smin)/(Sopt - Sx))
        elif S >= 5:
            fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b           
    elif S >= Sopt:
        Sx = Smax
        b = 4.4 # found by solver in fs file
        fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b

    # empirically defined losses
    losses = losses20 * teta ** (Temp - 20)

    # limiting factors:
    g = min(fN,fI,fP) * fT * fS
    # print(f"t={t} g= {g} fN = {fN}, fI = {fI} fP = {fP} fT = {fT} fS = {fS}")

    
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    
    # Reactor Next -> Nint -> m and feedback

    dNextdt = - Neff * uN * m  - dNextoutdt * Next + dNextindt # [umol N/l/h] # added *Next
    dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * g # units: [%g N/g DW/h]

    if fI == 0:
        losses = 0
    else: 
        losses = losses20 * teta ** (Temp - 20)
    dmdt = (miu * g - losses) * m #units: [g DW/l/h]
    
    #dNextdt = -Neff * uN * m - dNextoutdt + dNextindt # units: [umol N/l/h]
    #dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * fN  #units: [%g N/g DW/h]
    #dmdt = miu * fN * m - dmoutdt #units: [g DW/l/h]
    
    return [dNextdt,dNintdt,dmdt]




def constant_N_new(x,t,Nintmax,Nintmin,Vmax,Ks,dNextoutdt,dNextindt,miu,dmoutdt,Nintcrit,S,Z,KI,K0,Ka,Topt,
                Tmin,Tmax,losses20,teta,Sopt,Smin,Smax,n,umol_to_percent_DW,fTemp0,fI0,VA):


    Next = x[0] # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin)  # units: [ ]
    if Next <= 0:
        uN = 0
    else: 
        uN = Vmax * Next / (Ks + Next)  # units: [umol N/g DW/h]
    if Nint >= Nintcrit:
        fN = 1
    else:
        fN = ((Nint - Nintmin)/Nint) / ((Nintcrit - Nintmin)/Nintcrit)  # units: [ ]
    fP = 1  #(N:P < 12)


    # density - light penetration effects:
    SD = m * VA # Stocking Density. units: [g DW/ m^2]
    #I_average = (I0 / (Ka * SD)) * (1 - np.exp(-(Ka * SD))) # units: [umul photons/(m^2 s)] #light measured in water, so no need in K0
    I0 = fI0(t)
    #print(I0)
    I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
    fI = I_average / (I_average + KI) # units: [-]
    
        # Temperature effects:         
    Temp = fTemp0(t)
    if Temp <= Topt:
        Tx = Tmin
    else:
        Tx = Tmax

    fT = np.exp(-2.3 * ((Temp - Topt) / (Tx - Topt))**n) # Temp from temperature data

    # S (salinity) effects
    if S < Sopt:
        Sx = Smin
        b = 2.5
        if S < 5:
            fS = ((S - Smin)/(Sopt - Sx))
        elif S >= 5:
            fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b           
    elif S >= Sopt:
        Sx = Smax
        b = 4.4 # found by solver in fs file
        fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b

    # empirically defined losses
    losses = losses20 * teta ** (Temp - 20)

    # limiting factors:
    g = min(fN,fI,fP) * fT * fS

    
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    
    # Reactor Next -> Nint -> m and feedback

    dNextdt = 0 # [umol N/l/h]
    dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]

    if fI == 0:
        losses = 0
    dmdt = (miu * g - losses) * m #units: [g DW/l/h]
        
    return [dNextdt,dNintdt,dmdt]



def controlled_N_normalized(x,t,Nintmax,Nintmin,Vmax,Ks,dNextoutdt,dNextindt,miu,dmoutdt,Nintcrit,S,Z,KI,K0,Ka,Topt,
                Tmin,Tmax,losses20,teta,Sopt,Smin,Smax,n,umol_to_percent_DW,Temp,I0):


    Next = x[0] # units: [umol N/l]
    Nint = x[1] # units: [% g N/g DW]
    m = x[2] # units: [g DW/l]
    
    Neff = (Nintmax - Nint)/(Nintmax - Nintmin)  # units: [ ]
    if Next <= 0:
        uN = 0
    else: 
        uN = Vmax * Next / (Ks + Next)  # units: [umol N/g DW/h]
    if Nint >= Nintcrit:
        fN = 1
    else:
        fN = ((Nint - Nintmin)/Nint) / ((Nintcrit - Nintmin)/Nintcrit)  # units: [ ]
    fP = 1  #(N:P < 12)


    # density - light penetration effects:
    SD = m * 5 / (0.2 * 0.178) # Stocking Density. units: [g DW/ m^2]
    I_average = (I0(t) / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
    fI = I_average / (I_average + KI) # units: [-]
    
        # Temperature effects:         
    if Temp <= Topt:
        Tx = Tmin
    else:
        Tx = Tmax

    fT = np.exp(-2.3 * ((Temp - Topt) / (Tx - Topt))**n) # Temp from temperature data

    # S (salinity) effects
    if S < Sopt:
        Sx = Smin
        b = 2.5
        if S < 5:
            fS = ((S - Smin)/(Sopt - Sx))
        elif S >= 5:
            fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b           
    elif S >= Sopt:
        Sx = Smax
        b = 4.4 # found by solver in fs file
        fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b

    # empirically defined losses
    losses = losses20 * teta ** (Temp - 20)

    # limiting factors:
    g = min(fN,fI,fP) * fT * fS

    
    umol_to_percent_DW = 100*14e-6 #[% g N/umol N] 
    
    # Reactor Next -> Nint -> m and feedback

    dNextdt = - Neff * uN * m  - dNextoutdt * Next + dNextindt # [umol N/l/h] # added *Next
    dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * g # units: [%g N/g DW/h]

    if fI == 0:
        losses = 0
    dmdt = (miu * g - losses) * m #units: [g DW/l/h]
    
    #dNextdt = -Neff * uN * m - dNextoutdt + dNextindt # units: [umol N/l/h]
    #dNintdt = umol_to_percent_DW * Neff * uN - Nint * miu * fN  #units: [%g N/g DW/h]
    #dmdt = miu * fN * m - dmoutdt #units: [g DW/l/h]
    
    return [dNextdt,dNintdt,dmdt]


def plot_result(t_model,Next_model,Nint_model,m_model,Nint=None, yerrNint=None, tNint=None,m=None, yerrm=None, tm=None, Next=None, yerrNext=None, tNext=None):
    """ Plot time series of the results """
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(12,10))
    ax[0].plot(t_model,Next_model,'r.--')
    if Next is not None:
        ax[0].errorbar(tNext[::],Next[::],yerr=yerrNext,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        ax[0].plot()
        # ax[0].plot(tNext[::],Next[::],'ko')
    ax[0].set_ylabel('Next \n [umol / l]', fontsize=12, weight='bold')
    # plt.axis([0, 180, 0, 430])
    ax[1].plot(t_model,Nint_model,'gx--')
    ax[1].set_ylabel('Nint \n [% g N / g DW]', fontsize=12, weight='bold')
    if Nint is not None:
        ax[1].errorbar(tNint[::],Nint[::],yerr=yerrNint,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        # Add Continuous Errors? https://jakevdp.github.io/PythonDataScienceHandbook/04.03-errorbars.html
        ax[1].plot()
        # ax[1].plot(tNint[::],Nint[::],'ko')
        
    # plt.axis([0, 180, 2.0, 4.0])
    ax[2].plot(t_model,m_model,'b+--')
    ax[2].set_xlabel('time [hours]')
    if m is not None:
        ax[2].errorbar(tm[::],m[::],yerr=yerrm,fmt='.',color='black',ecolor='lightgray',elinewidth=3, capsize=0)
        ax[2].plot()
        # ax[2].plot(tm[::],m[::],'ko')
    ax[2].set_ylabel('m \n [g DW $L^{-1}$]', fontsize=12, weight='bold')
    # plt.axis([0, 180, 0.03, 0.3])
    # plot_result(np.hstack(T),np.hstack(NEXT),np.hstack(NINT),np.hstack(M))
    # ax[0].plot(df.iloc[:][1::2],df.iloc[:][2::2],'--s')

    

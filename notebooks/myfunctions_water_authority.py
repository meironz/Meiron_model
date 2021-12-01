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

    
# Growth function with ex-situ (IMS) data and Reading parameters
def Bottles_IMS(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, f1,f0,dilution,n,umol_to_percent_DW,dNextoutdt):
    """
    This is a second version of Reading_val. The difference is that this function works with IMS data and not HOBO data
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    dydt = []
    
    # every 4 state variables are a single reactor
    n_reactors = 1
    
    Temp = f1(t)
    #Temp = 20 #https://isramar.ocean.org.il/isramar_data/TimeSeries.aspx
    I0 = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs
        
        # state variables: 
        
        i = 3*i_reactor # temporary counter of the state variables
        
        Next = y[i]  
        Nint = y[i+1]  
        m = y[i+2]
        
        # Nutrient consumption:
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
        SD = m * 0.001 * 0.085 # Stocking Density. units: [g DW/ m^2] #### Check bottle dimensions
        I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects: 
       
        if Temp <= Topt:
            Tx = Tmin
        else:
            Tx = Tmax
        
        n = 2
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
        
        # Reactor Next -> Nint -> m and feedback
        dNext = - Neff * uN * m - dNextoutdt * Next # [umol N/l/h]
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]

        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNext,dNint,dm])
        
    return dydt
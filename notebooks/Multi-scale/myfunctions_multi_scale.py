import matplotlib.pyplot as plt
import numpy as np

# growth function with HOBO data and Reading parameters
def Reading_val(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream,f1,f0,dilution,n,umol_to_percent_DW):
    """
    This is a multi-reactor, algae farm, simulation base function that will be 
    solved in time. The multiple reactors are built along the streamwise (x) direction and numbered
    in the order of appearance, i.e. x=0 is the first reactor that meets the original nutrient stream, Sea 
    Every following reactor will receive a different nutrient concentration, Sea(x) due to the dilution 
    by the previous reactor. In addition, we could introduce the dilution by the stream itself, if distances
    between the reactors are large for instance. 
    
    The main construction here is that each reactor is 4 ODEs: 
    1. dNsea/dt
    2. dNext/dt
    3. dNint/dt, and 
    4. dm/dt
    
    we develop a list of such using a loop and in the loop we create the coupling by propagating the Sea[i] 
    from the first reactor to the following ones.
    
    Nsea :  nutrient content in the stream of water, [umol N/l]
    Qsea : flow rate defined, [l/h]
    
    Qp : flow rate of the airlift pump [l/h], transfers Sea to the reactor to increase Next, 
    due to conservation, the same flow rate is overflows from the reactor back to the Qsea, reducing the
    concentration in the Qsea stream, thus affecting the following reactors
    
    Next : nutrient content in the reactor [umol N/l]
    
    Nint : nutrient content in the algae biomass [% g N/g DW]
    
    m : biomass weight [g DW/l]
    
    all the inputs are floating numbers
    
    *args: additional arguments:
    
    Nintcrit :  
    Nintmax : 
    Nintmin : 
    Vmax : 
    Ks :  
    KN :  
    miu : 
    S : 
    Z : 
    KI : 
    K0 : 
    Ka : 
    Topt : 
    Tmin : 
    Tmax : 
    losses20 : 
    teta : 
    Sopt : 
    Smin : 
    Smax :
    Qp[n] : n (number of reactors) - airlift pumps, constant flow rate, different per reactor
    Qsea : Stream flow rate in the sea
    Nsea_upstream : concentration of N upstream the reactors
   
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    dydt = []
    
    # let's figure out the number of reactors, every 4 state variables are a single reactor
    n_reactors = 1
       
    # print(f"{n_reactors} reactors")
    # print(f"{Qp} Qp")
    
    Temp = f1(t)
    
    I_average = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs

        # print(f'reactor #{i_reactor}')
        
        # state variables: 
        
        i = 4*i_reactor # temporary counter of the state variables
        
        Nsea = y[i]   
        Next = y[i+1]  
        Nint = y[i+2]  
        m = y[i+3]
        
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
        SD = m * 1.785 / 2 * 1000 # Stocking Density. units: [g DW/ m^2]
        #I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects: 
        #Temp = T_interp[time]
        
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

        #print('I_average = ' + str(I_average))

        
        # Qsea to reactor and back, per reactor, input is from the last one:
        if i_reactor > 0:
            Nsea_upstream = y[i-4] # y is a vector state, 4 state variables backwards is previous Nsea
            
        # print(i_reactor, Nsea,Nsea_upstream)
            
        #dNsea = (Qsea * (Nsea_upstream - Nsea) - Qp[i_reactor] * (Nsea - Next)) / (1.785 * 1000)
        dNsea = (Qsea * (Nsea_upstream * (1 - dilution) - Nsea) - Qp * (Nsea - Next)) / (1.785 * 1000)
        
        # Reactor Next -> Nint -> m and feedback
        #dNext = Qp[i_reactor] / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNext = Qp / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h] 
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]
        #if miu * g < losses:   #testing if limiting to positive or zero growth will help
        #    dm = 0
        #else:
        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNsea,dNext,dNint,dm])
        
    return dydt


# Growth function with IMS data and Reading parameters
def Reading_val_IMS(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream,f1,f0,dilution,n,umol_to_percent_DW):
    """
    This is a second version of Reading_val. The difference is that this function works with IMS data and not HOBO data
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    dydt = []
    
    # let's figure out the number of reactors, every 4 state variables are a single reactor
    n_reactors = 1
       
    # print(f"{n_reactors} reactors")
    # print(f"{Qp} Qp")
    
    Temp = f1(t)
    #Temp = 20 #https://isramar.ocean.org.il/isramar_data/TimeSeries.aspx
    I0 = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs
        # print(f'reactor #{i_reactor}')
        
        # state variables: 
        
        i = 4*i_reactor # temporary counter of the state variables
        
        Nsea = y[i]   
        Next = y[i+1]  
        Nint = y[i+2]  
        m = y[i+3]
        
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
        SD = m * 1.785 / 2 * 1000 # Stocking Density. units: [g DW/ m^2]
        #I0 = I[time]
        I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects: 
        #Temp = T_interp[time]
        
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
        
        #print('I0 = ' + str(I0))
        #print('I_average = ' + str(I_average))
        
        # Qsea to reactor and back, per reactor, input is from the last one:
        if i_reactor > 0:
            Nsea_upstream = y[i-4] # y is a vector state, 4 state variables backwards is previous Nsea
            
        # print(i_reactor, Nsea,Nsea_upstream)
            
        #dNsea = (Qsea * (Nsea_upstream - Nsea) - Qp[i_reactor] * (Nsea - Next)) / (1.785 * 1000)
        dNsea = (Qsea * (Nsea_upstream * (1 - dilution) - Nsea) - Qp * (Nsea - Next)) / (1.785 * 1000)
        
        # Reactor Next -> Nint -> m and feedback
        #dNext = Qp[i_reactor] / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNext = Qp / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]
        #if miu * g < losses:   #testing if limiting to positive or zero growth will help
        #    dm = 0
        #else:
        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNsea,dNext,dNint,dm])
        
    return dydt

# Growth function with IMS data and Reading parameters
def Reading_val_IMS_high(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream,f1,f0,dilution,n,umol_to_percent_DW):
    """
    This is a second version of Reading_val. The difference is that this function works with IMS data and not HOBO data
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    dydt = []
    
    # let's figure out the number of reactors, every 4 state variables are a single reactor
    n_reactors = 1
       
    # print(f"{n_reactors} reactors")
    # print(f"{Qp} Qp")
    
    Temp = f1(t)
    #Temp = 20 #https://isramar.ocean.org.il/isramar_data/TimeSeries.aspx
    I0 = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs
        # print(f'reactor #{i_reactor}')
        
        # state variables: 
        
        i = 4*i_reactor # temporary counter of the state variables
        
        Nsea = y[i]   
        Next = y[i+1]  
        Nint = y[i+2]  
        m = y[i+3]
        
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
        SD = m * 1.785 / 2 * 1000 # Stocking Density. units: [g DW/ m^2]
        #I0 = I[time]
        I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects: 
        #Temp = T_interp[time]
        
        if Temp <= Topt:
            Tx = Tmin
        else:
            Tx = Tmax

        fT = np.exp(-2.3 * ((Temp - Topt) / (Tx - Topt))**n) # Temp from temperature data

        # S (salinity) effects
        fS = 1
        #if S < Sopt:
        #    Sx = Smin
        #    n = 2.5
        #    if S < 5:
        #        fS = ((S - Smin)/(Sopt - Sx))
        #    elif S >= 5:
        #        fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** n           
        #elif S >= Sopt:
        #    Sx = Smax
        #    n = 4.4 # found by solver in fs file
        #    fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** n

        # empirically defined losses
        losses = losses20 * teta ** (Temp - 20)


        # limiting factors:
        g = min(fN,fI,fP) * fT * fS
        
        #print('I0 = ' + str(I0))
        #print('I_average = ' + str(I_average))
        
        # Qsea to reactor and back, per reactor, input is from the last one:
        if i_reactor > 0:
            Nsea_upstream = y[i-4] # y is a vector state, 4 state variables backwards is previous Nsea
            
        # print(i_reactor, Nsea,Nsea_upstream)
            
        #dNsea = (Qsea * (Nsea_upstream - Nsea) - Qp[i_reactor] * (Nsea - Next)) / (1.785 * 1000)
        dNsea = (Qsea * (Nsea_upstream * (1 - dilution) - Nsea) - Qp * (Nsea - Next)) / (1.785 * 1000)
        
        # Reactor Next -> Nint -> m and feedback
        #dNext = Qp[i_reactor] / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNext = Qp / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]
        #if miu * g < losses:   #testing if limiting to positive or zero growth will help
        #    dm = 0
        #else:
        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNsea,dNext,dNint,dm])
        
    return dydt


# let's prepare a simple prototype and explain the concept
def multi_N_f_un(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream,f1,f0,dilution,n,umol_to_percent_DW):
    """
    This is a multi-reactor, algae farm, simulation base function that will be 
    solved in time. The multiple reactors are built along the streamwise (x) direction and numbered
    in the order of appearance, i.e. x=0 is the first reactor that meets the original nutrient stream, Sea 
    Every following reactor will receive a different nutrient concentration, Sea(x) due to the dilution 
    by the previous reactor. In addition, we could introduce the dilution by the stream itself, if distances
    between the reactors are large for instance. 
    
    The main construction here is that each reactor is 4 ODEs: 
    1. dNsea/dt
    2. dNext/dt
    3. dNint/dt, and 
    4. dm/dt
    
    we develop a list of such using a loop and in the loop we create the coupling by propagating the Sea[i] 
    from the first reactor to the following ones.
    
    Nsea :  nutrient content in the stream of water, [umol N/l]
    Qsea : flow rate defined, [l/h]
    
    Qp : flow rate of the airlift pump [l/h], transfers Sea to the reactor to increase Next, 
    due to conservation, the same flow rate is overflows from the reactor back to the Qsea, reducing the
    concentration in the Qsea stream, thus affecting the following reactors
    
    Next : nutrient content in the reactor [umol N/l]
    
    Nint : nutrient content in the algae biomass [% g N/g DW]
    
    m : biomass weight [g DW/l]
    
    all the inputs are floating numbers
    
    *args: additional arguments:
    
    Nintcrit :  
    Nintmax : 
    Nintmin : 
    Vmax : 
    Ks :  
    KN :  
    miu : 
    S : 
    Z : 
    KI : 
    K0 : 
    Ka : 
    Topt : 
    Tmin : 
    Tmax : 
    losses20 : 
    teta : 
    Sopt : 
    Smin : 
    Smax :
    Qp[n] : n (number of reactors) - airlift pumps, constant flow rate, different per reactor
    Qsea : Stream flow rate in the sea
    Nsea_upstream : concentration of N upstream the reactors

    
    
    """
    
    #print( t, y)
    # print(Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
    #        Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream)
    dydt = []
    
    # let's figure out the number of reactors, every 4 state variables are a single reactor
    n_reactors = len(y)//4
       
    # print(f"{n_reactors} reactors")
    # print(f"{Qp} Qp")
    
    Temp = f1(t)
    I0 = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs
        #time = times[index]
        # print(f'reactor #{i_reactor}')
        
        # state variables: 
        
        i = 4*i_reactor # temporary counter of the state variables
        
        Nsea = y[i]   
        Next = y[i+1]  
        Nint = y[i+2]  
        m = y[i+3]
        
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
        SD = m * 1.785 / 2 * 1000 # Stocking Density. units: [g DW/ m^2]
        #I0 = I[time]
        I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
        #I_average = 80 #cancle after I is imported from file
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects: 
        #Temp = T_interp[time]
        
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
        
        #print('fI = ' + str(fI))
        #print('g = ' + str(g))
        #print('I0 = ' + str(I0))
        
        # Qsea to reactor and back, per reactor, input is from the last one:
        if i_reactor > 0:
            Nsea_upstream = y[i-4] # y is a vector state, 4 state variables backwards is previous Nsea
            
        # print(i_reactor, Nsea,Nsea_upstream)
            
        #dNsea = (Qsea * (Nsea_upstream - Nsea) - Qp[i_reactor] * (Nsea - Next)) / (1.785 * 1000)
        dNsea = (Qsea * (Nsea_upstream * (1 - dilution) - Nsea) - Qp * (Nsea - Next)) / (1.785 * 1000)
        
        # Reactor Next -> Nint -> m and feedback
        #dNext = Qp[i_reactor] / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNext = Qp / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]
        #if miu * g < losses:   #testing if limiting to positive or zero growth will help
            #dm = 0
        #else:
        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNsea,dNext,dNint,dm])
        #index = index + 1
        # print(dydt)
        
    return dydt

# plotting NSEA, NEXT, NINT and M vs time - each reactor (or 10th reactor) gets a color
def plot_result_un(T,NSEA,NEXT,NINT,M,resolution,n_reactors,Nsea):
    """ Plot time series of the results, it's organized as a list of solution structures, sol.t,sol.y 
    
    We will plot 4 axes each for the separate quantity, every reactor gets its line
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig,ax = plt.subplots(4,1,sharex=True,figsize=(12,10))
    
    xlabels = ['day 0', 'day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7', 'day 8', 'day 9','day 10', 'day 11',
               'day 12','day 13', 'day 14']
    t = T[0]
    I = []
    
    for i in range(0,n_reactors,resolution):
        ax[0].plot(t,NSEA[i],'.',markersize=2)
        ax[0].set_ylim([0,Nsea*1.2])
        ax[1].plot(t,NEXT[i],'.', markersize=2)
        ax[1].set_ylim([0,Nsea*1.2])
        ax[2].plot(t,NINT[i],'.', markersize=2)
        ax[2].set_ylim([0.5,3.5])
        ax[3].plot(t,M[i],'.', markersize=2)
        I.append(i)
        
    ax[0].set_ylabel('Nsea \n [umol / l]',fontsize=10, weight="bold")
    ax[1].set_ylabel('Next \n [umol / l]',fontsize=10, weight="bold")
    ax[2].set_ylabel('Nint \n [% g N / g DW]',fontsize=10, weight="bold")
    ax[3].set_ylabel('m \n [g DW / l]',fontsize=10, weight="bold")
    ax[3].set_xlabel('Time',fontsize=14, weight="bold")
    
    ax[3].set_xlabel('time [hours]')
    
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    ax[2].set_xticklabels([])
    ax[3].set_xticklabels([])
    ax[0].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[1].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[2].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[3].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[3].set_xlabel('Time',fontsize=14, weight="bold")
    ax[3].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=10, weight="bold")
    
    ax[0].legend(I)
### A-plots
# TRphasediagram.py

import pickle

############ OMEGA PLOTTING FUNCTION ######################
def plotomegaPT(hug,label=False,branch=1):
    omegapmin=85 # Erik's recommendation for minimum real data
    omegapmax=111
    #usmax=9.5
    usmin = 9.23
    usmax = 9.5
    mss = 5
    aaa=0.35
    if not((branch==4) or (branch==5) or (branch==7)):
        col9='darkgoldenrod'
        aaa=0.2
        zzz = np.where((hug.P > omegapmin) & (hug.us < usmin))[0]
        if label:
            plt.plot(hug.P[zzz[0]],hug.T[zzz[0]],'.',color=col9,markersize=mss,alpha=aaa)#,label='Decaying shock-solid branch (Omega)')
            plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)#,label='Decay: Transparent front')#,label='Solid branch (Omega)')
            plt.plot(hug.P[zzz[0]],hug.T[zzz[0]],'.',color=col9,markersize=mss,alpha=1,label='Decaying Transparent front')
        else:
            plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)

    if (branch==1) or (branch==2) or (branch==7):
        col9='grey'
        zzz = np.where((hug.us > usmin) & (hug.us < usmax))[0]
        if label:
            plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)#,label='Decay: Transition')
            plt.plot(hug.P[zzz[0]],hug.T[zzz[0]],'.',color=col9,markersize=mss,alpha=1,label='Decaying Transition')
            #plt.plot(hug.P[zzz[0]],hug.T[zzz[0]],'.',color=col9,markersize=6,alpha=aaa,label='Shock front becomes opaque (Omega)')
            #plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)#,label='Shock front becomes opaque (Omega)')
        else:
            plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)
    if (branch==1) or (branch==4) or (branch==5) or (branch==6) or (branch==7):
        if not(branch==5 or branch==6):
            #mss = 2
            col9='yellowgreen'
#            col9='green'
            #zzz = np.where((hug.us >= usmax) & (hug.rho <=5.2))[0]
            zzz = np.where((hug.us >= usmax) & (hug.T <= 6200))[0]
            if label:
                plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)#,label='Decay: Opaque front')#,label='Opaque liquid shock front layer (Omega)')
                plt.plot(hug.P[zzz[0]],hug.T[zzz[0]],'.',color=col9,markersize=mss,alpha=1,label='Decaying Opaque front')#,label='Opaque liquid shock front layer (Omega)')
            else:
                plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)#,label='Disordered phases')
        if not(branch==7):
            rhotrans=5.2
            #mss = 2
            col9='green'
            #zzz = np.where(hug.rho > rhotrans)[0]
            zzz = np.where(hug.T > 6200)[0]
            if label:
                plt.plot(hug.P[zzz[0]],hug.T[zzz[0]],'.',color=col9,markersize=mss,alpha=1,label='Decaying Liquid Field')#,label='Decaying shock-liquid branch (Omega)')
                plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)#,label='Liquid branch (Omega)')
            else:
                plt.plot(hug.P[zzz],hug.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)#,label='Liquid')

############ MAIN PLOTS ######################

#for pind in range(28):
if 1:
    fig1=False
    fig2=False
    fig3=False
    fig4=False
    fig5=False
    fig6=False  # with release isentropes
    fig15=False
    fig7=False
    fig8=False # Kaitlyn data presentation
    fig9=False # Kaitlyn data presentation
    fig10=False # Luo et al. 2002 comparison
    fig11=False # Luo et al. 2003 comparison
    fig12=False # very high pressure melt curve
    fig13=False # compare to Andrault discussion point
    fig14=False # compare to Andrault discussion point
    fig15=False # compare to Andrault discussion point
    fig16=False # compare to Andrault discussion point
    fig17=False # quartz discussion
    fig18=False # FS discussion
    fig19=False
    fig20=False
    fig21=False
    fig22=False
    fig23=False
    fig24=False
    fig25=False
    fig26=False
    fig27=False
    fig28=False
    fig29=False
    fig30=False
    figtest=False
    if pind == 1:
        fig1=True # ab initio and static phase boundary data
    if pind == 2:
        fig15=True
    if pind == 3:
        fig2=True
    if pind == 4:
        fig3=True
    if pind == 5:
        fig4=True
    if pind == 6:
        fig5=True
    if pind == 7:
        fig6=True
    if pind == 8:
        fig7=True
    if pind == 9: # Summarize new regions on the fused silica Hugoniot
        fig8=True
    if pind == 10: # Summarize new regions on the Quartz Hugoniot
        fig9=True
    if pind == 11: # Luo et al. 2002 MD melt curve comparison
        fig10=True
    if pind == 12: # Luo et al. 2003 stv model comparison
        fig11=True
        fig3=True
    if pind == 13: # expand PT to high pressures with stv shock temp data LINEAR
        fig12=True  
        fig13=True
    if pind == 14: # expand PT to high pressures with stv shock temp data LOG
        fig12=True
        fig14=True
    if pind == 15: # expand PT to high pressures with stv shock temp data LOG
        fig9=True
        fig8=True
        fig13=True # add comparison to Andrault for discussion point
    if pind == 16: # plot omega correction
        fig3=True
    if pind == 17: # quartz discussion
        fig17=True
        fig21=True
    if pind == 18: # FS discussion
        fig18=True
    if pind == 19:
        fig19=True
    if pind ==20:
        fig20=True
    if pind==21: # metastable PT space
        fig20=True
    if pind == 22: # quartz discussion Z proposal
        fig17=True
        fig21=True
        fig22=True
    if pind == 23: # quartz paper fig1
        fig1=True
        fig23=True
    if pind == 24: # quartz paper fig2 -qtz
        fig3=True
        fig16=True
        fig24=True
    if pind == 124: # quartz paper fig2 -qtz
        fig3=True
        fig16=True
        fig24=True
    if pind == 25: # quartz paper fig 2-fs
        fig8=True
        fig25=True
    if pind == 125: # quartz paper fig 2-fs
        fig8=True
        fig25=True
    if pind == 26: # quartz transition front
        #fig17=True # omega colored data
        #fig21=True
        #fig22=True
        fig26=True
    if pind == 27: # fs transition front
        fig8=True # fs data
        #fig21=True
        #fig25=True
        fig27=True
    if pind == 28:
        #fig17=True # omega colored data
        #fig21=True
        #fig22=True
        #fig26=True
        fig28=True
    if pind == 29:
        #fig28=True
        fig29=True
    if pind == 30:
        fig30=True
    if pind == 2425:
        fig3=True
        fig16=True
        fig24=True       
        fig8=True
        fig25=True
        figtest=True
    if pind == 2627:
        fig26=True
        fig8=True #
        fig27=True        
    #########################################
    # default figure shape
    plt.rcParams.update(plt.rcParamsDefault)
    fig = plt.figure(figsize=(10, 6))
    plt.rcParams.update({'font.size': 14})
    if not(fig6):
        fig = plt.figure(figsize=(10, 6))
        plt.rcParams.update({'font.size': 12})
    else:
        fig = plt.figure(figsize=(10, 8))  
        plt.rcParams.update({'font.size': 12})
    if fig20:
        fig = plt.figure(figsize=(8, 6))
        plt.rcParams.update({'font.size': 14})
    if fig21:
        fig = plt.figure(figsize=(9, 6))
        plt.rcParams.update({'font.size': 14})
    # MARCH 2025 REVISED FIGURES
    if fig23 or fig24 or fig25 or fig26 or fig27 or fig28 or fig29:
        fig = plt.figure(figsize=(10, 6))
        plt.rcParams.update({'font.size': 16})
    if fig23:
        fig = plt.figure(figsize=(10, 6))
        plt.rcParams.update({'font.size': 14})
    #if fig29:
    #    fig = plt.figure(figsize=(10, 6))
    #    plt.rcParams.update({'font.size': 26})
    if fig30:
        fig = plt.figure(figsize=(8, 6))
        plt.rcParams.update({'font.size': 14})  
    if pind==124:
        fig = plt.figure(figsize=(6, 6))
        plt.rcParams.update({'font.size': 14})  
    if pind==125:
        fig = plt.figure(figsize=(6, 6))
        plt.rcParams.update({'font.size': 14})  
    #########################################    
    # colors for solid-solid phase boundaries
    pbw=4
    pbcol='#C5C9C7'

    if not(fig12):
        # Gasparik phase diagram
        plt.plot(eosdat_list[G14mc['id']].P,eosdat_list[G14mc['id']].T,color="black")#,label='Low-P Phases')
        plt.plot(eosdat_list[G14lhqtz['id']].P,eosdat_list[G14lhqtz['id']].T,color="black")
        plt.plot(eosdat_list[G14trcr['id']].P,eosdat_list[G14trcr['id']].T,color="black")
        plt.plot(eosdat_list[G14crstv['id']].P,eosdat_list[G14crstv['id']].T,color="black")
        plt.plot(eosdat_list[G14trhq['id']].P,eosdat_list[G14trhq['id']].T,color="black")
        plt.plot(eosdat_list[G14crhq['id']].P,eosdat_list[G14crhq['id']].T,color="black")
        plt.plot(eosdat_list[G14lqcs['id']].P,eosdat_list[G14lqcs['id']].T,color="black")
        plt.plot(eosdat_list[G14hqcs['id']].P,eosdat_list[G14hqcs['id']].T,color="black")

    # plot background solid-solid phase boundaries GGA on most plots
    if 0:
        # plot GGA boundaries
        # stv-CaCl2 T04
        plt.plot(eosdat_list[T04pba['id']].P,eosdat_list[T04pba['id']].T,'-',color=pbcol,linewidth=pbw)
        plt.plot([np.max(eosdat_list[T04pba['id']].P), Pt2],[np.max(eosdat_list[T04pba['id']].T),Tt2],':',color=pbcol,linewidth=pbw)
        # CaCl2 -> alpha-PbO2
        plt.plot(eosdat_list[T04pbb['id']].P,eosdat_list[T04pbb['id']].T,'-',color=pbcol,linewidth=pbw)
        plt.plot([np.max(eosdat_list[T04pbb['id']].P), Pt3],[np.max(eosdat_list[T04pbb['id']].T),Tt3],':',color=pbcol,linewidth=pbw)
        zzz = np.where(eosdat_list[T04pbb['id']].T > 1000.)[0]
        res = np.polyfit(np.flip(eosdat_list[T04pbb['id']].P[zzz]),np.flip(eosdat_list[T04pbb['id']].T[zzz]),1)
        print('RES b=',res)
        print('4000 K pressure = ',np.interp(4000,np.flip(eosdat_list[T04pbb['id']].T[zzz]),np.flip(eosdat_list[T04pbb['id']].P[zzz])))
        #plt.plot(np.arange(110,170),res[1]+res[0]*np.arange(110,170),'k-')
        zzz = np.where(eosdat_list[T04pba['id']].T > 1000.)[0]
        res = np.polyfit(np.flip(eosdat_list[T04pba['id']].P[zzz]),np.flip(eosdat_list[T04pba['id']].T[zzz]),1)
        print('RES a=',res)
        print('4000 K pressure = ',np.interp(4000,np.flip(eosdat_list[T04pba['id']].T[zzz]),np.flip(eosdat_list[T04pba['id']].P[zzz])))
        #plt.plot(np.arange(60,100),res[1]+res[0]*np.arange(60,100),'k-')
        
        # alpha-PbO2->pyrite LDA from O05
 #       ind = np.where(eosdat_list[P25pbscan['id']].T<7250)[0]
        #plt.errorbar(eosdat_list[P25pbscan['id']].P[ind],eosdat_list[P25pbscan['id']].T[ind],fmt='-',color=pbcol,linewidth=pbw)
        ind = np.where(eosdat_list[P25pb['id']].T<6250)[0]
        plt.errorbar(eosdat_list[P25pb['id']].P[ind],eosdat_list[P25pb['id']].T[ind],fmt='-',color=pbcol,linewidth=pbw)
        #plt.errorbar(eosdat_list[O05pbc['id']].P,eosdat_list[O05pbc['id']].T,fmt=':',color=pbcol,linewidth=pbw)
    if not(fig23 or fig1):
        # plot Deng et al. 2026 D26
        colO05='teal'
        # add Deng et al. 2026 PBEsol boundaries
        # stv-CaCl2
        ptmp = np.arange(109-85)+85
        ttmp = 91.45*ptmp -3999
        #plt.plot(ptmp,ttmp,'--',color=colO05,label='D26 stv-CaCl₂ PBEsol')
        ind = np.where(eosdat_list[D26pb['id']].phase_solid =='stishovite/post-stishovite boundary (Eq.5)')
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'-',color=pbcol,linewidth=pbw)
        
        # CaCl2-siefertite
        ptmp = np.arange(129-111)+110
        ttmp = 132.98*ptmp -11011
        #plt.plot(ptmp,ttmp,'-.',color=colO05,label='D26 CaCl₂-pyrite PBEsol')
        ind = np.where(eosdat_list[D26pb['id']].phase_solid =='post-stishovite/seifertite boundary (Eq.5)')
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'-',color=pbcol,linewidth=pbw)
        
        # alpha-PbO2->pyrite LDA from O05
 #       ind = np.where(eosdat_list[P25pbscan['id']].T<7250)[0]
        #plt.errorbar(eosdat_list[P25pbscan['id']].P[ind],eosdat_list[P25pbscan['id']].T[ind],fmt='-',color=pbcol,linewidth=pbw)
        ind = np.where(eosdat_list[P25pb['id']].T<6250)[0]
        plt.errorbar(eosdat_list[P25pb['id']].P[ind],eosdat_list[P25pb['id']].T[ind],fmt='-',color=pbcol,linewidth=pbw)
        #plt.errorbar(eosdat_list[O05pbc['id']].P,eosdat_list[O05pbc['id']].T,fmt=':',color=pbcol,linewidth=pbw)


    #########################################   
    # plot proposed new melt curve with triple points
    if 1: #not(fig1 or fig7 or fig10):
        plt.plot(stsmcsmooth.P,stsmcsmooth.T,'m-',linewidth=3,label='Model melt curve')
        plt.fill_between(stsmcsmooth.P,stsmcsmooth.T-stsmcsmooth.T_std,stsmcsmooth.T+stsmcsmooth.T_std,color='magenta',alpha=0.1)
    if 0: #fig24 or fig25:
        # Plot Andrault et al. 2022 triple points
        plt.plot([Pt1],[Tt1],'mD',label='Est. 1st order triple pt.')
        plt.plot([Pt3],[Tt3],'mD')#,mfc='white') # estimate from extrapolated T04 GGA boundary; same value as Belonoshko et al. 1995 calculated point
        # Plot Andrault et al. 2022 triple points
        plt.plot([Pt2, Ptm],[Tt2,Ttm],'mD',mfc='white',label='Est. 2nd order triple pt.')
        #plt.plot([Pt3],[Tt3],'mD')#,mfc='white') # estimate from extrapolated T04 GGA boundary; same value as Belonoshko et al. 1995 calculated point
        #plt.plot([141],[6072],'mD',mfc='white') 
        #plt.plot([145,145],[0,10000])
    #if 1:
        # Plot Andrault et al. 2022 triple points
       # plt.plot([Pt1],[Tt1],'mD')
        #plt.plot([Pt3],[Tt3],'mD')#,mfc='white') # estimate from extrapolated T04 GGA boundary; same value as Belonoshko et al. 1995 calculated point
        # Plot Andrault et al. 2022 triple points 2nd order
        #plt.plot([Pt2, Ptm],[Tt2,Ttm],'mD',mfc='white')


    if fig30:
        # Plot Zhang et al. JGR 1993 cst-stv-L triple point
        plt.plot([13.7],[2800+273],'P',color='k',markersize=8,label='Z93 triple pt.')
        # Plot Andrault et al. 2022 triple points 
        plt.plot([15],[2900],'D',color='orange',mec='black',label='A22 triple pt.')
        plt.plot([40,85],[4800,5200],'D',color='orange',mfc='white',label='A22 2nd order triple pt.')
        # Plot Andrault et al. 2022 triple points
        #plt.plot([Pt1],[Tt1],'mD',label='Est. 1st order triple pt.')
       # estimate from extrapolated T04 GGA boundary; same value as Belonoshko et al. 1995 calculated point
#        plt.plot([Pt3],[Tt3],'mD',mec='k',label='Est. triple pt.') # estimate from extrapolated T04 GGA boundary; same value as Belonoshko et al. 1995 calculated point
        # Plot Andrault et al. 2022 triple points
#        plt.plot([Pt2],[Tt2],'mD',mfc='white',label='Est. 2nd order triple pt.')
        #plt.plot([Ptm],[Ttm],'mD',mfc='white',label='Est. 2nd order triple pt.')
        #plt.plot([Pt3],[Tt3],'mD')#,mfc='white') # estimate from extrapolated T04 GGA boundary; same value as Belonoshko et al. 1995 calculated point
        #plt.plot([141],[6072],'mD',mfc='white') 
        #plt.plot([145,145],[0,10000])
    #if 1:
        # Plot Andrault et al. 2022 triple points
       # plt.plot([Pt1],[Tt1],'mD')
        #plt.plot([Pt3],[Tt3],'mD')#,mfc='white') # estimate from extrapolated T04 GGA boundary; same value as Belonoshko et al. 1995 calculated point
        # Plot Andrault et al. 2022 triple points 2nd order
        #plt.plot([Pt2, Ptm],[Tt2,Ttm],'mD',mfc='white')

        # plot Deng et al. 2026 D26
        colO05='teal'
        # add Deng et al. 2026 PBEsol boundaries
        # stv-CaCl2
        tmpind = np.where(eosdat_list[D26pb['id']].phase_solid =='stishovite/post-stishovite boundary (Eq.5)')
        ind = np.max(tmpind)
        print(ind)
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'mD',mfc='white',label='D26 2nd order triple pt.')
        
        # CaCl2-siefertite
        tmpind = np.where(eosdat_list[D26pb['id']].phase_solid =='post-stishovite/seifertite boundary (Eq.5)')
        ind = np.max(tmpind)
        print(ind)
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'mD',label='D26 triple pt.')
        # sft-pyr Park et al. 2025
        plt.plot([165],[6370],'s',color='orange',mec='k',label='P26 triple pt.') #
    
    if 0:
        #plot Dong MC
        dong_dat = np.loadtxt('../data_silica_db/dong-mc.csv',skiprows=0,delimiter=',')
        plt.plot(dong_dat[:,0],dong_dat[:,1],'--',color='black',label='D25 MC')
        plt.plot(99.2,5867,'kD')
        #plt.plot([94,94],[0,8000])
        #plt.plot(63,4800,'cD',label='McQueen')
    #########################################   
    # draw experimental wedges for melt curve from 
    if not(fig1 or fig8 or fig9 or fig10 or fig12 or fig17 or fig18 or fig19 or fig20 or fig21 or fig23 or fig24 or fig25 or fig26 or fig28 or fig29):
        wcol = "#653700"
        wcol2 = "#FBDD7E"
        wcol3 = "#F97306"
        wop = 0.5
        wop2=0.75
        wop3=0.25
        ind = np.arange(0,11)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],label='melt curve exp',alpha=wop,color=wcol)
        ind = np.arange(12,16)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop,color=wcol)
        ind = np.arange(17,21)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop,color=wcol)
        #ind = np.arange(22,24)
        #print(stsmc.P[ind])
        #plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop2,color=wcol2,label='melt curve calc',ec='grey')
        #ind = np.arange(25,28)
        #print(stsmc.P[ind])
        #plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop2,color=wcol2,ec='grey')
    if fig13 or fig22:
        wcol = "#653700"
        wcol2 = "#FBDD7E"
        wcol3 = "#F97306"
        wop = 0.5
        wop2=0.75
        wop3=0.25
        ind = np.arange(0,11)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],label='S95+A22 melt curve',alpha=wop,color=wcol)
        ind = np.arange(12,16)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop,color=wcol)
        ind = np.arange(17,21)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop,color=wcol)
        ind = np.arange(22,24)
        #print(stsmc.P[ind])        
    if 0: #fig24 or fig25: # or fig26 or fig27: # 
        wcol = "#653700"
        wcol2 = "#FBDD7E"
        wcol3 = "#F97306"
        wop = 0.5
        wop2=0.75
        wop3=0.25
        ind = np.arange(0,11)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],label='S95+A22 melt curve',alpha=wop,color=wcol)
        ind = np.arange(12,16)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop,color=wcol)
        ind = np.arange(17,21)
        #print(stsmc.P[ind])
        plt.fill_between(stsmc.P[ind],stsmc.T_s[ind],stsmc.T_l[ind],alpha=wop,color=wcol)
        ind = np.arange(22,24)

    #########################################   
    # DAC melt curve points
    if fig1 or fig15 or fig10 or fig23:
        # stv melt curve Shen and Lazer 1995
        plt.errorbar(eosdat_list[S95mcstv['id']].P,eosdat_list[S95mcstv['id']].T,\
                     xerr=eosdat_list[S95mcstv['id']].P_std,yerr=eosdat_list[S95mcstv['id']].T_std,\
                     fmt='d',color='mediumblue',mfc='grey',markersize=5,label=S95mcstv['lab'])
        # Andrault et al. 2022 revised
        plt.errorbar(eosdat_list[A22mc['id']].P,eosdat_list[A22mc['id']].T,yerr=eosdat_list[A22mc['id']].T_std,\
                     fmt='d',mec='blueviolet',mfc='grey',markersize=5,color='blueviolet',label=A22mc['lab'])
        # do not plot old Andrault et al. 2020 points
        #plt.plot(eosdat_list[A20mc['id']].P,eosdat_list[A20mc['id']].T,'cd',mfc='orange',label=A20mc['lab'])

    #########################################   
    # Ab initio melt curve points\
    if fig12:
        # Gonzales et al. 2016 very high rpessure melt curve points
        plt.errorbar(eosdat_list[G16mc['id']].P[1:4],eosdat_list[G16mc['id']].T[1:4],fmt='d--',mfc='white',yerr=eosdat_list[G16mc['id']].T_std[1:4],label=G16mc['lab']) # pyrite structure
        plt.errorbar(eosdat_list[G16mc['id']].P[5:8],eosdat_list[G16mc['id']].T[5:8],fmt='s--',mfc='white',yerr=eosdat_list[G16mc['id']].T_std[5:8],label=G16mcb['lab']) # cotunnite structure
        plt.errorbar(eosdat_list[G16mc['id']].P[8:11],eosdat_list[G16mc['id']].T[8:11],fmt='o--',mfc='white',yerr=eosdat_list[G16mc['id']].T_std[8:11],label=G16mcc['lab']) # Fe2P structure

    if fig15:
        ind = np.where(eosdat_list[U10mcl['id']].P_liquid > 120)[0]
        plt.plot(eosdat_list[U10mcl['id']].P_liquid[ind],eosdat_list[U10mcl['id']].T_liquid[ind],'k+',label=U10mcl['lab'])
        ind = np.where(eosdat_list[U10mcl['id']].P_solid > 120)[0]
        plt.plot(eosdat_list[U10mcs['id']].P_solid[ind],eosdat_list[U10mcs['id']].T_solid[ind],'kx',label=U10mcs['lab'])
        plt.plot(eosdat_list[G16mc['id']].P[1:2],eosdat_list[G16mc['id']].T[1:2],'x',color='grey',label=G16mc['lab']) # pyrite structure
        #plt.errorbar(eosdat_list[G16mc['id']].P[0:4],eosdat_list[G16mc['id']].T[0:4],fmt='yd--',mfc='white',yerr=eosdat_list[G16mc['id']].T_std[0:4],label=G16mc['lab']) # pyrite structure
        plt.plot(eosdat_list[G24mcl['id']].P_liquid,eosdat_list[G24mcl['id']].T_liquid,'r+',label=G24mcl['lab'])
        plt.plot(eosdat_list[G24mcs['id']].P_solid,eosdat_list[G24mcs['id']].T_solid,'rx',label=G24mcs['lab'])
       
    if fig1 or fig10 or fig23: 
        # melt curve data
        plt.errorbar(eosdat_list[B95mc['id']].P,eosdat_list[B95mc['id']].T,fmt='X--',mfc='white',mec='sienna',color='sienna',yerr=eosdat_list[B95mc['id']].T_std,label=B95mc['lab'])
        plt.errorbar(eosdat_list[G16mc['id']].P[0:4],eosdat_list[G16mc['id']].T[0:4],fmt='d--',mfc='white',mec='sienna',color='sienna',yerr=eosdat_list[G16mc['id']].T_std[0:4],label=G16mc['lab']) # pyrite structure


        plt.plot(eosdat_list[D26mc['id']].P,eosdat_list[D26mc['id']].T,'kD',mfc='white',markersize=4,label=D26mc['lab']) # Deng et al. 2026

        #plt.errorbar(eosdat_list[G24mcl['id']].P,eosdat_list[G24mcl['id']].T,fmt='cd--',mfc='white',label=G24mcl['lab'])
        
        if fig10:
            # Luo data are offset from the rest - classical MD is not sufficient
            plt.errorbar(eosdat_list[L02mcstv['id']].P,eosdat_list[L02mcstv['id']].T,fmt='rs--',mfc='white',yerr=eosdat_list[L02mcstv['id']].T_std,label=L02mcstv['lab']+' 2phase')
            plt.errorbar(eosdat_list[L02mcstv['id']].P,eosdat_list[L02mcstv['id']].T_l,fmt='rs--',mfc='grey',label=L02mcstv['lab']+' 1phase')

        #plt.plot(eosdat_list[A22mc['id']].P,eosdat_list[A22mc['id']].T,'cd',label=A22mc['lab'])
        #plt.plot(eosdat_list[A22mc['id']].P,eosdat_list[A22mc['id']].T_liquid,'c+',label='A22 liquid')
        #plt.plot(eosdat_list[A22mc['id']].P,eosdat_list[A22mc['id']].T_solid,'cx',label='A22 solid')

        #plt.plot(eosdat_list[A22mcp['id']].P_liquid,eosdat_list[A22mcp['id']].T_liquid,'c+',label='A22plot L')
        #plt.plot(eosdat_list[A22mcp['id']].P_solid,eosdat_list[A22mcp['id']].T_solid,'cx',label='A22plot S')

        plt.plot(eosdat_list[U10mcl['id']].P_liquid,eosdat_list[U10mcl['id']].T_liquid,'k+',label=U10mcl['lab'])
        plt.plot(eosdat_list[U10mcs['id']].P_solid,eosdat_list[U10mcs['id']].T_solid,'kx',label=U10mcs['lab'])
        plt.plot(eosdat_list[G24mcl['id']].P_liquid,eosdat_list[G24mcl['id']].T_liquid,'r+',label=G24mcl['lab'])
        plt.plot(eosdat_list[G24mcs['id']].P_solid,eosdat_list[G24mcs['id']].T_solid,'rx',label=G24mcs['lab'])
        #plt.plot(eosdat_list[3].P[7::],eosdat_list[3].T[7::],'--',label='U10 CaCl2-PbO2')

        # stv melt curve Shen and Lazer 1995
        #plt.errorbar(eosdat_list[S95mcstv['id']].P,eosdat_list[S95mcstv['id']].T,\
        #             xerr=eosdat_list[S95mcstv['id']].P_std,yerr=eosdat_list[S95mcstv['id']].T_std,\
        #             fmt='d',color='orange',label=S95mcstv['lab'])
        # Zhang et al. JGR 1993 triple point Cs-Stv-L
        #plt.plot(13.7, 2800+273,'rs',label='Z93 Cs-Stv-L TP exp')
        # add literature points
        # plt.errorbar([500],[8300],xerr=[30],yerr=[300],fmt='bX',label='M15 T_shock plateau')

    if fig1 or fig10 or fig13 or fig14 or fig15: # or fig26 or fig27:
            #plot fit to Geng and Mohn 2024
            pg24arr = np.arange(5000)+50 # GPa
            tg24arr = 6.018*np.power((1+(pg24arr+202.338)/4.736e-04),0.5126)
            plt.plot(pg24arr,tg24arr,'--',color='orange',label='G24 2-phase PBE')
            if 0:
                for ii in range(len(pg24arr)):
                    print(pg24arr[ii],', ',tg24arr[ii])
            pg24arr = np.arange(5000)+50 # GPa
            tg24arr = 6.426*np.power((1+(pg24arr+256.011)/6.474e-04),0.5141)
            #plt.plot(pg24arr,tg24arr,'-',color='darkorange',label='G24 Z method GGA')

            ind = np.asarray([1,2,3,9,10,11,12,13])
            plt.plot(eosdat_list[P25mcpbe['id']].P[ind],eosdat_list[P25mcpbe['id']].T[ind],':',color='m',label='P26 2-phase PBEsol')
            ind = np.asarray([0,1,2,3,4,10,11,12])
            plt.plot(eosdat_list[P25mcscan['id']].P[ind],eosdat_list[P25mcscan['id']].T[ind],'-.',color='m',label='P26 2-phase SCAN')
            plt.plot(eosdat_list[B15mcpt['id']].P,eosdat_list[B15mcpt['id']].T,'--',color='dimgrey',label='B15 Gibbs model')

    if fig12 or fig24 or fig25:# or fig28:
        #plt.fill_between(omegastvp,omegastvt-omegastvt_err,omegastvt+omegastvt_err,label='M15 stv hug',alpha=0.5)
        p_millot = np.arange(6000)
        t_millot = 1968.5+307.8*np.power(p_millot,0.485)
        plt.plot(p_millot,t_millot,'--',color='dimgrey',linewidth=2,label='L83-M15 melt curve')
        #ind = np.where((p_millot > 450) & (p_millot < 520))
        #plt.plot(p_millot[ind],t_millot[ind],'-',color='black')
    
    if fig28:
        pg24arr = np.arange(440)+60 # GPa
        tg24arr = 6.426*np.power((1+(pg24arr+256.011)/6.474e-04),0.5141)
        #plt.plot(pg24arr,tg24arr,'--',color='gold',label='G24 melt curve (Z method GGA)',linewidth=2)
        pg24arr = np.arange(700)+50 # GPa
        tg24arr = 6.018*np.power((1+(pg24arr+202.338)/4.736e-04),0.5126)
        #plt.plot(pg24arr,tg24arr,'-',color='c',label='G24 melt curve 2phase PBE')

        #stv_dat = np.loadtxt('../data_silica_db/Millot-2015-Hug.csv',skiprows=1,delimiter=',')
        plt.plot(stv_hug_M15rev.P,stv_hug_M15rev.T,color='darkgoldenrod')
        ind = np.where(stv_hug_M15rev.P > stv_wiggle2.P*.995)
        plt.fill_between(stv_hug_M15rev.P[ind],stv_hug_M15rev.T[ind]-stv_hug_M15rev.T_std[ind],stv_hug_M15rev.T[ind]+stv_hug_M15rev.T_std[ind],alpha=0.4,color='darkorange',label='Rev. M15 stv decaying: Opaque front (vis)')
        ind = np.where(stv_hug_M15rev.P < stv_wiggle2.P*1.005)
        plt.fill_between(stv_hug_M15rev.P[ind],stv_hug_M15rev.T[ind]-stv_hug_M15rev.T_std[ind],stv_hug_M15rev.T[ind]+stv_hug_M15rev.T_std[ind],alpha=0.25,color='grey',label='Rev. M15 stv decaying: no P cal (vis)')

        
        
        #cst_dat = np.loadtxt('../data_silica_db/cst-plot-data.csv',skiprows=1,delimiter=',')
        #plt.plot(cst_dat[:,0],cst_dat[:,1],color='salmon')
        #plt.fill_between(cst_dat[:,0],cst_dat[:,1]-400,cst_dat[:,1]+400,alpha=0.25,color='salmon',label='F24 cst Hugoniot')

    if fig28:
        # columns 0 for P and 1 for Us and 9 for T and 10 for T_err
        #C84939 = np.loadtxt('../data_silica_db/cst-43939.csv',skiprows=1,delimiter=',')
        #C84941 = np.loadtxt('../data_silica_db/cst-43941.csv',skiprows=1,delimiter=',')
        #C84943 = np.loadtxt('../data_silica_db/cst-43943.csv',skiprows=1,delimiter=',')
        #plt.plot(C84939[:,0],C84939[:,9],'.',alpha=0.3)
        #plt.plot(C84941[:,0],C84941[:,9],'.',alpha=0.3)
        #plt.plot(C84943[:,0],C84943[:,9],'.',alpha=0.3)
        #plt.plot(cst_hug_liq.P,cst_hug_liq.T,'--',color='salmon',linewidth=4,label='New Cst Hug')
        if 0:
            # liquid
            ind = np.where(C43939.up>cst_hug_cross.up)
            plt.plot(C43939.P[ind],C43939.T[ind],'.',color='salmon',alpha=0.4)
            plt.plot(C43939.P[ind][0],C43939.T[ind][0],'.',color='salmon',alpha=0.4,markersize=12,label='Cst decay: liquid region')
            plt.plot(C43941.P,C43941.T,'.',color='salmon',alpha=0.4)
            ind = np.where(C43943.up>cst_hug_cross.up)
            plt.plot(C43943.P[ind],C43943.T[ind],'.',color='salmon',alpha=0.4)
            # opaque front
            colof = 'orchid'
            ind = np.where((C43943.up>cst_sol_cross.up) & (C43943.up<cst_hug_cross.up))
            plt.plot(C43939.P[ind],C43939.T[ind],'.',color=colof,alpha=0.4)
            plt.plot(C43939.P[ind][0],C43939.T[ind][0],'.',color=colof,alpha=0.4,markersize=12,label='Cst decay: opaque front')
            ind = np.where((C43943.up>cst_sol_cross.up) & (C43943.up<cst_hug_cross.up))
            plt.plot(C43943.P[ind],C43943.T[ind],'.',color=colof,alpha=0.3)
            if 0:
                # transition region        
                ind = np.where((C43939.up>uptrans1) & (C43939.up<cst_sol_cross.up))[0]
                plt.plot(C43939.P[ind],C43939.T[ind],'.',color='grey',alpha=0.3)#,label='Decay: temp. inversion') 
                plt.plot(C43939.P[ind][0],C43939.T[ind][0],'.',color='grey',alpha=0.4,markersize=12,label='Cst decay: temp. inversion') 
                ind = np.where((C43943.up>uptrans1) & (C43943.up<cst_sol_cross.up))[0]
                plt.plot(C43943.P[ind],C43943.T[ind],'.',color='grey',alpha=0.3)
                ind = np.where(C43939.up<uptrans1)
                plt.plot(C43939.P[ind],C43939.T[ind],'.',color='darkgoldenrod',alpha=0.3)#,label='Decay: transparent front')
                plt.plot(C43939.P[ind][0],C43939.T[ind][0],'.',color='darkgoldenrod',alpha=0.4,markersize=12,label='Cst decay: transparent front')

        
#        plt.fill_between(C43939.P,C43939.T-C43939.T_std,C43939.T+C43939.T_std,color='salmon',alpha=0.4,label='F24 cst decay: Opaque front (vis)')
        ind = np.where(C43939.us>cst_wiggle_us2)[0]
        plt.fill_between(C43939.P[ind],C43939.T[ind]-C43939.T_std[ind],C43939.T[ind]+C43939.T_std[ind],color='salmon',alpha=0.4,label='Rev. F24 cst decaying: Opaque front (vis)')
        plt.fill_between(C43941.P,C43941.T-C43941.T_std,C43941.T+C43941.T_std,color='salmon',alpha=0.4)
        plt.fill_between(C43943.P,C43943.T-C43943.T_std,C43943.T+C43943.T_std,color='salmon',alpha=0.4)
        ind = np.where(C43939.us<cst_wiggle_us2)[0]
        plt.fill_between(C43939.P[ind],C43939.T[ind]-C43939.T_std[ind],C43939.T[ind]+C43939.T_std[ind],color='darkgoldenrod',alpha=0.4,label='Rev. F24 cst decaying: Transparent front (vis)')
        #cst_pmc = [113., 190]
        #plt.plot(cst_pmc,np.interp(cst_pmc,stsmcsmooth.P,stsmcsmooth.T),'k^')

        # plot new fitted Coesite Hugoniot
#        ind = np.where(cst_hug_stv.P < cst_trans_cacl2_P)[0]
        col='black'
#        plt.plot(cst_hug_stv.P[ind],cst_hug_stv.T[ind],'-',color=col,linewidth=3,label='New cst Hugoniot')
        ind = np.where((cst_hug_cacl2.P > 77) & (cst_hug_cacl2.P < cst_solidus_P))[0]
        plt.plot(cst_hug_cacl2.P[ind],cst_hug_cacl2.T[ind],'-',color=col,linewidth=3,label='New cst Hugoniot')
        ind = np.where(cst_hug_liq.P > cst_liquidus_P)[0]
        plt.plot(cst_hug_liq.P[ind],cst_hug_liq.T[ind],'-',color=col,linewidth=3)

        ind = np.where((cst_hug_liq.P > 80) & (cst_hug_liq.P < cst_liquidus_P))[0]
        plt.plot(cst_hug_liq.P[ind],cst_hug_liq.T[ind],'--',color='hotpink',linewidth=3,label='Transition front')
        plt.text(150,3400,'opaque\nfront',rotation=47,size=12)
        plt.text(80,1600,'transparent\nfront',rotation=47,size=12)
        ind = np.where((stsmcsmooth.P > cst_solidus_P) & (stsmcsmooth.P <= cst_liquidus_P))[0]
        plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'-',color='hotpink',linewidth=3,label='Behind transition front')


        ind = np.where((stv_hug_liq.P > stv_wiggle2.P))[0]
        plt.plot(stv_hug_liq.P [ind],stv_hug_liq.T[ind],'--',color=col,linewidth=3,label='New stv Hugoniot')

    
    if fig29:
        # plot new qtz Hugoniot
        col='green'
#        ind = np.where((qtz_hug_stv.P >50) & (qtz_hug_stv.P < qtz_trans_cacl2_P))[0]
#        plt.plot(qtz_hug_stv.P[ind],qtz_hug_stv.T[ind],'-',color=col,linewidth=3,label='New qtz (z) Hugoniot')
#        ind = np.where((qtz_hug_cacl2.P > qtz_trans_cacl2_P) & (qtz_hug_cacl2.P < qtz_wiggle_P1))[0]
#        plt.plot(qtz_hug_cacl2.P[ind],qtz_hug_cacl2.T[ind],'-',color=col,linewidth=3)
        ind = np.where((qtz_hug_solid.P > 60) & (qtz_hug_solid.P < qtz_wiggle_P1*.97))[0]
        plt.plot(qtz_hug_solid.P[ind],qtz_hug_solid.T[ind],'-',color=col,linewidth=3)
        ind = np.where((stsmcsmooth.P > qtz_wiggle_P1) & (stsmcsmooth.P <= qtz_startliq_P))[0]
        plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'--',color='springgreen',linewidth=3)#,label='Behind transition front')
        ind = np.where(qtz_hug_liq.P > qtz_startliq_P)[0]
        plt.plot(qtz_hug_liq.P[ind],qtz_hug_liq.T[ind],'-',color=col,linewidth=3)

        # plot new FS Hugoniot
        ind = np.where(fs_hug_stv.P<fs_startsol_P*.95)
        plt.plot(fs_hug_stv.P[ind],fs_hug_stv.T[ind],'b-',linewidth=3,label='New FS Hugoniot')
        ind = np.where((stsmcsmooth.P > fs_wiggle_P2) & (stsmcsmooth.P < fs_startliq_P))[0]
        plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'--',color='deepskyblue',linewidth=3)#,label='Behind transition front')
        #plt.plot(fs_hug_mix1.P,fs_hug_mix1.T,'b-',linewidth=3)
        ind = np.where(fs_hug_liq.P>fs_startliq_P)[0]
        plt.plot(fs_hug_liq.P[ind],fs_hug_liq.T[ind],'b-',linewidth=3)

        col='salmon'
        ind = np.where((cst_hug_cacl2.P > 80) & (cst_hug_cacl2.P < cst_solidus_P))[0]
        plt.plot(cst_hug_cacl2.P[ind],cst_hug_cacl2.T[ind],'-',color=col,linewidth=3,label='New cst Hugoniot')
        ind = np.where(cst_hug_liq.P > cst_liquidus_P+2)[0]
        plt.plot(cst_hug_liq.P[ind],cst_hug_liq.T[ind],'-',color=col,linewidth=3)
        ind = np.where((stsmcsmooth.P > cst_solidus_P) & (stsmcsmooth.P <= cst_liquidus_P))[0]
        plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'--',color=col,linewidth=3)#,label='Behind transition front')

    if fig30:
        # plot new qtz Hugoniot
        col='green'
#        ind = np.where((qtz_hug_stv.P >50) & (qtz_hug_stv.P < qtz_trans_cacl2_P))[0]
#        plt.plot(qtz_hug_stv.P[ind],qtz_hug_stv.T[ind],'-',color=col,linewidth=3,label='New qtz (z) Hugoniot')
#        ind = np.where((qtz_hug_cacl2.P > qtz_trans_cacl2_P) & (qtz_hug_cacl2.P < qtz_wiggle_P1*.97))[0]
#        plt.plot(qtz_hug_cacl2.P[ind],qtz_hug_cacl2.T[ind],'-',color=col,linewidth=3)
        ind = np.where((qtz_hug_solid.P > 60) & (qtz_hug_solid.P < qtz_wiggle_P1*.97))[0]
        plt.plot(qtz_hug_solid.P[ind],qtz_hug_solid.T[ind],'-',color=col,linewidth=3,label='New qtz (z) Hugoniot')

        ind = np.where((stsmcsmooth.P > qtz_wiggle_P1) & (stsmcsmooth.P <= qtz_startliq_P))[0]
        #plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'--',color='springgreen',linewidth=3)#,label='Behind transition front')
        ind = np.where(qtz_hug_liq.P > qtz_startliq_P)[0]
        plt.plot(qtz_hug_liq.P[ind],qtz_hug_liq.T[ind],'-',color=col,linewidth=3)

        # plot new FS Hugoniot
        ind = np.where(fs_hug_stv.P<fs_startsol_P*.93)
        plt.plot(fs_hug_stv.P[ind],fs_hug_stv.T[ind],'b-',linewidth=3,label='New FS Hugoniot')
        ind = np.where((stsmcsmooth.P > fs_wiggle_P2) & (stsmcsmooth.P < fs_startliq_P))[0]
        #plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'--',color='deepskyblue',linewidth=3)#,label='Behind transition front')
        #plt.plot(fs_hug_mix1.P,fs_hug_mix1.T,'b-',linewidth=3)
        ind = np.where(fs_hug_liq.P>fs_startliq_P)[0]
        plt.plot(fs_hug_liq.P[ind],fs_hug_liq.T[ind],'b-',linewidth=3)

        col='salmon'
        ind = np.where((cst_hug_cacl2.P > 80) & (cst_hug_cacl2.P < cst_solidus_P))[0]
        plt.plot(cst_hug_cacl2.P[ind],cst_hug_cacl2.T[ind],'-',color=col,linewidth=3,label='New cst Hugoniot')
        ind = np.where(cst_hug_liq.P > cst_liquidus_P)[0]
        plt.plot(cst_hug_liq.P[ind],cst_hug_liq.T[ind],'-',color=col,linewidth=3)
        ind = np.where((stsmcsmooth.P > cst_solidus_P) & (stsmcsmooth.P <= cst_liquidus_P))[0]
        #plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'--',color=col,linewidth=3)#,label='Behind transition front')

        #col='firebrick'
        #ind = np.where(stv_hug_liq.P > stv_liquidus_P)[0]
        #plt.plot(stv_hug_liq.P[ind],stv_hug_liq.T[ind],'-',color=col,linewidth=3)
        #ind = np.where((stsmcsmooth.P > cst_solidus_P) & (stsmcsmooth.P <= cst_liquidus_P))[0]
        #plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'--',color=col,linewidth=3)#,label='Behind transition front')

    
    if False:
        # shock temperature data
        omegapmin=110
        zzz = np.where(M70811.P > omegapmin)[0]
        plt.fill_between(M70811.P[zzz],M70811.T[zzz]-M70811.T_std[zzz],M70811.T[zzz]+M70811.T_std[zzz],label='70811',alpha=0.1,color='green')
        zzz = np.where(M70813.P > omegapmin)[0]
        plt.fill_between(M70813.P[zzz],M70813.T[zzz]-M70813.T_std[zzz],M70813.T[zzz]+M70813.T_std[zzz],label='70813',alpha=0.1,color='green')
        zzz = np.where(M70817.P > omegapmin)[0]
        plt.fill_between(M70817.P[zzz],M70817.T[zzz]-M70817.T_std[zzz],M70817.T[zzz]+M70817.T_std[zzz],label='70817',alpha=0.1,color='green')
        zzz = np.where(M70820.P > omegapmin)[0]
        plt.fill_between(M70820.P[zzz],M70820.T[zzz]-M70820.T_std[zzz],M70820.T[zzz]+M70820.T_std[zzz],label='70820',alpha=0.1,color='green')
    if (pind==16):
        plt.errorbar(shockdat_list[H06fs['id']].P,shockdat_list[H06fs['id']].T,fmt='b+:',\
                 yerr=shockdat_list[H06fs['id']].T_std,label='H06 (decaying)')
        plt.plot(fs_hug_stv.P,fs_hug_stv.T,'b-',label='FS Hug NEW')
    if fig18:
        plt.plot(fs_hug_stv.P,fs_hug_stv.T,'b-',linewidth=3,label='FS Hugoniot NEW')
        #plt.plot(fs_hug_mix1.P,fs_hug_mix1.T,'b-',linewidth=3)
        ind = np.where(fs_hug_liq.P>fs_startliq_P)[0]
        plt.plot(fs_hug_liq.P[ind],fs_hug_liq.T[ind],'b-',linewidth=3)

    if fig3 or (pind==16):
        if not(fig24):
            omegapmin=85
            omegapmax=250
            mss = 2
            aaa = 0.2
            col9='grey'
            zzz = np.where((M70811.P > omegapmin) & (M70811.P < omegapmax))[0]
            plt.plot(M70811.P[zzz],M70811.T[zzz],'.',color=col9,markersize=mss,label='Omega qtz (z)',alpha=aaa)
            zzz = np.where((M70813.P > omegapmin) & (M70813.P < omegapmax))[0]
            plt.plot(M70813.P[zzz],M70813.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)
            zzz = np.where((M70817.P > omegapmin) & (M70817.P < omegapmax))[0]
            plt.plot(M70817.P[zzz],M70817.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)
            zzz = np.where((M70820.P > omegapmin) & (M70820.P < omegapmax))[0]
            plt.plot(M70820.P[zzz],M70820.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)
    if fig24:
        if 0:
            omegapmin=qtz_wiggle_P2
            omegapmax=250
            mss = 4
            aaa = 0.3
            col9='green'
            zzz = np.where((M70811.P > omegapmin) & (M70811.P < omegapmax))[0]
            plt.plot(M70811.P[zzz],M70811.T[zzz],'.',color=col9,markersize=mss,label='Omega qtz (z, decay)',alpha=aaa)
            zzz = np.where((M70813.P > omegapmin) & (M70813.P < omegapmax))[0]
            plt.plot(M70813.P[zzz],M70813.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)
            zzz = np.where((M70817.P > omegapmin) & (M70817.P < omegapmax))[0]
            plt.plot(M70817.P[zzz],M70817.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)
            zzz = np.where((M70820.P > omegapmin) & (M70820.P < omegapmax))[0]
            plt.plot(M70820.P[zzz],M70820.T[zzz],'.',color=col9,markersize=mss,alpha=aaa)

        qtz_us_fitmax = 35
        us_fitstart=qtz_wiggle_us2
        us_fitstart=8.5
        ind11 = np.where((M70811.us > us_fitstart) & (M70811.us < qtz_us_fitmax) )[0]
        ind13 = np.where((M70813.us > us_fitstart) & (M70813.us < qtz_us_fitmax) )[0]
        ind17 = np.where((M70817.us > us_fitstart) & (M70817.us < qtz_us_fitmax) )[0]
        ind20 = np.where((M70820.us > us_fitstart) & (M70820.us < qtz_us_fitmax) )[0]

        plt.fill_between(M70811.P[ind11],M70811.T[ind11]-M70811.T_std[ind11],M70811.T[ind11]+M70811.T_std[ind11],color='yellowgreen',alpha=0.2,label='Omega (decaying)')
        plt.fill_between(M70813.P[ind13],M70813.T[ind13]-M70813.T_std[ind13],M70813.T[ind13]+M70813.T_std[ind13],color='yellowgreen',alpha=0.2)
        plt.fill_between(M70817.P[ind17],M70817.T[ind17]-M70817.T_std[ind17],M70817.T[ind17]+M70817.T_std[ind17],color='yellowgreen',alpha=0.2)
        plt.fill_between(M70820.P[ind20],M70820.T[ind20]-M70820.T_std[ind20],M70820.T[ind20]+M70820.T_std[ind20],color='yellowgreen',alpha=0.2)

    if 0:
        mss = 2
        col9='grey'
        usmin = 9.23
        usmax = 9.5
        zzz = np.where((M70811.us > usmin) & (M70811.us < usmax))[0]
        ptrans = np.interp(M70811.T[zzz],qtz_hug_mix23.T,qtz_hug_mix23.P)
        plt.plot(ptrans,M70811.T[zzz],'.',color=col9,markersize=mss,label='Omega qtz (z) CaCl2-sft bnd')
        zzz = np.where((M70813.us > usmin) & (M70813.us < usmax))[0]
        ptrans = np.interp(M70813.T[zzz],qtz_hug_mix23.T,qtz_hug_mix23.P)
        plt.plot(ptrans,M70813.T[zzz],'.',color=col9,markersize=mss)
        zzz = np.where((M70817.us > usmin) & (M70817.us < usmax))[0]
        ptrans = np.interp(M70817.T[zzz],qtz_hug_mix23.T,qtz_hug_mix23.P)
        plt.plot(ptrans,M70817.T[zzz],'.',color=col9,markersize=mss)
        zzz = np.where((M70820.us > usmin) & (M70820.us < usmax))[0]
        ptrans = np.interp(M70820.T[zzz],qtz_hug_mix23.T,qtz_hug_mix23.P)
        plt.plot(ptrans,M70820.T[zzz],'.',color=col9,markersize=mss)
        omegapmin=np.max(M70820.P[zzz])
        omegapmax=250
        mss = 2
        col9='grey'
        zzz = np.where((M70811.P > omegapmin) & (M70811.P < omegapmax))[0]
        plt.plot(M70811.P[zzz],M70811.T[zzz],'.',color=col9,markersize=mss,label='Omega qtz (z)')
        zzz = np.where((M70813.P > omegapmin) & (M70813.P < omegapmax))[0]
        plt.plot(M70813.P[zzz],M70813.T[zzz],'.',color=col9,markersize=mss)
        zzz = np.where((M70817.P > omegapmin) & (M70817.P < omegapmax))[0]
        plt.plot(M70817.P[zzz],M70817.T[zzz],'.',color=col9,markersize=mss)
        zzz = np.where((M70820.P > omegapmin) & (M70820.P < omegapmax))[0]
        plt.plot(M70820.P[zzz],M70820.T[zzz],'.',color=col9,markersize=mss)
    if fig4 or fig5 or fig9:
        if fig5:
    #    if 0:
            omegapmin=85
            mss = 2
            zzz = np.where(M70811.P > omegapmin)[0]
            plt.plot(M70811.P[zzz],M70811.T[zzz],'g.',mfc='yellowgreen',markersize=mss,label='Omega qtz (z)')
            zzz = np.where(M70813.P > omegapmin)[0]
            plt.plot(M70813.P[zzz],M70813.T[zzz],'g.',mfc='white',markersize=mss)
            zzz = np.where(M70817.P > omegapmin)[0]
            plt.plot(M70817.P[zzz],M70817.T[zzz],'g.',mfc='yellow',markersize=mss)
            zzz = np.where(M70820.P > omegapmin)[0]
            plt.plot(M70820.P[zzz],M70820.T[zzz],'g.',mfc='tan',markersize=mss)

        if fig4:
            omegapmin=0
            zzz = np.where(M70811.P > omegapmin)[0]
            plt.plot(M70811.P[zzz],M70811.T[zzz],'g.',mfc='yellowgreen',label='Omega qtz')
            zzz = np.where(M70813.P > omegapmin)[0]
            plt.plot(M70813.P[zzz],M70813.T[zzz],'g.',mfc='white')
            zzz = np.where(M70817.P > omegapmin)[0]
            plt.plot(M70817.P[zzz],M70817.T[zzz],'g.',mfc='yellow')
            zzz = np.where(M70820.P > omegapmin)[0]
            plt.plot(M70820.P[zzz],M70820.T[zzz],'g.',mfc='tan')
    #if fig6:
    if 0:
            omegapmin=85
            mss = 2
            zzz = np.where(M70811.P > omegapmin)[0]
            plt.plot(M70811.P[zzz],M70811.T[zzz],'g.',mfc='yellowgreen',markersize=mss,label='Omega qtz')
            zzz = np.where(M70813.P > omegapmin)[0]
            plt.plot(M70813.P[zzz],M70813.T[zzz],'g.',mfc='white',markersize=mss)
            zzz = np.where(M70817.P > omegapmin)[0]
            plt.plot(M70817.P[zzz],M70817.T[zzz],'g.',mfc='yellow',markersize=mss)
            zzz = np.where(M70820.P > omegapmin)[0]
            plt.plot(M70820.P[zzz],M70820.T[zzz],'g.',mfc='tan',markersize=mss)

    if fig17 or fig26:
#        plotomegaPT(M70811,label=1,branch=6)
#        plotomegaPT(M70811,label=1,branch=7)
        plotomegaPT(M70811,label=1,branch=6)
        plotomegaPT(M70811,label=1,branch=7)
        plotomegaPT(M70813,branch=7)        
        plotomegaPT(M70813,branch=6)
        plotomegaPT(M70817,branch=7)        
        plotomegaPT(M70817,branch=6)
        plotomegaPT(M70820,branch=7)        
        plotomegaPT(M70820,branch=6)


    #############################################################################
    # NEW QUARTZ HUGONIOT
    # New quartz Hugoniot with no melt curve section
#    if fig24 or fig3 or fig17:
    if (fig3 or fig17) and not(fig24):
        ind = np.where(qtz_hug_stv.P < qtz_trans_cacl2_P)[0]
        col='black'
        plt.plot(qtz_hug_stv.P[ind],qtz_hug_stv.T[ind],'-',color=col,linewidth=4,label='New qtz (z) Hugoniot')
        ind = np.where((qtz_hug_cacl2.P > qtz_trans_cacl2_P) & (qtz_hug_cacl2.P < qtz_wiggle_P1))[0]
        plt.plot(qtz_hug_cacl2.P[ind],qtz_hug_cacl2.T[ind],'-',color=col,linewidth=4)
        ind = np.where(qtz_hug_liq.P > qtz_wiggle_P2)[0]
        plt.plot(qtz_hug_liq.P[ind],qtz_hug_liq.T[ind],'-',color=col,linewidth=4)

    # New quartz Hugoniot with transparent extension segment
    if fig26:
        ind = np.where(qtz_hug_stv.P < qtz_trans_cacl2_P)[0]
        col='black'
#        plt.plot(qtz_hug_stv.P[ind],qtz_hug_stv.T[ind],'-',color=col,linewidth=4,label='New qtz (z) Hugoniot')
#        ind = np.where((qtz_hug_cacl2.P > qtz_trans_cacl2_P) & (qtz_hug_cacl2.P < qtz_wiggle_P1))[0]
#        plt.plot(qtz_hug_cacl2.P[ind],qtz_hug_cacl2.T[ind],'-',color=col,linewidth=4)
        ind = np.where((qtz_hug_solid.P > 60) & (qtz_hug_solid.P < qtz_wiggle_P1))[0]
        plt.plot(qtz_hug_solid.P[ind],qtz_hug_solid.T[ind],'-',color=col,linewidth=4,label='New Qtz (z-cut) Hugoniot')
        ind = np.where(qtz_hug_liq.P > qtz_startliq_P)[0]
        plt.plot(qtz_hug_liq.P[ind],qtz_hug_liq.T[ind],'-',color=col,linewidth=4)
        plt.text(118,3800,'opaque\nfront',rotation=30,size=12)
        plt.text(80,2700,'transparent\nfront',rotation=30,size=12)
        ind = np.where((stsmcsmooth.P > 110) & (stsmcsmooth.P <= qtz_startliq_P))[0]
        plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'-',color='springgreen',linewidth=3,label='Behind transition front')
    
    # OLD SEGMENT
    if fig17 and (not fig22):
        ind = np.where((stsmcsmooth.P > Pt3) & (stsmcsmooth.P <= qtz_startliq_P))[0]
        plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'-',color='springgreen',linewidth=3,label='Shock states behind transition front')

    # NEW FUSED SILICA HUGONIOT
#    if fig25:
    if 0:
        ind = np.where(fs_hug_liq.P > fs_wiggle_P2)
        plt.plot(fs_hug_liq.P[ind],fs_hug_liq.T[ind],'k-',linewidth=3,label='New FS Hugoniot')
        #plt.plot(fs_hug_mix1.P,fs_hug_mix1.T,'b--',linewidth=2,label='FS Hug Mix NEW')
        ind = np.where(fs_hug_stv.P < fs_startsol_P)
        plt.plot(fs_hug_stv.P[ind],fs_hug_stv.T[ind],'k-',linewidth=3)#label='FS Hug Stv NEW')
    if fig27:
        ind = np.where(fs_hug_liq.P >  fs_startliq_P)
        plt.plot(fs_hug_liq.P[ind],fs_hug_liq.T[ind],'k-',linewidth=3,label='New FS Hugoniot')
        ind = np.where((fs_hug_liq.P > 45) & (fs_hug_liq.P < fs_startliq_P))
        plt.plot(fs_hug_liq.P[ind],fs_hug_liq.T[ind],'--',linewidth=3,color='deepskyblue',label='Transition front')
        #plt.plot(fs_hug_mix1.P,fs_hug_mix1.T,'b--',linewidth=2,label='FS Hug Mix NEW')
        ind = np.where(fs_hug_stv.P < fs_startsol_P)
        plt.plot(fs_hug_stv.P[ind],fs_hug_stv.T[ind],'k-',linewidth=3)#label='FS Hug Stv NEW')
        plt.text(78,3900,'opaque\nfront',rotation=44,size=12)
        plt.text(50,2300,'transparent\nfront',rotation=44,size=12)
        ind = np.where((stsmcsmooth.P > 70) & (stsmcsmooth.P <= fs_startliq_P))[0]
        plt.plot(stsmcsmooth.P[ind],stsmcsmooth.T[ind],'-',color='deepskyblue',linewidth=3,label='Behind transition front')

    #############################################################################
    ## PLOT DATA POINTS
    #############################################################################
    
    if fig21 or fig26:
        #ind = np.where(qtz_hug_liq.T>300)[0]
        #ind = np.where((qtz_hug_liq.T>3000) & (qtz_hug_liq.T<6000))[0]
        ind = np.where((qtz_hug_liq.P>70) & (qtz_hug_liq.P<qtz_startliq_P))[0]
        plt.plot(qtz_hug_liq.P[ind],qtz_hug_liq.T[ind],'--',color='springgreen',linewidth=2,label='Transition front')

    if fig2 or fig5 or fig6 or fig7 or fig8:
        if not(fig27):
            # Fused silica in blue
            #plt.errorbar(shockdat_list[L80fs['id']].P,shockdat_list[L80fs['id']].T,fmt='b+',xerr=shockdat_list[L80fs['id']].P_std,yerr=shockdat_list[L80fs['id']].T_std,label=L80fs['lab'])
            #plt.errorbar(shockdat_list[L83fs['id']].P,shockdat_list[L83fs['id']].T,fmt='bx',xerr=shockdat_list[L83fs['id']].P_std,yerr=shockdat_list[L83fs['id']].T_std,label=L83fs['lab'])
            #plt.errorbar(shockdat_list[B88RRfs['id']].P,shockdat_list[B88RRfs['id']].T,fmt='bo',mfc='red',yerr=shockdat_list[B88RRfs['id']].T_std,label=B88RRfs['lab'])
            if fig7:
                plt.errorbar(shockdat_list[B88fs['id']].P,shockdat_list[B88fs['id']].T,fmt='bo-',mfc='white',yerr=shockdat_list[B88fs['id']].T_std,label=B88fs['lab'])
            else:
                if not(fig25):
                    plt.errorbar(shockdat_list[B88fs['id']].P,shockdat_list[B88fs['id']].T,fmt='bo',mfc='white',yerr=shockdat_list[B88fs['id']].T_std,label=B88fs['lab'])
                ind = np.where(shockdat_list[B88RRfs['id']].us>0)
                plt.errorbar(shockdat_list[B88RRfs['id']].P[ind],shockdat_list[B88RRfs['id']].T[ind],fmt='bo',mfc='lightblue',yerr=shockdat_list[B88RRfs['id']].T_std[ind],label='L83R+B88R')#'L83R+'+B88RRfs['lab'])
                ind = np.where(shockdat_list[B88RRfs['id']].us==0)
                plt.errorbar(shockdat_list[B88RRfs['id']].P[ind],shockdat_list[B88RRfs['id']].T[ind],fmt='bo',mfc='white',yerr=shockdat_list[B88RRfs['id']].T_std[ind],label='B88 (no IM)')#B88fs['lab']+' (no IM)')
    
            plt.errorbar(shockdat_list[S80fs['id']].P,shockdat_list[S80fs['id']].T,fmt='bv',xerr=shockdat_list[S80fs['id']].P_std,yerr=shockdat_list[S80fs['id']].T_std,label='S80')#S80fs['lab'])
            plt.errorbar(shockdat_list[S81fs['id']].P,shockdat_list[S81fs['id']].T,fmt='bx',label='S81 (no IM)')#S81fs['lab']+' (no IM)')
            if not(fig5):
                plt.errorbar(shockdat_list[S82fs['id']].P,shockdat_list[S82fs['id']].T,fmt='b^',label='S82 620nm')#S82fs['lab']+'-620nm')
                ind = np.where(shockdat_list[S82fs['id']].T_b >0)[0]
                plt.errorbar(shockdat_list[S82fs['id']].P[ind],shockdat_list[S82fs['id']].T_b[ind],fmt='b^',mfc='white',label='S82 810nm')#S82fs['lab']+'-810nm')
        
            # Hicks decaying shock
            #plt.plot(shockdat_list[H06fs['id']].P,shockdat_list[H06fs['id']].T,'b+:',label='H06 FS')
            plt.errorbar(shockdat_list[H06fs['id']].P,shockdat_list[H06fs['id']].T,fmt='b+-',\
                     yerr=shockdat_list[H06fs['id']].T_std,label='H06 (decaying)')
        if not(fig25 or fig27):
            plt.plot(fs_hug_liq.P,fs_hug_liq.T,'b-',label='FS Hug Liq NEW')
            #plt.plot(fs_hug_mix1.P,fs_hug_mix1.T,'b--',linewidth=2,label='FS Hug Mix NEW')
            plt.plot(fs_hug_stv.P,fs_hug_stv.T,'b--',label='FS Hug Stv NEW')
        if not(fig7 or fig8 or fig5):
            ind = np.where(shockdat_list[S81fs['id']].P <55)
            fstpfit = np.polyfit(shockdat_list[S81fs['id']].P[ind],shockdat_list[S81fs['id']].T[ind],1)
            #plt.plot(shockdat_list[S81fs['id']].P[ind],fstpfit[1]+fstpfit[0]*shockdat_list[S81fs['id']].P[ind],'b-',label='S81 fs fit')
            #plt.plot(np.arange(50,68,1),fstpfit[1]+fstpfit[0]*np.arange(50,68,1),'b--')
            #print('Fused silica S81 P-T fit = ',fstpfit)
        #ind = np.where((shockdat_list[B88fs['id']].P >69) & (shockdat_list[B88fs['id']].P <90))
        #fstpfit = np.polyfit(shockdat_list[B88fs['id']].P[ind],shockdat_list[B88fs['id']].T[ind],1)
        #plt.plot(shockdat_list[B88fs['id']].P[ind],fstpfit[1]+fstpfit[0]*shockdat_list[B88fs['id']].P[ind],'b-')


        
        if not(fig7 or fig8 or fig5):
            plt.errorbar(shockdat_list[A23fs['id']].P,shockdat_list[A23fs['id']].T,fmt='bD',mfc='white',xerr=shockdat_list[A23fs['id']].P_std,yerr=shockdat_list[A23fs['id']].T_std,label=A23fs['lab'])
            #[5027.3122   6703.340637 4573.230527 4200.1755   4483.578691 4611.0908 -999.       4670.915    5239.       5170.      
            #irind = np.asarray([4,4,4,4,2,4,0,4,1,1]) # 4 is 4.8 microns, 2 is 2.3 microns, 1 is visible
            #itmp = np.where(irind >=2)[0]
            #plt.errorbar(shockdat_list[A23fs['id']].P[itmp],shockdat_list[A23fs['id']].T[itmp],fmt='bD',mfc='pink',xerr=shockdat_list[A23fs['id']].P_std[itmp],yerr=shockdat_list[A23fs['id']].T_std[itmp],label='This Work (IR)')
        
    if fig2:
        plt.plot(eosdat_list[G24mcl['id']].P_liquid,eosdat_list[G24mcl['id']].T_liquid,'r+',label=G24mcl['lab'])
        plt.plot(eosdat_list[G24mcs['id']].P_solid,eosdat_list[G24mcs['id']].T_solid,'rx',label=G24mcs['lab'])


    

    if fig3 or fig4 or fig5 or fig6 or fig7 or fig9:
        if not(fig24):
            # alpha quartz in green
            #plt.errorbar(shockdat_list[L80qtz['id']].P,shockdat_list[L80qtz['id']].T,fmt='g+',xerr=shockdat_list[L80qtz['id']].P_std,yerr=shockdat_list[L80qtz['id']].T_std,label=L80qtz['lab'])
            #plt.errorbar(shockdat_list[L83qtz['id']].P,shockdat_list[L83qtz['id']].T,fmt='gx',xerr=shockdat_list[L83qtz['id']].P_std,yerr=shockdat_list[L83qtz['id']].T_std,label=L83qtz['lab'])
            if fig7:
                plt.errorbar(shockdat_list[B88qtz['id']].P,shockdat_list[B88qtz['id']].T,fmt='go-',mfc='white',yerr=shockdat_list[B88qtz['id']].T_std,label=B88qtz['lab'])
            else:
                plt.errorbar(shockdat_list[B88qtz['id']].P,shockdat_list[B88qtz['id']].T,fmt='go',mfc='white',yerr=shockdat_list[B88qtz['id']].T_std,label=B88qtz['lab']+'(no IM)')
                ind = np.where(shockdat_list[B88RRqtz['id']].us>0)
                #plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='red',yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label=B88RRqtz['lab'])
                #        ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'z')[0]
                plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='yellow',\
                     yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label='L83R+'+B88RRqtz['lab']) # z cut
                #ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'pz')[0]
                #plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='orange',\
                #     yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label=B88RRqtz['lab']+' (pz)')        
    if fig3:
        plt.errorbar(shockdat_list[A23qtz['id']].P,shockdat_list[A23qtz['id']].T,fmt='gD',mfc='yellow',xerr=shockdat_list[A23qtz['id']].P_std,yerr=shockdat_list[A23qtz['id']].T_std,label=A23qtz['lab'])

    if fig24:

      #plt.errorbar(shockdat_list[B88qtz['id']].P,shockdat_list[B88qtz['id']].T,fmt='go',mfc='black',yerr=shockdat_list[B88qtz['id']].T_std,label=B88qtz['lab']+' (noim)',markersize=10)
        
        print(shockdat_list[B88qtz['id']].P)
        print(shockdat_list[B88RRqtz['id']].P)
        print(shockdat_list[B88RRqtz['id']].T)
        print(shockdat_list[B88qtz['id']].rem)
        #ind = np.where(shockdat_list[B88RRqtz['id']].us>0)
        #plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='red',yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label=B88RRqtz['lab'])
        #ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'z')[0]
        #ind = np.where((shockdat_list[B88RRqtz['id']].orientation == 'z') & (shockdat_list[B88RRqtz['id']].T > 0))[0]
        ind=[1,2,4]
        plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='yellow',\
             yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label='L83R+B88R')#'+B88RRqtz['lab']) # z cut
        ind = [1] # no revised IM z
        plt.errorbar(shockdat_list[B88qtz['id']].P[ind],shockdat_list[B88qtz['id']].T[ind],fmt='gs',mfc='white',label='B88 (no IM)')#label=B88qtz['lab']+' (no IM)')
  
        #ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'pz')[0]
        #ind = np.where((shockdat_list[B88RRqtz['id']].orientation == 'pz') & (shockdat_list[B88RRqtz['id']].T>0))[0]
        #plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='orange',\
             #yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label=B88RRqtz['lab']+' (pz)')        
        ind = [0,3]
        plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='orange',\
             yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label='L83R+B88R (pz)')#+B88RRqtz['lab']+' (pz)')              
        ind = [0,3,5] # no revised IM pz
        plt.errorbar(shockdat_list[B88qtz['id']].P[ind],shockdat_list[B88qtz['id']].T[ind],fmt='gs',mfc='bisque',label='B88 (pz, no IM)')#B88qtz['lab']+' (pz, no IM)')

        plt.errorbar(shockdat_list[L04qtz['id']].P,shockdat_list[L04qtz['id']].T,fmt='gv',mfc='orange',\
                     xerr=shockdat_list[L04qtz['id']].P_std,yerr=shockdat_list[L04qtz['id']].T_std,label='L04 (x)')#L04qtz['lab']+' (x)')
        
        plt.errorbar(shockdat_list[S81qtz['id']].P,shockdat_list[S81qtz['id']].T,fmt='gx',label='S81 (-, no IM)')#S81qtz['lab']+' (-, no cal)')
        # Luo et al. 2003 quartz Hugoniot model shock and release
        if fig11:
            plt.errorbar(shockdat_list[L03qtz['id']].P,shockdat_list[L03qtz['id']].T,fmt='r>',label=L03qtz['lab'])
        
        # Hicks decaying shock
        #plt.plot(qtz_hug_Hicks06.P,qtz_hug_Hicks06.T,'g+:',label='H06 qtz')
        #if fig7:
        #    plt.errorbar(qtz_hug_Hicks06.P,qtz_hug_Hicks06.T,fmt='g+:',\
        #                 yerr = qtz_hug_Hicks06.T_std,label='H06 qtz')
        #else:
        #plt.errorbar(qtz_hug_Hicks06.P[3::],qtz_hug_Hicks06.T[3::],fmt='g+',\
        #                 yerr = qtz_hug_Hicks06.T_std[3::],label='H06 qtz (z, decay)')
        #plt.fill_between(qtz_hug_Hicks06.P[3::],qtz_hug_Hicks06.T[3::]-qtz_hug_Hicks06.T_std[3::],qtz_hug_Hicks06.T[3::]+qtz_hug_Hicks06.T_std[3::],color='green',alpha=0.2,label='H06 qtz(z)')
        #plt.fill_between(qtz_hug_Hicks06.P,qtz_hug_Hicks06.T-qtz_hug_Hicks06.T_std,qtz_hug_Hicks06.T+qtz_hug_Hicks06.T_std,color='green',alpha=0.2,label='H06 qtz(z)')
        if not(fig7 or fig8 or fig5 or fig9):
            ind = np.where((shockdat_list[B88qtz['id']].P >75) & (shockdat_list[B88qtz['id']].P <110))
            qtztpfit = np.polyfit(shockdat_list[B88qtz['id']].P[ind],shockdat_list[B88qtz['id']].T[ind],2)
            #print(shockdat_list[B88qtz['id']].P[ind],shockdat_list[B88qtz['id']].T[ind])
            #print('quartz Boslough P-T fit 75-110 GPa: ',qtztpfit)
            #plt.plot(shockdat_list[B88qtz['id']].P[ind],qtztpfit[1]+qtztpfit[0]*shockdat_list[B88qtz['id']].P[ind],'g:',label='Likely d-NiAs')
            qtzptmp = np.arange(75,110,1)
            qtzttmp = qtztpfit[2]+qtztpfit[1]*qtzptmp+qtztpfit[0]*qtzptmp*qtzptmp
            #plt.plot(qtzptmp,qtzttmp,'g--',label='CaCl2 region P-T fit')
#        if not(fig7 or fig8 or fig5 or fig9):
           # plt.errorbar(shockdat_list[A23qtz['id']].P,shockdat_list[A23qtz['id']].T,fmt='gD',mfc='white',xerr=shockdat_list[A23qtz['id']].P_std,yerr=shockdat_list[A23qtz['id']].T_std,label=A23qtz['lab'])
       #if fig8 or fig9:
            if not(fig24):
                plt.errorbar(shockdat_list[A23qtz['id']].P,shockdat_list[A23qtz['id']].T,fmt='gD',mfc='yellow',xerr=shockdat_list[A23qtz['id']].P_std,yerr=shockdat_list[A23qtz['id']].T_std,label='This work qtz')
            #plt.errorbar([84.709,102.44],[4484,5384],yerr=[212,264],fmt='gD',label='New qtz(z)')
            #plt.errorbar([84.709,98.03],[4484,5384],yerr=[212,264],fmt='md',label='Shots 47,48 VIS Us')

    if fig24:
        ind = np.where(shockdat_list[H06qtz['id']].T>0)
        plt.errorbar(shockdat_list[H06qtz['id']].P[ind],shockdat_list[H06qtz['id']].T[ind],fmt='g+-',\
                         yerr = shockdat_list[H06qtz['id']].T_std[ind],label='H06 (decaying)')
        ind = np.where(shockdat_list[A23qtz['id']].P > qtz_singleshock_P)[0]
        plt.errorbar(shockdat_list[A23qtz['id']].P[ind],shockdat_list[A23qtz['id']].T[ind],fmt='gD',mfc='yellow',markersize=7,xerr=shockdat_list[A23qtz['id']].P_std[ind],yerr=shockdat_list[A23qtz['id']].T_std[ind],label='This work (vis)')
        ind = np.where(shockdat_list[A23qtz['id']].P < qtz_singleshock_P)[0]
        plt.errorbar(shockdat_list[A23qtz['id']].P[ind],shockdat_list[A23qtz['id']].T[ind],fmt='gD',mfc='white',markersize=7,xerr=shockdat_list[A23qtz['id']].P_std[ind],yerr=shockdat_list[A23qtz['id']].T_std[ind],label='This work (2-wave)')

    
    if fig26:
        #plt.errorbar(shockdat_list[B88qtz['id']].P,shockdat_list[B88qtz['id']].T,fmt='go',mfc='white',yerr=shockdat_list[B88qtz['id']].T_std,label=B88qtz['lab']+'(no IM)')
        indtmp = np.where((shockdat_list[B88RRqtz['id']].orientation == 'z') & (shockdat_list[B88RRqtz['id']].T > 0))[0]
        print(shockdat_list[B88RRqtz['id']].P[indtmp])
        ind=indtmp[0:3]
        plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='yellow',\
             yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label='Steady L83R+'+B88RRqtz['lab'])#+' (z)')
        
        
        #plt.errorbar(shockdat_list[A23qtz['id']].P,shockdat_list[A23qtz['id']].T,fmt='gD',mfc='yellow',markersize=6,xerr=shockdat_list[A23qtz['id']].P_std,yerr=shockdat_list[A23qtz['id']].T_std,label='Steady This work qtz')
        ind = np.where( (shockdat_list[A23qtz['id']].us < 9.5) & (shockdat_list[A23qtz['id']].us>7.7))[0]
        plt.errorbar(shockdat_list[A23qtz['id']].P[ind],shockdat_list[A23qtz['id']].T[ind],fmt='gD',mfc='yellow',label='Steady This work',\
                    xerr=shockdat_list[A23qtz['id']].P_std[ind],yerr=shockdat_list[A23qtz['id']].T_std[ind])
        ind = np.where( (shockdat_list[A23qtz['id']].us < 7.7))[0]
        plt.errorbar(shockdat_list[A23qtz['id']].P[ind],shockdat_list[A23qtz['id']].T[ind],fmt='gD',mfc='white',label='Steady This work (2-wave)',\
                    xerr=shockdat_list[A23qtz['id']].P_std[ind],yerr=shockdat_list[A23qtz['id']].T_std[ind])

        #ind = np.where((shockdat_list[S81qtz['id']].P <90))
        #NiAstpfit = np.polyfit(shockdat_list[S81qtz['id']].P[ind],shockdat_list[S81qtz['id']].T[ind],1)
        #print(shockdat_list[S81qtz['id']].P[ind],shockdat_list[S81qtz['id']].T[ind])
        #print('NiAs P-T fit < 90 GPa: ',NiAstpfit)
        #plt.plot(shockdat_list[S81qtz['id']].P[ind],NiAstpfit[1]+NiAstpfit[0]*shockdat_list[S81qtz['id']].P[ind],'r-',label='d-NiAs region fit')


        #overlaps with Omega data
        #ind = np.where(shockdat_list[B88qtz['id']].P >115)
        #qtztpfit = np.polyfit(shockdat_list[B88qtz['id']].P[ind],shockdat_list[B88qtz['id']].T[ind],1)
        #plt.plot(shockdat_list[B88qtz['id']].P[ind],qtztpfit[1]+qtztpfit[0]*shockdat_list[B88qtz['id']].P[ind],'g-')

    if fig8 and not(fig27):
        #plt.errorbar(shockdat_list[A23fs['id']].P,shockdat_list[A23fs['id']].T,fmt='bD',mfc='white',xerr=shockdat_list[A23fs['id']].P_std,yerr=shockdat_list[A23fs['id']].T_std,label=A23fs['lab'])
        plt.errorbar(shockdat_list[A23fs['id']].P,shockdat_list[A23fs['id']].T,fmt='bD',mfc='lightblue',markersize=7,xerr=shockdat_list[A23fs['id']].P_std,yerr=shockdat_list[A23fs['id']].T_std,label='This work (vis)')
        #[5027.3122   6703.340637 4573.230527 4200.1755   4483.578691 4611.0908 -999.       4670.915    5239.       5170.      
        irind = np.asarray([4,4,4,4,2,4,0,4,1,1]) # 4 is 4.8 microns, 2 is 2.3 microns, 1 is visible
        itmp = np.where(irind >=2)[0]
        plt.errorbar(shockdat_list[A23fs['id']].P[itmp],shockdat_list[A23fs['id']].T[itmp],fmt='bD',mfc='plum',xerr=shockdat_list[A23fs['id']].P_std[itmp],yerr=shockdat_list[A23fs['id']].T_std[itmp],label='This Work (IR)')
        
#    if fig26 or fig27:
    if 0:
        plt.plot([qtz_wiggle_P1,qtz_wiggle_P1],[0,10000])
        plt.plot([qtz_wiggle_P2,qtz_wiggle_P2],[0,10000])

    #        coldecay='deepskyblue'
#        col48 = 'cornflowerblue'

    if fig27:
        coldecayt='darkgoldenrod'
        with open('shot50/objs.pkl', 'rb') as f:  # Python 3: open(..., 'wb')
            s50time, s50temp, s50temperr, s50pres, s50us, s50up, ind50, isort50,s50tempbackout,Tsteady,Tsteady_err,Psteady=pickle.load(f)
        with open('shot60/objs.pkl', 'rb') as f:  # Python 3: open(..., 'wb')
            s60time, s60temp, s60temperr, s60pres, s60us, s60up, ind60, isort60=pickle.load(f)
        plt.plot(s50pres[ind50],s50temp[ind50],'-',color=coldecayt)
        plt.fill_between(s50pres[ind50[isort50]],s50temp[ind50[isort50]]-s50temperr[ind50[isort50]],s50temp[ind50[isort50]]+s50temperr[ind50[isort50]], color=coldecayt, alpha = 0.35,label='Decaying Transparent front (vis)')
        # plot the shot 50 IM solution point
        #68.23569784953546 1.4210854715202004e-14
        #5170.920498064298 19.364571422900298
        plt.errorbar(68.2,5170,fmt='bD',mfc='lightblue',yerr=19,label='Steady This work')#,label='Vis decaying shock')
        #
        plt.plot(s60pres[ind60],s60temp[ind60],'-',color=coldecay)
        plt.fill_between(s60pres[ind60[isort60]],s60temp[ind60[isort60]]-s60temperr[ind60[isort60]],s60temp[ind60[isort60]]+s60temperr[ind60[isort60]], color=coldecay, alpha = 0.25,label='Decaying Opaque front (vis)')
        # plot the shot 60 IM solution point
        plt.errorbar(87.7,5250,fmt='bD',mfc='lightblue',xerr=.3,yerr=27)#,label='Vis decaying shock')
        plt.plot(s50pres[ind50[40::]],s50tempbackout[ind50[40::]],'-',color=col48)
        plt.fill_between(s50pres[ind50[40::]],s50tempbackout[ind50[40::]]-Tsteady_err,s50tempbackout[ind50[40::]]+Tsteady_err, color=col48, alpha = 0.25,label='Decaying Opaque front (4.8 μm)')

        
        #plt.errorbar(shockdat_list[S81fs['id']].P,shockdat_list[S81fs['id']].T,fmt='bx',label='Steady '+S81fs['lab'])
        ind = np.where((shockdat_list[B88RRfs['id']].T>0) & (shockdat_list[B88RRfs['id']].us>0))[0]
        plt.errorbar(shockdat_list[B88RRfs['id']].P[ind],shockdat_list[B88RRfs['id']].T[ind],fmt='bo',mfc='lightblue',yerr=shockdat_list[B88RRfs['id']].T_std[ind],label='Steady L83R+'+B88RRfs['lab'])
        #plt.errorbar(shockdat_list[B88fs['id']].P[ind],shockdat_list[B88fs['id']].T[ind],fmt='bo',mfc='white',yerr=shockdat_list[B88fs['id']].T_std[ind],label=B88fs['lab'])
        # Hicks decaying shock
        #plt.plot(shockdat_list[H06fs['id']].P,shockdat_list[H06fs['id']].T,'b+:',label='H06 FS')
        plt.errorbar(shockdat_list[H06fs['id']].P,shockdat_list[H06fs['id']].T,fmt='b+',\
                 yerr=shockdat_list[H06fs['id']].T_std,label='Decay H06 FS')
        #ind=np.where(shockdat_list[A23fs['id']].P<55.)
        #ind=7
        #plt.errorbar(shockdat_list[A23fs['id']].P[ind],shockdat_list[A23fs['id']].T[ind],fmt='kD',mfc='blue',markersize=7,xerr=shockdat_list[A23fs['id']].P_std[ind],yerr=shockdat_list[A23fs['id']].T_std[ind],label='Gas gun FS')
        if 0:
            plt.plot([fs_wiggle_P1,fs_wiggle_P1],[0,10000])
            plt.plot([fs_wiggle_P2,fs_wiggle_P2],[0,10000])
            plt.plot([80,80],[0,10000])
            plt.plot([90,90],[0,10000])
    
    if fig3 or fig6:
        if not(fig24):
            #ind = np.where(qtz_hug_KD.P > Pt3)[0]
            #plt.plot(qtz_hug_KD.P[ind],qtz_hug_KD.T[ind],'g-',label='Qtz Z fit',linewidth=3)
            ind = np.where(qtz_hug_liq.P > qtz_startliq_P)[0]
            plt.plot(qtz_hug_liq.P[ind],qtz_hug_liq.T[ind],'g-',label='Qtz Hug Liq NEW')

    if 0: #fig30:
        colO05='teal'
        # add Deng et al. 2026 PBEsol boundaries
        # stv-CaCl2
        ptmp = np.arange(109-85)+85
        ttmp = 91.45*ptmp -3999
        #plt.plot(ptmp,ttmp,'--',color=colO05,label='D26 stv-CaCl₂ PBEsol')
        ind = np.where(eosdat_list[D26pb['id']].phase_solid =='stishovite/post-stishovite boundary (Eq.5)')
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'--',color=colO05,label='D26 stv-CaCl₂ PBEsol')
        
        # CaCl2-siefertite
        ptmp = np.arange(129-111)+110
        ttmp = 132.98*ptmp -11011
        #plt.plot(ptmp,ttmp,'-.',color=colO05,label='D26 CaCl₂-pyrite PBEsol')
        ind = np.where(eosdat_list[D26pb['id']].phase_solid =='post-stishovite/seifertite boundary (Eq.5)')
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'-.',color=colO05,label='D26 CaCl₂-sft PBEsol')
    # plot solid-solid phase boundaries original calculations
    if fig1 or fig23:
        colT04='turquoise'
        colO05='teal'
        # add Deng et al. 2026 PBEsol boundaries
        # stv-CaCl2
        ptmp = np.arange(109-85)+85
        ttmp = 91.45*ptmp -3999
        #plt.plot(ptmp,ttmp,'--',color=colO05,label='D26 stv-CaCl₂ PBEsol')
        ind = np.where(eosdat_list[D26pb['id']].phase_solid =='stishovite/post-stishovite boundary (Eq.5)')
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'--',color=colO05,label='D26 stv-CaCl₂ PBEsol')
        
        # CaCl2-siefertite
        ptmp = np.arange(129-111)+110
        ttmp = 132.98*ptmp -11011
        #plt.plot(ptmp,ttmp,'-.',color=colO05,label='D26 CaCl₂-pyrite PBEsol')
        ind = np.where(eosdat_list[D26pb['id']].phase_solid =='post-stishovite/seifertite boundary (Eq.5)')
        plt.plot(eosdat_list[D26pb['id']].P[ind],eosdat_list[D26pb['id']].T[ind],'-.',color=colO05,label='D26 CaCl₂-sft PBEsol')
                
        #
        # stv-CaCl2
        #if not(fig23):
            # Fischer et al. 2018 lab
        #plt.errorbar(eosdat_list[F18pb['id']].P,eosdat_list[F18pb['id']].T,fmt='3-',label=F18pb['lab'],color='orange')
        # Tsuchiya et al. 2004
        plt.errorbar(eosdat_list[T04pba['id']].P,eosdat_list[T04pba['id']].T,fmt='3-',label=T04pba['lab'],color=colT04)
        # Oganov et al. 2005
        #plt.errorbar(eosdat_list[O05pba['id']].P,eosdat_list[O05pba['id']].T,fmt='3-',label=O05pba['lab'],color=colO05)     
        #
        # CaCl2 -> alpha-PbO2
        plt.errorbar(eosdat_list[T04pbb['id']].P,eosdat_list[T04pbb['id']].T,fmt='2--',label=T04pbb['lab'],color=colT04)
        #plt.errorbar(eosdat_list[O05pbb['id']].P,eosdat_list[O05pbb['id']].T,fmt='2--',label=O05pbb['lab'],color=colO05)
        #
        # alpha-PbO2 -> pyrite
        plt.errorbar(eosdat_list[P25pb['id']].P,eosdat_list[P25pb['id']].T,fmt=':',color=colO05,label=P25pb['lab'])
        #plt.errorbar(eosdat_list[O05pbc['id']].P,eosdat_list[O05pbc['id']].T,fmt='4:',label=O05pbc['lab'],color=colO05)
        ind = np.where(eosdat_list[P25pbscan['id']].T<7250)[0]
        plt.errorbar(eosdat_list[P25pbscan['id']].P[ind],eosdat_list[P25pbscan['id']].T[ind],fmt='-.',color='paleturquoise',label='P26 sft-pyrite SCAN')

        # Fischer et al. 2018 lab
        #plt.errorbar(eosdat_list[F18pb['id']].P,eosdat_list[F18pb['id']].T,fmt='3-',label=F18pb['lab'],color='tan')
        
    # stv-CaCl2
    if 0:
        if fig1 and not(fig23):
            # Fischer et al. 2018
            plt.errorbar(eosdat_list[F18pb['id']].P,eosdat_list[F18pb['id']].T,fmt='m3-',label=F18pb['lab'])
        #parr = np.arange(75)+60 # GPa
        #tarr = 64.6*parr-2830
        #ind = np.where(tarr < 3000.)[0]
        #    plt.plot(parr[ind],tarr[ind] ,'m:', label='F18 stv-CaCl2 fit')
        #    plt.plot(parr[ind],tarr[ind] ,'m', label='F18 stv-CaCl2 exp')
        #plt.plot(parr,tarr ,'m:')
        # Tsuchiya et al. 2004
        plt.errorbar(eosdat_list[T04pba['id']].P,eosdat_list[T04pba['id']].T,fmt='3-',label=T04pba['lab'])
        plt.errorbar(eosdat_list[O05pba['id']].P,eosdat_list[O05pba['id']].T,fmt='3-',label=O05pba['lab'])        
    # CaCl2 -> alpha-PbO2
    if 0:
    #if fig1 or fig10 or fig23:
        plt.errorbar(eosdat_list[T04pbb['id']].P,eosdat_list[T04pbb['id']].T,fmt='2--',label=T04pbb['lab'])
        plt.errorbar(eosdat_list[O05pbb['id']].P,eosdat_list[O05pbb['id']].T,fmt='2--',label=O05pbb['lab'])
    #plotomegaPT(M70811,label=1,branch=1)
    if 0:#if fig21 or fig26:
#    if not(fig22): # EXTEND BOUNDARY
        #plt.plot([np.max(eosdat_list[O05pbb['id']].P+Pbbshift), Pt3],[4000,Tt3],':',color=pbcol,linewidth=pbw)
        #plt.plot([Pt2,130],[4900,4900],'m--',label='Metastable stv liquidus')
        if 0:
            ind = np.where(qtz_hug_liq.T>300)[0]
            ind = np.where((qtz_hug_liq.T>3000) & (qtz_hug_liq.T<6300))[0]
            plt.plot(qtz_hug_liq.P[ind],qtz_hug_liq.T[ind],'c--',label='Shock front transition layer')
        if 0:
            #plt.arrow(80,3500,5,800,head_width=5,head_length=25,shape='full')
            plt.annotate("", xy=(83, 4300), xytext=(80, 3500),
                arrowprops=dict(arrowstyle="->"))
            plt.annotate("", xy=(99, 5100), xytext=(95, 4000),
                arrowprops=dict(arrowstyle="->"))
            plt.annotate("", xy=(123, 6200), xytext=(120, 4800),
                arrowprops=dict(arrowstyle="->"))
            plt.annotate("", xy=(132, 6200), xytext=(130, 5500),
                arrowprops=dict(arrowstyle="->"))
            plt.text(90,2700,'Shock front transition layer',rotation=38,size=12)
            #plt.annotate("", xy=(145, 4250), xytext=(125, 3470),rotation=38,size=12,
            #    arrowprops=dict(arrowstyle="->"))
        ind2 = np.where(qtz_hug_glass.T <=6200)[0]
        # find the liquid densities that correspond to the P and T of the Hugoniot data
        rnew = []
        tnew = []
        pnew = []
        for iii in range(len(ind2)):
            karkipvt.isotherm_liquid(qtz_hug_glass.T[ind2[iii]])
            #rnew.append(np.interp(qtz_hug_glass.T[ind2[iii]],karkipvt.isotherm_p,karkipvt.rhoarr))
            pnew.append(np.interp(qtz_hug_glass.P[ind2[iii]],karkipvt.isotherm_p,karkipvt.rhoarr))
            tnew.append(qtz_hug_glass.T[ind2[iii]])


    if 0:
        # also new data in GM2024
        # Usui 2020 SiO RCN
        plt.text(14,4100,'Si-O\nRCN',horizontalalignment='center',     verticalalignment='center')
        plt.text(23,3900,'5.1',horizontalalignment='center',     verticalalignment='center')
        plt.text(50,5200,'5.6',horizontalalignment='center',     verticalalignment='center')
        plt.text(77,5600,'5.6',horizontalalignment='center',     verticalalignment='center')
        plt.text(108,5800,'6.3',horizontalalignment='center',     verticalalignment='center')
        plt.text(145,6000,'6.3',horizontalalignment='center',     verticalalignment='center')
        #plt.text(150,6000,'6.5',horizontalalignment='center',     verticalalignment='center') # Karki
        plt.text(62,4700,'stv')
        plt.text(100,5250,r'CaCl$_2$')
        plt.text(140,5600,r'$\alpha$-PbO$_2$')
        plt.text(212,6050,'pyrite')
        

    #############################################################################
    
    if 0:#fig26 or fig27:
        plt.errorbar(eosdat_list[B95mc['id']].P,eosdat_list[B95mc['id']].T,fmt='yX--',mfc='white',yerr=eosdat_list[B95mc['id']].T_std,label=B95mc['lab'])
 
    if fig3 or fig13:
        colarr = ['b','g','r','c','m','y','k','darkorange','limegreen','grey','brown','steelblue','darkviolet']

    if fig21:
        if 0:
            isort = np.argsort(lqtz_pz_list)
            first=1
            for i in isort:
                col='k'
                mfcz='yellow'
                if lqtz_z_list[i].P[0]>111:
                    mfcz='yellow'#'springgreen'    
                if first:
                    plt.errorbar(lqtz_z_list[i].P[0],lqtz_z_list[i].T[0],yerr=lqtz_z_list[i].T_std[0],fmt='o',markersize=7,ecolor=col,mec=col,mfc=mfcz,label='Steady shock - Lyzenga et al. 1983 (gas gun)')
                    first=0
                else:
                    plt.errorbar(lqtz_z_list[i].P[0],lqtz_z_list[i].T[0],yerr=lqtz_z_list[i].T_std[0],fmt='o',markersize=7,ecolor=col,mec=col,mfc=mfcz)
                #print(lqtz_z_list[i].ref[0],lqtz_z_list[i].rho[0],lqtz_z_list[i].rho_std[0],lqtz_z_list[i].P[0],lqtz_z_list[i].T[0])

        ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'z')[0][1:3]
        plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='yellow',\
                     label=B88RRqtz['lab']+' (z)')
        print('FIT BOSLOUGH Z CUT')
        print(shockdat_list[B88RRqtz['id']].us[ind])
        ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'pz')[0]
        plt.errorbar(shockdat_list[B88RRqtz['id']].P[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='orange',\
                     label=B88RRqtz['lab']+' (pz)')

        colarr = ['b','g','r','c','m','y','k','darkorange','limegreen','grey','brown','steelblue','darkviolet']
        
        isort = np.argsort(qtz_pgood_list)
        first=1
        mfcz2='cyan'
        for i in isort:
            col='k'
            mfcz='yellow'
            if qtz_good_list[i].P[0]>111:
                mfcz='springgreen'  
                mfcz='yellow'

            if first:
                plt.errorbar(qtz_good_list[i].P[0],qtz_good_list[i].T[0],xerr=qtz_good_list[i].P_std[0],yerr=qtz_good_list[i].T_std[0],fmt='D',markersize=7,ecolor=col,mec=col,mfc=mfcz,label='Steady shock - This Work (gas gun)')
                #plt.errorbar(qtz_good_list[i].P[0],qtz_good_list[i].T_rel[0],xerr=qtz_good_list[i].P_std[0],yerr=qtz_good_list[i].T_rel_std[0],fmt='d',markersize=7,ecolor=col,mec=col,mfc=mfcz2,label='Release temp - This Work (gas gun)')
                first=0
            else:
                plt.errorbar(qtz_good_list[i].P[0],qtz_good_list[i].T[0],xerr=qtz_good_list[i].P_std[0],yerr=qtz_good_list[i].T_std[0],fmt='D',markersize=7,ecolor=col,mec=col,mfc=mfcz)
            #plt.errorbar(qtz_good_list[i].P[0],qtz_good_list[i].T_rel[0],xerr=qtz_good_list[i].P_std[0],yerr=qtz_good_list[i].T_rel_std[0],fmt='d',markersize=7,ecolor=col,mec=col,mfc=mfcz2)
            print(qtz_good_list[i].ref[0],qtz_good_list[i].rho[0],qtz_good_list[i].rho_std[0],qtz_good_list[i].P[0],qtz_good_list[i].T[0])
        isort = np.argsort(qtz_p_list)
        for i in isort:
            print(qtz_list[i].ref[0],qtz_list[i].P_std[0])
            if (qtz_list[i].ref[0] == 'shot28'):
                plt.errorbar(qtz_list[i].P[0],qtz_list[i].T[0],xerr=qtz_list[i].P_std[0],yerr=qtz_list[i].T_std[0],fmt='D',markersize=7,ecolor=col,mec=col,mfc=mfcz)
                print('found shot 28')
                #plt.errorbar(qtz_list[i].P[0],qtz_list[i].T_rel[0],xerr=qtz_list[i].P_std[0],yerr=qtz_list[i].T_rel_std[0],fmt='d',markersize=7,ecolor=col,mec=col,mfc=mfcz)
                
    tval=4900
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 100) & (karkipvt.isotherm_p<160))[0]
    #plt.plot(karkipvt.rhoarr[ind],karkipvt.isotherm_p[ind],'m-.')
    p_start = 112
    rl_start = np.interp(p_start,karkipvt.isotherm_p[ind],karkipvt.rhoarr[ind])
    tval=6150
    ind = np.where((karkipvt.isotherm_p > 100) & (karkipvt.isotherm_p<160))[0]
    karkipvt.isotherm_liquid(tval)
    #plt.plot(karkipvt.rhoarr[ind],karkipvt.isotherm_p[ind],'m-.')
    p_end = 146
    rl_end = np.interp(p_end,karkipvt.isotherm_p[ind],karkipvt.rhoarr[ind])

    tval=4900
    ind = np.where((stvpvt.isotherm_p > 100) & (stvpvt.isotherm_p<160))[0]
    stvpvt.isotherm_stv(tval)
    #plt.plot(stvpvt.rhoarr[ind],stvpvt.isotherm_p[ind],'-.',color=colsolidus)
    rs_start = np.interp(p_start,stvpvt.isotherm_p[ind],stvpvt.rhoarr[ind])
    tval=6150
    ind = np.where((stvpvt.isotherm_p > 100) & (stvpvt.isotherm_p<160))[0]
    stvpvt.isotherm_stv(tval)
    rs_end = np.interp(p_end,stvpvt.isotherm_p[ind],stvpvt.rhoarr[ind])
    #plt.plot([p_start,p_end],[4900,6150],'m-.',label='Metastable stv-mc liquidus')
    #plt.plot([Pt2,p_start,p_end],[4900,4900,6150],'-.',color=colsolidus,label='Metastable stv-mc solidus')

    if fig24:
        plt.errorbar(120.5,5580,yerr=100,xerr=0.3,fmt='gD',mfc='tomato',markersize=7,label='This work (IR)')
   # if fig25:
   #     plt.plot(63,4800,'kD') # McQueen & Fritz FS 'phase change'

    ######## PLOT ISENTROPES ##############
    if fig6:
        #plt.plot(stvpvt.sp3000,stvpvt.st3000)
        ind = np.where((stvpvt.sp3500 > 14)&(stvpvt.sp3500 <= 100))[0]
        plt.plot(stvpvt.sp3500[ind],stvpvt.st3500[ind],'m--',label='stv S='+str(int(np.round(stvpvt.svalsarr[1]))))   
        ind = np.where((stvpvt.sp4000 > 14)&(stvpvt.sp4000 <= 100))[0]
        plt.plot(stvpvt.sp4000[ind],stvpvt.st4000[ind],'b--',label='stv S='+str(int(np.round(stvpvt.svalsarr[2]))))
        ind = np.where((stvpvt.sp4500 > 20)&(stvpvt.sp4500 <= 100))[0]
        plt.plot(stvpvt.sp4500[ind],stvpvt.st4500[ind],'--',color='yellowgreen',label='stv S='+str(int(np.round(stvpvt.svalsarr[3]))))
        ind = np.where((stvpvt.sp5000 > 30)&(stvpvt.sp5000 <= 100))[0]
        plt.plot(stvpvt.sp5000[ind],stvpvt.st5000[ind],'--',color='orange',label='stv S='+str(int(np.round(stvpvt.svalsarr[4]))))
        ind = np.where((stvpvt.sp5500 > 30)&(stvpvt.sp5500 <= 100))[0]
        plt.plot(stvpvt.sp5500[ind],stvpvt.st5500[ind],'r--',label='stv S='+str(int(np.round(stvpvt.svalsarr[5]))))
        #print(stvpvt.svalsarr,stvpvt.stempsarr)
    if fig6:
        ind = np.where(karkipvt.sp3000 <= 15)[0]
        plt.plot(karkipvt.sp3000[ind],karkipvt.st3000[ind],'m:',label='liq S='+str(int(np.round(karkipvt.svalsarr[0]))))   
        ind = np.where(karkipvt.sp3500 <= 15)[0]
        plt.plot(karkipvt.sp3500[ind],karkipvt.st3500[ind],'b:',label='liq S='+str(int(np.round(karkipvt.svalsarr[1]))))   
        ind = np.where(karkipvt.sp4000 <= 15)[0]
        plt.plot(karkipvt.sp4000[ind],karkipvt.st4000[ind],':',color='yellowgreen',label='liq S='+str(int(np.round(karkipvt.svalsarr[2]))))   
        #ind = np.where(karkipvt.sp4500 <= 15)[0]
        #plt.plot(karkipvt.sp4500[ind],karkipvt.st4500[ind],':',color='orange',label='liq S='+str(int(np.round(karkipvt.svalsarr[3]))))   
        ind = np.where(karkipvt.sp4500 > 0)[0]
        plt.plot(karkipvt.sp4500[ind],karkipvt.st4500[ind],':',color='orange',label='liq S='+str(int(np.round(karkipvt.svalsarr[3]))))   
        #ind = np.where(karkipvt.sp5000 <= 20)[0]
        #plt.plot(karkipvt.sp5000[ind],karkipvt.st5000[ind],'r:',label='liq S='+str(int(np.round(karkipvt.svalsarr[4]))))   
        ind = np.where(karkipvt.sp5000 >= 0)[0]
        plt.plot(karkipvt.sp5000[ind],karkipvt.st5000[ind],'r:',label='liq S='+str(int(np.round(karkipvt.svalsarr[4]))))   

    if 0:
        # Kraus et al. 2012 Entropy Tie Points for SiO2
        # density T
        #plt.plot([1/0.2125,1/0.2125,1/0.2153,1/0.1983],[334,4500,4500,4800],'co',label='K12 S tie pts')
        # P-T
        plt.plot([1.e-4,70,70,113],[334,4500,4500,4800],'co',label='K12 S tie pts')

    #################### TITLES AND FILENAMES #############################
    if fig1:
        #plt.title('SiO$_2$ Phase Diagram Data')
        figname = 'A1-SiO2-PTphase-diagram'
    if fig15:
        plt.errorbar(eosdat_list[G24mcl['id']].P,eosdat_list[G24mcl['id']].T,fmt='cd--',mfc='white',label=G24mcl['lab'])
       
        plt.title('SiO$_2$ Phase Diagram with New Proposed Melt Curve')
        figname = 'A15-SiO2-PTphase-diagram'
    if fig2:
        plt.title('SiO$_2$ Phase Diagram with Fused Silica Hugoniot')
        figname = 'A2-SiO2-PTphase-diagram-with-fs'
    if fig3:
        plt.title('SiO$_2$ Phase Diagram with Quartz Hugoniot')
        figname = 'A3-SiO2-PTphase-diagram-with-qtz'
    if fig4:
        plt.title('SiO$_2$ Phase Diagram with Quartz Hugoniot')
        figname = 'A4-SiO2-PTphase-diagram-with-qtz-omega'
    if fig5:    
        #plt.title('SiO$_2$ Phase Diagram with Fused Silica and Quartz Hugoniots')
        #plt.title('SiO$_2$ Phase Diagram with Fused Silica (blue) and Quartz Hugoniots (green)')
        figname = 'A5-SiO2-PTphase-diagram-with-both-noA23'
    if fig6:    
        plt.title('SiO$_2$ Phase Diagram with Hugoniots and Isentropes')
        #plt.title('SiO$_2$ Phase Diagram with Fused Silica (blue) and Quartz Hugoniots (green)')
        figname = 'A6-SiO2-PTphase-diagram-with-release'
    if fig7:    
        #plt.title('SiO$_2$ Phase Diagram with Shock Data')
        #plt.title('SiO$_2$ Phase Diagram with Fused Silica (blue) and Quartz Hugoniots (green)')
        figname = 'A7-SiO2-PTphase-diagram-with-prev-shock-data'
    if fig8:    
        #plt.title('SiO$_2$ Phase Diagram with Shock Data')
        #plt.title('SiO$_2$ Phase Diagram with Fused Silica (blue) and Quartz Hugoniots (green)')
        figname = 'A8-SiO2-PTphase-diagram-with-fs-new-interp'
    if fig9:    
        #plt.title('SiO$_2$ Phase Diagram with Shock Data')
        #plt.title('SiO$_2$ Phase Diagram with Fused Silica (blue) and Quartz Hugoniots (green)')
        figname = 'A9-SiO2-PTphase-diagram-with-qtz-new-interp'
    if fig10:    
        #plt.title('SiO$_2$ Phase Diagram with Shock Data')
        plt.title('SiO$_2$ Melt Curve data compared to Luo et al. 2002')
        figname = 'A10-SiO2-PTphase-diagram-Luo2002-comparison'
    if fig11:    
        #plt.title('SiO$_2$ Phase Diagram with Shock Data')
        plt.title('SiO$_2$ Qtz Hugoniot compared to Luo et al. 2003 Qtz Model')
        figname = 'A11-SiO2-PTphase-diagram-Luo2003-comparison'
    if fig12 and (pind==12):    
        #plt.title('SiO$_2$ Phase Diagram with Shock Data')
        figname = 'A12-SiO2-PTphase-diagram-higher-pressures'
    if fig12 and (pind==13):    
        #plt.title('SiO$_2$ Phase Diagram with Shock Data')
        figname = 'A13-SiO2-PTphase-diagram-higher-pressures-log'        
    if fig13:
        figname = 'A14-SiO2-PTphase-diagram-compare-to-static'        
    if pind==16:
        plt.title('New SiO$_2$ Hugoniot Models (8/30/2023)')
        figname = 'A16-SiO2-PTphase-diagram-new-mc'        
    if fig17:
        plt.title('Quartz Hugoniot States behind Shock Front')
        figname = 'A17-SiO2-PTphase-diagram-new-qtz-hug'        
    if fig18:
        figname = 'A18-SiO2-PTphase-diagram-new-fs-hug'              
    if fig19:
        figname = 'A19-SiO2-PTphase-check'
    if fig20:
        figname = 'book-logo.png'
    if fig21:
        plt.title('')
        figname = 'A21-SiO2-PT-metastable-shock-front'
    if fig22:
        plt.title('')
        figname = 'A22-SiO2-PT'
    if fig23:
        plt.title('')
        figname = 'A23-paper-fig1a'
    if fig24:
        plt.title('')
        figname = 'A24-paper-fig2b'
    if fig25:
        plt.title('')
        figname = 'A25-paper-fig2a'
    if fig26:
        plt.title('')
        figname = 'A26-paper-fig7b'
    if fig27:
        plt.title('')
        figname = 'A27-paper-fig7a'    
    if fig28:
        plt.title('')
        figname = 'A28-paper-fig7c'    
    if fig29:
        plt.title('')
        figname = 'A29-paper-fig9'   
    if fig30:
        plt.title('')
        figname = 'A30-paper-fig9' 
    if pind == 2425:
        plt.title('')
        figname = 'Fig2-testcombo' 
    if pind == 2627:
        plt.title('')
        figname = 'Fig7-testcombo' 
    ################## EXTRA ANNOTATIONS ###########################
    
    if not(fig12 or fig21 or fig28) and (pind != 124) and (pind !=125):
        ts = 14
        if fig22 or fig23 or fig24 or fig25 or fig26 or fig27 or fig29:
            ts=14
        plt.text(30,6600,'liquid',fontsize=ts)
    #    plt.text(20,1250,'stishovite')
        plt.text(11,950,'stishovite\nrutile-type',fontsize=ts)
        if fig23:
            plt.text(70,945,'CaCl₂-\ntype',fontsize=ts)
        else:
            plt.text(70,950,'CaCl₂-\ntype',fontsize=ts)
            
        #plt.text(70,950,'type',fontsize=ts)
        #plt.text(130,950,'seifertite',fontsize=ts)
        plt.text(130,950,'seifertite\nα-PbO₂-type',fontsize=ts)
        if not(fig8 or fig9) and not(pind==124) and not(pind==125):
            plt.text(205,990,'pyrite-\ntype',fontsize=ts)
        if fig27:
            plt.text(205,950,'pyrite-type',fontsize=ts)
        if fig23:
            plt.text(12,4800,'LDM',fontsize=ts)
            plt.text(55,6000,'HDM',fontsize=ts)
    if fig21:
        ts = 10
        yy = 750
        yy2 = 950
        #plt.text(40,6600,'liquid',fontsize=ts)
    #    plt.text(20,1250,'stishovite')
        plt.text(11,yy,'stishovite \nrutile-type',fontsize=ts)
        plt.text(67,yy2,r'CaCl$_2$-',fontsize=ts)
        plt.text(67,yy,'type',fontsize=ts)
        plt.text(107,yy2,r'seifertite',fontsize=ts)
        plt.text(107,yy,'α-PbO$_2$-type',fontsize=ts)
        if not(fig8 or fig9):
            plt.text(205,yy,'pyrite-type',fontsize=ts)

    ## add annotations
    #if fig24:
    #    plt.text(5,9500,'B. α-quartz',fontsize=14)
    if fig25 and not(fig27) and not(pind==125):
        #plt.text(5,9500,'A. Fused silica',fontsize=14)
        plt.text(205,950,'pyrite-type',fontsize=ts)
    #if fig26:
    #    plt.text(10,7500,'A. α-quartz',fontsize=16)

    ################### Axis labels and dimensions ###################################
    plt.xlim(0,250)
    plt.ylim(0,8000)

    if fig2 or fig3:
        plt.xlim(0,250)
        plt.ylim(0,10000)
    if fig7:
        plt.xlim(0,250)
        plt.ylim(0,10000)
    if fig8:
        plt.xlim(0,200)
        plt.ylim(0,8000)
    if fig12 and (pind==13):
        plt.xlim(0,5700)
        plt.ylim(0,25000)
    if fig12 and (pind ==14):
        plt.xlim(10,5700)
        plt.ylim(0,25000)
        plt.semilogx()
    if fig24 or fig25 or fig26 or fig27:
        plt.xlim(0,250)
        plt.ylim(0,10000)
    if (pind==124) or (pind==125):
        plt.xlim(0,250)
        plt.ylim(1000,9000)
    if fig28:
        plt.xlim(0,800)
        plt.ylim(0,13000)
    if fig29:
        plt.xlim(0,250)
        plt.ylim(0,10000)    
    #if fig30:
    #    plt.xlim(0,750)
    #    plt.ylim(0,15000) 
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Temperature (K)')
    
    ################### LEGEND ###################################
    ### Add items to legend
    if fig23:
        plt.text(10,7500,'A',fontsize=16)
        plt.plot([],[],' ',label=' ')#r"\textbf{Melt Curve Exp.}")
        plt.plot([],[],' ',label=' ')#'Melt Curve Calc.')
        plt.plot([],[],' ',label=' ')#'Solid boundaries')
    if fig24 or fig25:
        plt.plot([],[],' ',label=' ')#'This work')
        plt.plot([],[],' ',label=' ')#'Previous work')
        
    #get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()
    #specify order of items in legend
    order = range(len(handles))
    if not(fig20): # no legend on figure 20
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='center left', bbox_to_anchor=(1, 0.5),fontsize=10) 

    if fig23:
        #plt.plot([],[],' ',label='Solid boundaries Exp.')
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        order = range(len(handles))
        #handles.insert(0,'test')
        #order = [7,10,0,8,11,1,2,12,3,4,5,6,9,13,14,15,16,17]
        #order = [0,9,12,1,10,2,3,13,14,4,5,6,7,8,11,15,16,17,18,19,20,21]
        #order = [0,10,13,1,11,2,3,14,9,15,4,5,6,7,8,12,16,17,18,19,20,21,22]
        #order = [0,12,13,1,2,14,8,15,3,4,5,6,7,16,17,18,19,20,21,22] # before revision
        order = [0,15,16,17,2,3,9,18,4,5,6,7,8,1,19,20,21,22,10,11]
#                 12,13,1,2,14,8,15,3,4,5,6,7,16,17,18,19,20,21,22] # after revision
        print(len(order))
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12) 

    if fig24:
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        order = range(len(handles))
        #handles.insert(0,'test')
        #order = [5,13,3,4,0,6,1,2,11,7,8,9,10,12]
        #order = [5,13,3,4,0,6,1,2,11,8,9,7,10,12]
        #order = [5,7,3,4,0,6,1,2,12,9,10,8,11, 13]
        # added triple points
        #order = [ 7 , 9 , 5 , 6, 0,1, 2 , 8 , 3 , 4 ,14 ,11 ,12 ,10, 13, 15]
        #order = [ 3, 12, 2, 0, 4, 1, 10,7, 8, 6, 9,11]
        #order = [ 3, 12, 13, 2, 0, 4, 1, 10,7, 8, 6, 9,11]
        order = [13, 15, 14, 2, 0, 1, 11,6,7, 8, 9,10,12]
        print(len(order))
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right', fontsize=10,framealpha=1,title='B. α-quartz')

    if pind==124:
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right', fontsize=10,framealpha=1,title='B. α-quartz')
    if fig25:
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        #order = range(len(handles))
        #handles.insert(0,'test')
        #order = [4,12,3,0,5,1,2,7,8,9,10,6,11]
        #order = [4,13,3,0,5,1,2,8,9,10,11,6,7,12]
        # added triple points
        #order = [ 6, 15,  5, 0,1, 2,  7,  3,  4, 10, 11, 12, 13,  8,  9, 14]
        #order = [ 2, 11, 0, 3, 1, 6, 7, 8, 9, 4, 5]
        #order = [11, 0, 1, 6, 7, 8, 9, 4, 5,10] # before revision
        order = [11, 12, 0, 1, 6, 7, 8, 9, 4, 5,10] # before revision
        #order = [9, 10, 0, 1, 6, 7, 4, 5, 8]
#                 7, 8, 9, 4, 5,10] # after revision
        print(len(order))  
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='center right', fontsize=10,framealpha=1, title='A. Fused silica') 
        #plt.plot([80,80],[0,10000])
        #plt.plot([84,84],[0,10000])
    if pind==125:
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right', fontsize=10,framealpha=1, title='A. Fused silica') 

    if fig26:
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        order = range(len(handles))
        #handles.insert(0,'test')
        #order = [0,1,6,7,8,9,10,3,5,4,2]#2,3,4,5,6,7]
        #order = [0,  5,  6,  7,  8,  9,  2 , 4 , 3 , 1]
        order = [0,5,7,6,8,9,10,2,4,3,1]
        print(len(order))
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right', fontsize=10,framealpha=1,title='B. α-Quartz z-cut') 

    if fig27:
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        order = range(len(handles))
        #handles.insert(0,'test')
        #order = [0,1,2,4,3,9,10,11,8,6,5,7]#  8,9,10,7,5,6]
        #order = [0,  1 , 3 , 2,  8,  9, 10,  7,  5,  4,  6]
        order = [0,1,2,3,8,9,7,5,6,4]
        print(len(order))
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='upper right', fontsize=10,framealpha=1,title='A. Fused silica') 

    if fig28 and not(fig29):
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        order = range(len(handles))
        #handles.insert(0,'test')
        #order = [5,13,3,4,0,6,1,2,11,7,8,9,10,12]
        #order = [0,6,8,7,1,2,3,4,5] # removed 2
#        order = [0,5,7,6,1,2,3,4] #before adding no P annotation
        #order = [0,6,8,7,1,2,3,4,5]
        order = [0,5,6,7,3,4,8,1,2]
        print(len(order))
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right', fontsize=10,framealpha=1,title='Coesite & Stishovite') 

    if fig29:
        plt.legend( fontsize=10,framealpha=1)
    #plt.plot([0,250],[4500,4500])
    if 0: #fig28: # remember cst and stv intersections iwth melt curve
        plt.plot([cst_solidus_P,cst_solidus_P],[0,12000]) # cst enters melt curve
        plt.plot([cst_liquidus_P,cst_liquidus_P],[0,12000]) # cst exits melt curve
        plt.plot([470,470],[0,12000]) # stv enters melt curve
        plt.plot([540,540],[0,12000]) # stv exits melt curve
    if fig30:
        plt.text(10,4500,'LDM')
        plt.text(50,5500,'HDM')
        plt.legend(fontsize=10)     

    if pind == 2425:
        plt.legend('', frameon=False)
    if pind == 2627:
        plt.legend('', frameon=False)
    ################### END LEGEND ###################################

    print('figsize=',plt.rcParams['figure.figsize'])
    try:
        figname
    except NameError:
        print('figname not set! pind=',pind) 
    plt.tight_layout()
    if fig20:
        plt.savefig(outputfigpath+figname) # png file for book logo
    else:
        plt.savefig(outputfigpath+figname+'.pdf')
    if fig23 or fig24 or fig25 or fig26 or fig27 or fig28 or fig29 or fig30:
        plt.savefig(figurefigpath+figname+'.pdf')
    print('Saved figure: '+figname+'.pdf')
    plt.show()
    ################### END PLOT ###################################
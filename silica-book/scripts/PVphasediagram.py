# PVphasediagram.py

def plotomegaRP(hug,label=False):
    omegapmin=85
    omegapmax=111
    mss = 4
    col9='olive'
    zzz = np.where((hug.P > omegapmin) & (hug.P < omegapmax))[0]
    if label:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss,label='CaCl$_2$-type')
    else:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss)#,label='CaCl$_2$-type')
        
    col9='black'
    usmin = 9.23
    usmax = 9.5
    zzz = np.where((hug.us > usmin) & (hug.us < usmax))[0]
    if label:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss,label='Transition')
    else:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss)#,label='Transition')
        
    #mss = 2
    col9='limegreen'
    zzz = np.where((hug.us >= usmax) & (hug.rho <=5.2))[0]
    if label:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss,label='Disordered phases')
    else:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss)#,label='Disordered phases')
        
    rhotrans=5.2
    #mss = 2
    col9='green'
    zzz = np.where(hug.rho > rhotrans)[0]
    if label:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss,label='Liquid')
    else:
        plt.plot(hug.rho[zzz],hug.P[zzz],'.',color=col9,markersize=mss)#,label='Liquid')

print(pind)

if 1:
#for pind in [1,2,3,4,5,6]:
    fig1=False
    fig2=False
    fig3=False
    fig4=False
    fig5=False
    fig6=False
    fig7=False
    fig8=False
    fig9=False
    fig10=False
    fig11=False
    fig12=False
    fig13=False
    fig14=False
    fig15=False
    if pind == 1:
        fig1=True
    if pind == 2:
        fig2=True
    if pind == 3:
        fig3=True
    if pind == 4:
        fig4=True
    if pind == 5:
        fig5=True
    if pind == 6:
        fig6=True
    if pind == 7:
        fig7=True
    if pind == 8:
        fig8=True
        #fig3=True
    if pind == 9:
        fig9=True
    if pind == 10:
        fig10=True
    if pind == 11:
        fig11=True
    if pind == 12:
        fig10 = True
        fig12 = True # add isotherms to Fig 10 
    if pind == 13:
        #fig10 = True
        #fig12 = True # add isotherms to Fig 10 
        fig13=True
    if pind == 14:
        fig14=True
    if pind == 15: ## attempt to shift P instead of rho for DFT offset - does not work
        fig15=True        
    plt.figure(figsize=(10, 6))
    plt.rcParams.update({'font.size': 16})
    #fig = plt.figure(figsize=(8, 6.5))
    #plt.rcParams.update({'font.size': 12})
    # ------------------------------------------------------------------------
    # PLOT PHASE BOUNDARIES
    widthsol  = 3
    widthsol2 = 1
    # plot liquidus
    ind = np.where(stsmcsmooth.rho < 10)[0]
    plt.plot(stsmcsmooth.rho[ind],stsmcsmooth.P[ind],'m-',linewidth=widthsol,label='Model liquidus')
    if fig15:
        plt.plot(stsmcsmooth.rho[ind]-lda_liquid_density_shift,stsmcsmooth.P[ind],'m--',linewidth=widthsol,label='LDA orig')
        plt.plot(stsmcsmooth.rho[ind]-lda_liquid_density_shift,stsmcsmooth.P[ind]+10,'c--',linewidth=widthsol,label='LDA orig + 10 GPa')
        
    #ind = np.where(stsmcsmooth.rho > 5.5)[0]
    #plt.plot([stsmcsmooth.rho[ind][0],6.0],[stsmcsmooth.P[ind][0],295],'m:',linewidth=widthsol,label='Estimated liquidus')
    # plot the triple points
    #ppp = [Pt1, 40, Pt2, Pt3] #40 is LDM to HDM transition from Andrault et al. 2022
    #plt.plot(np.interp(ppp,stsmcsmooth.P,stsmcsmooth.rho),ppp,'mD')
    ppp = [Pt1, Pt3] #40 is LDM to HDM transition from Andrault et al. 2022
    #plt.plot(np.interp(ppp,stsmcsmooth.P,stsmcsmooth.rho),ppp,'mD')
    #ppp = [40, Pt2] #40 is LDM to HDM transition from Andrault et al. 2022
    #plt.plot(np.interp(ppp,stsmcsmooth.P,stsmcsmooth.rho),ppp,'mD',mfc='white')

    #if fig1:
    #    ppp = [Pt1, Pt3] #40 is LDM to HDM transition from Andrault et al. 2022
    #    plt.plot(np.interp(ppp,stsmcsmooth.P,stsmcsmooth.rho),ppp,'mD',label='Est. 1st order triple pt.')
        #ppp = [40, Pt2] #40 is LDM to HDM transition from Andrault et al. 2022
        #plt.plot(np.interp(ppp,stsmcsmooth.P,stsmcsmooth.rho),ppp,'mD',mfc='white',label='Est. 2nd order triple pt.')

    # plot solidus
    ind = np.where((stsmcsmooth.rho_s < 10) & (stsmcsmooth.rho_s > 0))[0]
    plt.plot(stsmcsmooth.rho_s[ind],stsmcsmooth.P[ind],'-',color=sol_col,linewidth=widthsol,label='Model solidus')
    #plt.plot(np.interp(ppp,stsmcsmooth.P,stsmcsmooth.rho_s),ppp,'D',color=sol_col)
    plt.plot(max(p23bnd.rho_b),max(p23bnd.P),'D',color=sol_col)
    plt.plot(max(p23bnd.rho_a),max(p23bnd.P),'D',color=sol_col)
    plt.plot(np.interp(Pt3,stsmcsmooth.P,stsmcsmooth.rho),Pt3,'mD')#,label='Triple point')
    # sft-pyrite triple point
    plt.plot(np.interp(165,stsmcsmooth.P,stsmcsmooth.rho),165,'mD')#,label='Triple point')
    ind = max(np.where(stsmcsmooth.P<165)[0])
    plt.plot(stsmcsmooth.rho_s[ind],165,'D',color=sol_col)
    ind = min(np.where(stsmcsmooth.P>165)[0])
    plt.plot(stsmcsmooth.rho_s[ind],165,'D',color=sol_col)

    # cst-stv-L
    plt.plot(np.interp(Pt1,stsmcsmooth.P,stsmcsmooth.rho),Pt1,'mD')
    #plt.plot(np.interp(Pt1,stsmcsmooth.P,stsmcsmooth.rho)*1.01,Pt1,'D',color=sol_col) # Cst lower density than qtz near TP
    plt.plot(np.interp(Pt1+1,stsmcsmooth.P,stsmcsmooth.rho_s),Pt1,'D',color=sol_col)
    # plot phase boundary CaCl2 to seifertite; calculated in load_solid_models.py
    #plt.plot([rho0_23_sft, rhos_23_sft],[Pt3_0k,Pt3],'--',color=sol_col,linewidth=widthsol2)
    #plt.plot([rho0_23_stv, rhos_23_stv],[Pt3_0k,Pt3],'--',color=sol_col,linewidth=widthsol2)
    # plot phase boundary from stishovite to CaCl2
    # print(rho0_12_stv,rhos_12_stv)
    #plt.plot([rho0_12_stv, rhos_12_stv],[Pt2_0k,Pt2],'--',color=sol_col,linewidth=widthsol2)
    # Stv to coesite
    #plt.plot([rho0_01_stv,rhos_01_stv],[Pt1_0k,Pt1],'--',color=sol_col,linewidth=widthsol2)

    if 0:
        ind = np.where(p23bnd.rho_a> 0)[0]
        plt.plot(p23bnd.rho_a[ind],stsmcsmooth.P[ind],'r-',linewidth=2,label='Calculated mixed hug')
        r2tie = np.interp(111,qtz_hug_cacl2.P,qtz_hug_cacl2.rho)
        r3tie = np.interp(4800,p23bnd.T,p23bnd.rho_b)
        p3tie = np.interp(4800,p23bnd.T,p23bnd.P)
        #plt.plot([r2tie,r3tie],[111,p3tie],'r--',label='Estimated sft Hugoniot')
        #plt.plot([r3tie,np.max(hugrho[ind])],[p3tie,np.min(stsmcsmooth.P[ind])],'r--')
        ind = np.where(p23bnd.rho_a > 0)[0]
        plt.plot(p23bnd.rho_a[ind],p23bnd.P[ind],'r-',linewidth=4,label='Calculated mixed hug')
        #ind = np.where(hugrhofs > 0)[0]
        #plt.plot(hugrhofs[ind],stsmcsmooth.P[ind],'b-',linewidth=4,label='Calculated mixed hug')
        

    #if not(fig1 or fig4 or fig7):
    if 0:
        # multi-wave stv-cacl2 transition region
        ind = np.where( (stsmcsmooth.P >= Pt2a) & (stsmcsmooth.P<= Pt2b))[0]
        temppa = np.zeros(len(ind)+1)+Pt2a
        temppb = np.concatenate((stsmcsmooth.P[ind],np.asarray([Pt2b])))
        tempr = np.concatenate((stsmcsmooth.rho_s[ind],np.asarray([4.835])))
        plt.fill_between(tempr,temppa,temppb,color=sol_col,alpha=0.3,label=r'Stv-CaCl$_2$ shock transition')

    #plt.plot(p12bnd.rho_a,p12bnd.P,'c--',linewidth=widthsol2,label='Stv-CaCl2 bnd')
#    plt.plot(p23bnd.rho_b,p23bnd.P,'--')

    #plt.plot(p23bnd.rho_a,p23bnd.P,'--',linewidth=widthsol2,color=sol_col,label='CaCl2-sft boundary')
    #plt.plot(p23bnd.rho_b,p23bnd.P,'--',linewidth=widthsol2,color=sol_col)

    
    # Plot 0 K isotherms
    amu=sio2mol
    colvinet='black'
    # Vinet
    par = stv
    rho0 = amu*1.6605391/par[0]
    #print(rho0)
    rhoarr = np.arange(rho0,rho0*2,.0001)
    eta = np.power(rho0/rhoarr,0.333)
    parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
#    ind = np.where((parr <= 53) & (parr >= 7))[0]
    ind = np.where((parr <= Pt2_0k) & (parr >= 7))[0]
    #ind = np.where((parr <= 50))[0]
    plt.plot(rhoarr[ind],parr[ind],color=colvinet,linewidth=widthsol,label='T=300K')
    #plt.plot(min(rhoarr[ind]),min(parr[ind]),'kD')
    #plt.plot(max(rhoarr[ind]),max(parr[ind]),'kD')

    # Vinet
    par = cacl2
    rho0 = amu*1.6605391/par[0]
    #print(rho0)
    rhoarr = np.arange(rho0,rho0*2,.0001)
    eta = np.power(rho0/rhoarr,0.333)
    parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
    ind = np.where((parr <= Pt3_0k) & (parr > Pt2_0k))[0]
    plt.plot(rhoarr[ind],parr[ind],color=colvinet,linewidth=widthsol)
    #print(np.max(rhoarr[ind]),np.max(parr[ind]))
    #plt.plot(max(rhoarr[ind]),max(parr[ind]),'kD')

    par = pbo2
    rho0 = amu*1.6605391/par[0]
    #print(rho0)
    rhoarr = np.arange(rho0,rho0*2,.001)
    eta = np.power(rho0/rhoarr,0.333)
    parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
#    ind = np.where((parr >= Pt3_0k) & (rhoarr <= 10))[0]
    ind = np.where((parr >= Pt3_0k) & (parr < Pt4_0k))[0]
    plt.plot(rhoarr[ind],parr[ind],color=colvinet,linewidth=widthsol)
    #plt.plot(min(rhoarr[ind]),min(parr[ind]),'kD')
    #print(np.min(rhoarr[ind]),np.min(parr[ind]))
    #print('PbO2 0k rho = ',rhoarr[ind[0]])

    # add pyrite from Das et al.
    par = pyrite
    rho0 = sio2mol/12.97 # g/cm3
    rhoarr = np.arange(rho0,rho0*2,.001)
    eta = np.power(rho0/rhoarr,0.333)
    parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
    ind = np.where((parr >= Pt4_0k) & (rhoarr <= 10))[0]
    plt.plot(rhoarr[ind],parr[ind],color=colvinet,linewidth=widthsol)    
# ------------------------------------------------------------------------
    if 0:
        plt.plot([4,6],[63,63])
    # plot wiggle region
    if fig10:
        #plt.fill_between([4.97,5.07],[qtz_wiggle_P1,qtz_wiggle_P1],[qtz_wiggle_P2,qtz_wiggle_P2],color='grey',alpha=0.3,label='Temp. #inversion')
        plt.fill_between([qtz_wiggle_rho1,qtz_wiggle_rho2],[qtz_wiggle_P1,qtz_wiggle_P1],[qtz_wiggle_P2,qtz_wiggle_P2],color='grey',alpha=1,label='Temp. inversion')
        #plt.plot([0,10],[141,141])
    if fig9:
#        plt.fill_between([4.5,4.7],[fs_wiggle_P1,fs_wiggle_P1],[fs_wiggle_P2,fs_wiggle_P2],color='grey',alpha=0.3,label='Temp. inversion')
        plt.fill_between([fs_wiggle_rho1,fs_wiggle_rho2],[fs_wiggle_P1,fs_wiggle_P1],[fs_wiggle_P2,fs_wiggle_P2],color='grey',alpha=0.7,label='Temp. inversion')
        #plt.plot([0,10],[80,80])
                                
    
    # ------------------------------------------------------------------------
    # FUSED SILICA SHOCK HUGONIOT DATA
    if fig3:
        plt.errorbar(fsdf['rho'][indfs],fsdf['P'][indfs],fmt='b+',label='RB FS')

    if fig2 or fig3 or fig5 or fig9:
        #plt.errorbar(shockdat_list[L83fshug['id']].rho,shockdat_list[L83fshug['id']].P,fmt='bo',mfc='white',\
        #             yerr=shockdat_list[L83fshug['id']].P_std,label=L83fshug['lab'])
        plt.errorbar(shockdat_list[W62fs['id']].rho,shockdat_list[W62fs['id']].P,fmt='bx',label='W62',\
                     #W62fs['lab'],\
                     xerr=shockdat_list[W62fs['id']].rho_std,yerr=shockdat_list[W62fs['id']].P_std)
        plt.errorbar(shockdat_list[J68fs['id']].rho,shockdat_list[J68fs['id']].P,fmt='v',color='blue',mfc='lightblue',label='J68',\
                     #J68fs['lab'],\
                     xerr=shockdat_list[J68fs['id']].rho_std,yerr=shockdat_list[J68fs['id']].P_std)
        plt.errorbar(shockdat_list[B88RRfs['id']].rho,shockdat_list[B88RRfs['id']].P,fmt='bo',mfc='lightblue',yerr=shockdat_list[B88RRfs['id']].P_std,\
                     xerr=shockdat_list[B88RRfs['id']].rho_std,label='L83R+B88R')#'+B88RRfs['lab'])

        #plt.errorbar(shockdat_list[T18fs['id']].rho,shockdat_list[T18fs['id']].P,fmt='b<',mfc='white',\
        #             xerr = shockdat_list[T18fs['id']].rho_std, yerr=shockdat_list[T18fs['id']].P_std, label=T18fs['lab'])
        #plt.errorbar(shockdat_list[S80fs['id']].rho,shockdat_list[S80fs['id']].P,fmt='bv',xerr=shockdat_list[S80fs['id']].rho_std,yerr=shockdat_list[S80fs['id']].P_std,label=S80fs['lab'])
        plt.errorbar(shockdat_list[T18fs['id']].rho,shockdat_list[T18fs['id']].P,fmt='<',color='blue',mfc='white',label='T18',\
                     #T18fs['lab'],\
                     xerr=shockdat_list[T18fs['id']].rho_std,yerr=shockdat_list[T18fs['id']].P_std)
        plt.errorbar(shockdat_list[R19fs['id']].rho,shockdat_list[R19fs['id']].P,fmt='bd',\
                 yerr=shockdat_list[R19fs['id']].P_std,xerr=shockdat_list[R19fs['id']].rho_std,label='R19')#R19fs['lab'])
        plt.errorbar(shockdat_list[O25fs['id']].rho,shockdat_list[O25fs['id']].P,fmt='b^',\
                 yerr=shockdat_list[O25fs['id']].P_std,xerr=shockdat_list[O25fs['id']].rho_std,label='O25')#O25fs['lab'])
        plt.errorbar(shockdat_list[B19fs['id']].rho,shockdat_list[B19fs['id']].P,fmt='b>',mfc='white',label='B19',\
                     #B19fs['lab'],\
                     xerr=shockdat_list[B19fs['id']].rho_std,yerr=shockdat_list[B19fs['id']].P_std)
        #if not(fig5):
            #plt.errorbar(shockdat_list[A23fs['id']].rho,shockdat_list[A23fs['id']].P,fmt='bD',mfc='lightblue',xerr=shockdat_list[A23fs['id']].rho_std,yerr=shockdat_list[A23fs['id']].P_std,label='This work FS')
                        #A23fs['lab'])
#        plt.errorbar(shockdat_list[J79fs['id']].rho,shockdat_list[J79fs['id']].P,fmt='r<',yerr=shockdat_list[J79fs['id']].P_std,label=J79fs['lab'])
        #plt.errorbar(shockdat_list[B88RRfs['id']].rho,shockdat_list[B88RRfs['id']].P,fmt='bo',mfc='lightblue',yerr=shockdat_list[B88RRfs['id']].P_std,label=B88RRfs['lab'])
        #plt.errorbar(shockdat_list[J79Rfs['id']].rho,shockdat_list[J79Rfs['id']].P,fmt='^',color='grey',mfc='grey',label=J79Rfs['lab'],\
        #             xerr=shockdat_list[J79Rfs['id']].rho_std,yerr=shockdat_list[J79Rfs['id']].P_std,ecolor='grey')
        #plt.errorbar(shockdat_list[A23fs['id']].rho,shockdat_list[A23fs['id']].P,fmt='bD',mfc='lightblue',\
        #         xerr=shockdat_list[A23fs['id']].rho_std,yerr=shockdat_list[A23fs['id']].P_std,label=A23fs['lab'])
        indA23fs = np.arange(0,9) # EOS points with error bars
#        plt.errorbar(np.asarray(shockdat_list[A23fs['id']].rho[indA23fs]),shockdat_list[A23fs['id']].P[indA23fs],fmt='bD',mfc='lightblue',\
#            xerr=shockdat_list[A23fs['id']].rho_std[indA23fs],yerr=shockdat_list[A23fs['id']].P_std[indA23fs], alpha = 0.4,label=A23fs['lab'])
        plt.errorbar(np.asarray(shockdat_list[A23fs['id']].rho[indA23fs]),shockdat_list[A23fs['id']].P[indA23fs],fmt='bD',mfc='lightblue',\
            xerr=shockdat_list[A23fs['id']].rho_std[indA23fs],yerr=shockdat_list[A23fs['id']].P_std[indA23fs], alpha = 0.4)#,label='This work')
        plt.plot(np.asarray(shockdat_list[A23fs['id']].rho),shockdat_list[A23fs['id']].P,'bD',mfc='lightblue',label='This work')
        

    # ------------------------------------------------------------------------
    if fig9 or fig10 or fig15:
        res = np.polyfit(eosdat_list[G24mclrho['id']].rho_liquid,eosdat_list[G24mclrho['id']].P_liquid,2)
        #plt.plot(eosdat_list[G24mclrho['id']].rho_liquid+.13,eosdat_list[G24mclrho['id']].P_liquid,'-',color='red',label=G24mclrho['lab']+' (2phase, shift)')

        xrho = np.arange(3.5,6.75,.01)
        ypres = res[2]+res[1]*xrho+res[0]*xrho*xrho
        #plt.plot(xrho+pbe_liquid_density_shift,ypres,'-',color='red',label=G24mclrho['lab']+'\n(2phase, shift)')
        #plt.plot(xrho,ypres-10,'-',color='cyan',label=G24mclrho['lab']+'\n(2phase, Pshift)')
        if fig15:
            plt.plot(xrho,ypres,'--',color='red',label=G24mclrho['lab']+'\n(2phase, orig)')

        #plt.plot([0,10],[fs_startliq_P,fs_startliq_P],label='fs PT cross')
        #plt.plot([0,10],[qtz_startliq_P,qtz_startliq_P],label='qtz PT cross')

    # ------------------------------------------------------------------------
    
    # QUARTZ SHOCK HUGONIOT DATA

    if fig2 or fig3 or fig5 or fig10:
        #plt.errorbar(qtzdf['rho'][indqtz],qtzdf['P'][indqtz],fmt='g+',label='M80 qtz')
        #plt.errorbar(shockdat_list[L83qtzhug['id']].rho,shockdat_list[L83qtzhug['id']].P,fmt='go',mfc='white',\
        #             yerr=shockdat_list[L83qtzhug['id']].P_std,label=L83qtzhug['lab'])
        ind=np.where(np.asarray(shockdat_list[B88RRqtz['id']].orientation) == 'z')[0]
        #plt.errorbar(shockdat_list[B88RRqtz['id']].rho[ind],shockdat_list[B88RRqtz['id']].P[ind],fmt='go',mfc='yellow',yerr=shockdat_list[B88RRqtz['id']].P_std[ind],label=B88RRqtz['lab']+' (z)')
        print('A=',np.asarray(shockdat_list[B88RRqtz['id']].P)[ind])
        ind=np.where(np.asarray(shockdat_list[B88RRqtz['id']].orientation) == 'pz')[0]
        #plt.errorbar(shockdat_list[B88RRqtz['id']].rho[ind],shockdat_list[B88RRqtz['id']].P[ind],fmt='go',mfc='orange',yerr=shockdat_list[B88RRqtz['id']].P_std[ind],label=B88RRqtz['lab']+'-pz')
        print('B=',np.asarray(shockdat_list[B88RRqtz['id']].P)[ind])
        for i in [2]:
            print(shockdat_list[B88RRqtz['id']].P[i],shockdat_list[B88RRqtz['id']].orientation[i])
            plt.errorbar(shockdat_list[B88RRqtz['id']].rho[i],shockdat_list[B88RRqtz['id']].P[i],fmt='go',mfc='yellow',\
                     yerr=shockdat_list[B88RRqtz['id']].P_std[i],xerr=shockdat_list[B88RRqtz['id']].rho_std[i],label='L83R+B88R')#+B88RRqtz['lab'])
        for i in [4]:
            print(shockdat_list[B88RRqtz['id']].P[i],shockdat_list[B88RRqtz['id']].orientation[i])
            plt.errorbar(shockdat_list[B88RRqtz['id']].rho[i],shockdat_list[B88RRqtz['id']].P[i],fmt='go',mfc='yellow',\
                     yerr=shockdat_list[B88RRqtz['id']].P_std[i],xerr=shockdat_list[B88RRqtz['id']].rho_std[i])#,label=B88RRqtz['lab']+'-IM')

        ind = np.where(shockdat_list[W62qtz['id']].orientation==' z')[0]
        plt.errorbar(shockdat_list[W62qtz['id']].rho[ind],shockdat_list[W62qtz['id']].P[ind],fmt='gx',label='W62')#W62qtz['lab'])

#plt.errorbar(shockdat_list[W62qtz['id']].rho,shockdat_list[W62qtz['id']].P,fmt='gx',label=W62qtz['lab'])
        #plt.errorbar(shockdat_list[L03qtz['id']].rho,shockdat_list[L03qtz['id']].P,fmt='r>',label=L03qtz['lab'])
        #plt.errorbar(shockdat_list[T20qtz['id']].rho,shockdat_list[T20qtz['id']].P,fmt='g<',mfc='white',\
        #             xerr = shockdat_list[T20qtz['id']].rho_std, yerr=shockdat_list[T20qtz['id']].P_std, label=T20qtz['lab'])
        #plt.errorbar(shockdat_list[T20qtz['id']].rho,shockdat_list[T20qtz['id']].P,fmt='<',color='grey',\
        #             xerr = shockdat_list[T20qtz['id']].rho_std, yerr=shockdat_list[T20qtz['id']].P_std, label=T20qtz['lab'])
        #print(shockdat_list[T20qtz['id']].up,shockdat_list[T20qtz['id']].us)
        #res = np.polyfit(shockdat_list[T20qtz['id']].up,shockdat_list[T20qtz['id']].us,1)
        #print(res)
        #plt.errorbar(shockdat_list[L04qtz['id']].rho,shockdat_list[L04qtz['id']].P,fmt='gv',mfc='white',\
        #         xerr=shockdat_list[L04qtz['id']].rho_std,yerr=shockdat_list[L04qtz['id']].P_std,label=L04qtz['lab'])
        if not(fig5):
            ind = np.where(shockdat_list[A23qtz['id']].P>qtz_singleshock_P)[0]
            plt.errorbar(shockdat_list[A23qtz['id']].rho[ind],shockdat_list[A23qtz['id']].P[ind],fmt='gD',mfc='yellow',xerr=shockdat_list[A23qtz['id']].rho_std[ind],yerr=shockdat_list[A23qtz['id']].P_std[ind],label='This work')
                        #A23qtz['lab']+' (z)')


    # ------------------------------------------------------------------------

    # plot new best fits
    if fig9: # FS
        ind = np.where(fs_hug_liq.P > fs_wiggle_P2)[0]
        plt.plot(fs_hug_liq.rho[ind],fs_hug_liq.P[ind],'b--',label='New Hugoniot: Liquid region')
        ind = np.where((fs_hug_stv.P > 34) & (fs_hug_stv.P < fs_wiggle_P1))[0]
        plt.plot(fs_hug_stv.rho[ind],fs_hug_stv.P[ind],'b-',label='New Hugoniot: Solid region')

    #if fig13:
    #    plt.plot(cst_hug_F24.rho,cst_hug_F24.P,'-',color='salmon',linewidth=4,label='F24 Cst Hug: solid region')
    
    if fig10: # Qtz
        ind = np.where(qtz_hug_liq.P > qtz_wiggle_P2)[0]
        plt.plot(qtz_hug_liq.rho[ind],qtz_hug_liq.P[ind],'g--',label='New Hugoniot Liquid region')
        ind = np.where((qtz_hug_cacl2.P > fs_wiggle_P1) & (qtz_hug_cacl2.P < qtz_wiggle_P1))[0]
        #plt.plot(qtz_hug_cacl2.rho[ind],qtz_hug_cacl2.P[ind],'g:',label='New Qtz Hug: CaCl$_2$ mc region')
        #ind = np.where((qtz_hug_stv.P > 34) & (qtz_hug_stv.P < fs_wiggle_P1))[0]
        #plt.plot(qtz_hug_stv.rho[ind],qtz_hug_stv.P[ind],'g-',label='New Qtz Hug: Stv region')
        ind = np.where((qtz_hug_solid.P > 60) & (qtz_hug_solid.P < qtz_wiggle_P1))[0]
        plt.plot(qtz_hug_solid.rho[ind],qtz_hug_solid.P[ind],'g-',label='New Hugoniot: Solid region')
    if fig13:
        uptrans1 = 4.089998285871261 # km/s # max solid
        uptrans2 = 4.950818032171942 # km/s # start opaque
        uptrans3 = 5.603566588729885 # km/s # start liquid
        #plt.fill_between([5.,5.35], [113.13,113.13],[159,159],color='grey',alpha=0.2,label='Temp. inversion')
        #plt.fill_between([5.14,5.18], [114,114],[134,134],color='grey',alpha=0.2,label='Temp. inversion')
        plt.fill_between([cst_wiggle_rho1,cst_wiggle_rho2], [cst_wiggle_P1,cst_wiggle_P1],[cst_wiggle_P2,cst_wiggle_P2],color='grey',alpha=1,label='Temp. inversion')
        # Coesite fit Feng et al. 2024
        #        Us = 2.097(±0.099) + 1.815(±0.022) Up – 0.020(±0.001) U2p, Us ≥ 10.4 km/s.
        ind = np.where(cst_hug_F24.us > 10.75)
        #plt.plot(cst_hug_F24.rho[ind],cst_hug_F24.P[ind],'--', color='darkkhaki',linewidth=3,label='F24 Cst Hug: Liquid region')
        #plt.plot(eosdat_list[G24mclrho['id']].rho_liquid,eosdat_list[G24mclrho['id']].P_liquid,'-',color='red',label=G24mclrho['lab']+' (two-phase)')

        #res = np.polyfit(eosdat_list[G24mclrho['id']].rho_liquid,eosdat_list[G24mclrho['id']].P_liquid,2)
        #xrho = np.arange(3.5,6.75,.01)
        #ypres = res[2]+res[1]*xrho+res[0]*xrho*xrho
        #plt.plot(xrho+gengshift,ypres,'-',color='red',label=G24mclrho['lab']+'\n(2phase, shift)')

        
        ind = np.where(cst_hug_liq.P > cst_wiggle_P2)
        plt.plot(cst_hug_liq.rho[ind],cst_hug_liq.P[ind],'--', color='salmon',linewidth=3,label='New Hugoniot: Liquid region')
        ind = np.where(cst_hug_cacl2.P < cst_solidus_P)
        plt.plot(cst_hug_cacl2.rho[ind],cst_hug_cacl2.P[ind],'-',color='salmon',linewidth=3,label='New Hugoniot: Solid region')
        #plt.plot(np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.rho),cst_liquidus_P,'kD')
        #ind = np.where(qtz_hug_cacl2.P<qtz_wiggle_P1)
        #ind = np.where((qtz_hug_cacl2.P > fs_wiggle_P1) & (qtz_hug_cacl2.P < qtz_wiggle_P1))[0]
        #plt.plot(qtz_hug_cacl2.rho[ind],qtz_hug_cacl2.P[ind],':',color='green',linewidth=2,label='New Qtz Hug: CaCl$_2$ mc region')
        #plt.plot(qtz_hug_stv.rho,qtz_hug_stv.P,'-',color='green',linewidth=2,label='New Qtz Hug: Stv region')
        #ind = np.where(C43939.up>uptrans2)
        #plt.plot(C43939.rho[ind],C43939.P[ind],'.',color='red',label='Decay: opaque front')
        #ind = np.where(C43943.up>uptrans2)
        #plt.plot(C43943.rho[ind],C43943.P[ind],'.',color='red')
        #ind = np.where((C43939.up>uptrans1) & (C43939.up<uptrans2))[0]
        #plt.plot(C43939.rho[ind],C43939.P[ind],'.',color='black')#,label='Decay: temp. inversion') 
        #plt.plot(C43939.rho[ind][0],C43939.P[ind][0],'.',color='black',markersize=10,label='Decay: temp. inversion')
        #ind = np.where(C43939.up<uptrans1)
        #plt.plot(C43939.rho[ind],C43939.P[ind],'.',color='gold')#,label='Decay: transparent front')
        #plt.plot(C43939.rho[ind][0],C43939.P[ind][0],'.',color='gold',markersize=10,label='Decay: transparent front')
        #ind = np.where(C43943.up<uptrans1)
        #plt.plot(C43943.rho[ind],C43943.P[ind],'.',color='gold')#,label='Decay: transparent front')
        # plot Cst data near transition
        #139.593	1.800 	5.2570
        #111.165	1.757 	5.1305
        #149 ± 7 5.45 ± 0.26
        #plt.plot([	4.8057, 4.9912, 5.0758, 5.2570],[66.529,94.261,98.235,139.593],'X',color='tomato',markersize=8,label='L02')
        #plt.plot([	3.2540, 3.4866,3.6795,4.0633,4.7595,5.0195,5.1305],[12.392,22.137,30.986,43.187,61.799,84.545,111.165],'+',color='tomato',markersize=8,label='P81')
        #plt.plot([4.7595,5.0195,5.1305],[61.799,84.545,111.165],'d',color='tomato',markersize=8,label='P81')
        #cst_startliq_P = 190
        #cst_startliq_rho = np.interp(cst_startliq_P,cst_hug_liq.P,cst_hug_liq.rho)
        #plt.plot(cst_startliq_rho,cst_startliq_P,'k^',mfc='tomato',markersize=12,label='Cst Hug transition to liquid field')
        #ind = np.where(C43943.us<9.5)
        #plt.plot(C43943.rho[ind],C43943.P[ind],'.',color='gold')

        
        plt.errorbar(shockdat_list[P81cst['id']].rho,shockdat_list[P81cst['id']].P,fmt='kX',label='P81',\
                     #P81cst['lab'],\
                    xerr=shockdat_list[P81cst['id']].rho_std,yerr=shockdat_list[P81cst['id']].P_std)
        plt.errorbar(shockdat_list[L02cst['id']].rho,shockdat_list[L02cst['id']].P,fmt='ko',mfc='white',label='L02',\
                     #L02cst['lab'],\
                    xerr=shockdat_list[L02cst['id']].rho_std,yerr=shockdat_list[L02cst['id']].P_std)
        plt.errorbar(shockdat_list[F24cst['id']].rho,shockdat_list[F24cst['id']].P,fmt='kD',mfc='white',label='F24',\
                     #F24cst['lab'],\
                    xerr=shockdat_list[F24cst['id']].rho_std,yerr=shockdat_list[F24cst['id']].P_std)

        #plt.plot(eosdat_list[G24mclrho['id']].rho_liquid,eosdat_list[G24mclrho['id']].P_liquid,'-+',color='red',label=G24mclrho['lab'])
        plt.plot(cst_liquidus_rho,cst_liquidus_P,'kx',label='PT liquidus intersection')
        
        
 
    # PLOT Z HUGONIOTS
    if fig2 or fig3 or fig5 or fig8:
        ind = np.where(qtz_hug_KD.P > 140)[0]
        plt.plot(qtz_hug_KD.rho[ind],qtz_hug_KD.P[ind],'g--',label='Qtz Z fit',linewidth=1)
        plt.plot(qtz_hug_liq.rho,qtz_hug_liq.P,'g-',label='Qtz liq NEW')
        #tmpx = np.interp(111,qtz_hug_KD.P,qtz_hug_KD.rho)
        #plt.plot(tmpx,111,'rP')
        #ind = np.where((qtz_hug_KD.P > 117) & (qtz_hug_KD.P < 150.))[0]
        #plt.plot(qtz_hug_KD.rho[ind],qtz_hug_KD.P[ind],'g--')

        #ind = np.where((qtz_hug_KD.P > 117) & (qtz_hug_KD.P < 150.))[0]
        # manual estimate of PbO2 Hugoniot
        #plt.plot([5.18,5.25],[118,150],'g--')

        if not(fig5):
            #ind = np.where(fs_hug_RTK19.P > 74.)[0]
            #plt.plot(fs_hug_RTK19.rho[ind],fs_hug_RTK19.P[ind],'b--',label='FS Z fit',linewidth=1)
            plt.plot(fs_hug_liq.rho,fs_hug_liq.P,'b-',label='FS liq NEW')
        #ind = np.where(fs_hug_RTK19.P > 72.)[0]
        #plt.plot(fs_hug_RTK19.rho[ind],fs_hug_RTK19.P[ind],'b--',linewidth=2)
        #plt.plot([4.71,4.8],[74,110],'r--',linewidth=2)
        # following liquidus region?
        #plt.plot(liquidusrho[4:6],  liquidusP[4:6],'b--',linewidth=3)

    #if 0:#fig1:
        #plt.plot(M70813.rho,M70813.P,'r--',linewidth=2,label='Omega Hugoniot')
        #plt.plot(M70813.rho,M70813.P,'g--',linewidth=3,label='Est. Qtz Hugoniot')
        #plt.plot([5.2,5.2],[0,250])
        #plt.plot([5.22,5.22],[0,250])
        #plt.plot([0,10],[145,145])
#        plt.plot([0,10],[140,140])
        #plt.plot(qtz_hug_liq.rho,qtz_hug_liq.P,'g-',label='Qtz liq NEW')

    if fig2 or fig3 or fig5 or fig6 or fig8:
        plt.plot(fs_hug_stv.rho,fs_hug_stv.P,'-',color='blue',linewidth=2,label='FS Hug NEW')
        #plt.plot(np.interp(fs_liquidus_P,fs_hug_liq.P,fs_hug_liq.rho),fs_liquidus_P,'kD')

    if fig2 or fig3 or fig5 or fig6 or fig8:
        plt.plot(qtz_hug_stv.rho,qtz_hug_stv.P,'-',color='green',linewidth=2,label='Qtz Hug NEW')
        ind = np.where((qtz_hug_cacl2.P<112))[0]
        #plt.plot(qtz_hug_cacl2.rho[ind],qtz_hug_cacl2.P[ind],'-',color='green',linewidth=2)
        #plt.plot(np.interp(qtz_liquidus_P,qtz_hug_liq.P,qtz_hug_liq.rho),qtz_liquidus_P,'kD')
        #ind = np.where((qtz_hug_sft.P>112) & (qtz_hug_sft.P<140))[0]
        #plt.plot(qtz_hug_sft.rho[ind],qtz_hug_sft.P[ind],'r--',linewidth=2,label='Qtz Hug Sft NEW')

    if fig10:
        plt.errorbar(shockdat_list[K13qtz['id']].rho,shockdat_list[K13qtz['id']].P,fmt='gd',markersize=6,\
             xerr=shockdat_list[K13qtz['id']].rho_std,yerr=shockdat_list[K13qtz['id']].P_std,label='K13')#K13qtz['lab'])
    if fig5:
        # add knudson data
        plt.errorbar(shockdat_list[K13qtz['id']].rho,shockdat_list[K13qtz['id']].P,fmt='gd',markersize=6,\
             xerr=shockdat_list[K13qtz['id']].rho_std,yerr=shockdat_list[K13qtz['id']].P_std,label=K13qtz['lab'])
        # plot new qtz data points from UCD gas gun lab
        plt.plot([4.945,5.03885],[84.709,102.44],'rP',label='Shots 47,48 MGR')
        for iii in range(len(qtz_list)):
            plt.errorbar(qtz_list[iii].rho,qtz_list[iii].P,fmt='D',mfc='yellow',markersize=6,\
                xerr=qtz_list[iii].rho_std,yerr=qtz_list[iii].P_std)#,label='UCD')
        for iii in range(len(fs_list)):
            plt.errorbar(fs_list[iii].rho,fs_list[iii].P,fmt='D',mfc='cyan',markersize=6,\
                xerr=fs_list[iii].rho_std,yerr=fs_list[iii].P_std)#,label='UCD')
        # Altshuler
        palt = 2.65*np.asarray([3.13,3.92,6.2])*np.asarray([7.18,8.54,12.01])
        #print('palt = ',palt)
        rhoalt = 1./((1./2.65)*(1.-np.asarray([3.13,3.92,6.2])/np.asarray([7.18,8.54,12.01])))
        plt.plot(rhoalt,palt,'gP',label='qtz Altshuler 1965')

    # plot ab initio information
    if fig1:

        ind = np.where(eosdat_list[U10mcl['id']].rho_liquid < 6.5)[0]
        #plt.plot(eosdat_list[U10mcl['id']].rho_liquid[ind],eosdat_list[U10mcl['id']].P_liquid[ind],'+',color='grey',label=U10mcl['lab'])
        plt.plot(eosdat_list[U10mcl['id']].rho_liquid[ind]+lda_liquid_density_shift,eosdat_list[U10mcl['id']].P_liquid[ind],'k+',label=U10mcl['lab']+' shift')
        ind = np.where(eosdat_list[U10mcl['id']].rho_solid < 6.5)[0]
        #plt.plot(eosdat_list[U10mcs['id']].rho_solid[ind],eosdat_list[U10mcs['id']].P_solid[ind],'x',color='grey',label=U10mcs['lab'])
        plt.plot(eosdat_list[U10mcs['id']].rho_solid[ind]+lda_liquid_density_shift,eosdat_list[U10mcs['id']].P_solid[ind],'kx',label=U10mcs['lab']+' shift')

        #plt.errorbar(eosdat_list[G16mc['id']].rho_l[1:4],eosdat_list[G16mc['id']].P[1:4],fmt='d--',mfc='white',color='#FFA500',label=G16mc['lab']) # pyrite structure
        #plt.plot([3.731,4.199,4.55],  [14.75,37.45,60.1],'d',color='#FFA500',label='B95 stv MD')


    if fig1:
        gengshift=pbe_liquid_density_shift

        #ind = np.where(eosdat_list[G24mcs['id']].P_solid<200)        
        ind = np.where((eosdat_list[G24mcs['id']].rho_solid<5.6)&(eosdat_list[G24mcs['id']].P_solid<165))
        plt.plot(eosdat_list[G24mcs['id']].rho_solid[ind]+gengshift,eosdat_list[G24mcs['id']].P_solid[ind],'x',color='red',label='G24 solid PBE shift')#label=G24mcs['lab']+' shift')

        ind = np.where(eosdat_list[G24mcl['id']].P_liquid<200)
        plt.plot(eosdat_list[G24mcl['id']].rho_liquid[ind]+gengshift,eosdat_list[G24mcl['id']].P_liquid[ind],'+',color='red',label=G24mcl['lab']+ ' shift')        #plt.plot(eosdat_list[G24mclrho['id']].rho_liquid,eosdat_list[G24mclrho['id']].P_liquid,'-',color='red',label=G24mclrho['lab']+' (two-phase)')
        res = np.polyfit(eosdat_list[G24mclrho['id']].rho_liquid,eosdat_list[G24mclrho['id']].P_liquid,2)
        plt.plot(eosdat_list[G24mclrho['id']].rho_liquid+pbe_liquid_density_shift,eosdat_list[G24mclrho['id']].P_liquid,'o',markersize=6,color='red',label=G24mclrho['lab']+' shift')

        xrho = np.arange(3.5,6.75,.01)
        ypres = res[2]+res[1]*xrho+res[0]*xrho*xrho
        #plt.plot(xrho+gengshift,ypres,'-',color='red',label=G24mclrho['lab']+'\n(2phase, shift)')

        # Park & Deng 2025 PBEsol and SCAN
        ind = np.asarray([0,1,2,3,9,10,11,12]) #PBEsol
        plt.plot(eosdat_list[P25mcpbe['id']].rho_l[ind],eosdat_list[P25mcpbe['id']].P[ind],':',color='m',label='P25 liquidus PBEsol')#P25mcpbe['lab'])
        ind = np.asarray([0,1,2,3])#PBEsol
        plt.plot(eosdat_list[P25mcpbe['id']].rho_s[ind],eosdat_list[P25mcpbe['id']].P[ind],':',color=sol_col,label='P25 solidus PBEsol')#P25mcpbe['lab'])
        ind = np.asarray([9,10,11,12])#PBEsol
        plt.plot(eosdat_list[P25mcpbe['id']].rho_s[ind],eosdat_list[P25mcpbe['id']].P[ind],':',color=sol_col)#P25mcpbe['lab'])
        ind = np.asarray([1,2,3,4,5,10,11,12,13])#SCAN

        plt.plot(eosdat_list[P25mcscan['id']].rho_l[ind],eosdat_list[P25mcscan['id']].P[ind],'-.',color='m',label='P25 liquidus SCAN')#P25mcscan['lab'])
#        plt.plot(eosdat_list[P25mcscan['id']].rho_l[ind],eosdat_list[P25mcscan['id']].P[ind],'--',color='plum',label='P25 liquidus SCAN')#P25mcscan['lab'])

        
        for icn in range(len(eosdat_list[G24mclcn['id']].CN_liquid)):
            cn = str(eosdat_list[G24mclcn['id']].CN_liquid[icn])[0:3]
            pp = eosdat_list[G24mclcn['id']].P_liquid[icn]
            irho = np.interp(pp,eosdat_list[G24mclrho['id']].P_liquid,eosdat_list[G24mclrho['id']].rho_liquid)#+0.13)
            if pp < 250:
                plt.text(irho-.06,pp,cn,color='red',fontsize=12,horizontalalignment='center', verticalalignment='center')

    # ------------------------------------------------------------------------
    # arrows and labels
    if 0: #not(fig6 or fig7 or fig1 or fig9 or fig10 or fig13):
        fs=9
        plt.text(5.72,5,'cr/h-qtz/cs',fontsize=fs)
        plt.annotate('',xy=(5.7,0), xytext=(5.7,Pt1), arrowprops=dict(arrowstyle='<-'))
        plt.text(5.72,40,'stishovite',fontsize=fs)
        plt.annotate('',xy=(5.7,Pt1), xytext=(5.7,Pt2), arrowprops=dict(arrowstyle='<->'))
        plt.text(5.72,105,r'CaCl$_2$-type',fontsize=fs)
        plt.annotate('',xy=(5.7,Pt2), xytext=(5.7,Pt3), arrowprops=dict(arrowstyle='<->'))
        plt.text(5.72,160,r'seifertite',fontsize=fs)
        plt.annotate('',xy=(5.7,Pt3), xytext=(5.7,Pt4), arrowprops=dict(arrowstyle='<->'))
        plt.text(5.72,195,'Melt\ncurve\nphase\npyrite-type',fontsize=fs)
        plt.annotate('',xy=(5.7,Pt4), xytext=(5.7,225), arrowprops=dict(arrowstyle='->'))


    if fig7:
        fs=10
        plt.text(4.7,52,'stishovite',color='black',fontsize=fs)
        plt.text(4.6,120,'liquid',color='black',fontsize=fs)
        plt.text(5.25,107,r'$\alpha$PbO$_2$',color='black',fontsize=fs)
        plt.text(5.01,82,r'CaCl$_2$',horizontalalignment='left',color='black',fontsize=fs)


    ## add Andrault stv melt curve estimates 2020
    sio2mol = 28.085+ 2*15.999
    #print(sio2mol/14.6916,sio2mol/19.459)
    #plt.plot([sio2mol/14.6916,sio2mol/19.459],[15,15],'rd',label='A20 stv')
    #plt.plot([sio2mol/13.6434,sio2mol/14.2008,sio2mol/14.3883],[45,45,45],'rd')
    #print(sio2mol/13.6434,sio2mol/14.2008,sio2mol/14.3883)

    if 0:
        ## add Andrault stv melt curve estimates 2020
        sio2mol = 28.085+ 2*15.999
        print(sio2mol/14.6916,sio2mol/19.459)
        plt.plot([sio2mol/14.6916,sio2mol/19.5],[15,15],'rd',label='A22 stv')
        plt.plot([sio2mol/13.63,sio2mol/14.1,sio2mol/14.2],[40,40,40],'rd')
        print('A2022: ',sio2mol/13.63,sio2mol/14.1,sio2mol/14.2)


    # Kraus et al. 2012 Entropy Tie Points for SiO2
    # density T
    # plt.plot([1/0.2125,1/0.2125,1/0.2153,1/0.1983],[334,4500,4500,4800],'co',label='K12 S tie pts')
    # P-T
    # plt.plot([1.e-4,70,70,113],[334,4500,4500,4800],'co',label='K12 S tie pts')
    # rho-P
    #plt.plot([1/0.2125,1/0.2125,1/0.2153,1/0.1983],[1.e-4,70,70,113],'co',label='K12 S tie pts')
    #print('kraus',1/0.2125,1/0.1983)


    if 0:
        par = pyrite
        rho0 = amu*1.6605391/par[0]
        #print(rho0)
        rhoarr = np.arange(rho0,rho0*2.2,.1)
        eta = np.pow-er(rho0/rhoarr,0.333)
        parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
        #plt.plot(rhoarr,parr,label='pyrite')




    # ------------------------------------------------------------------------
    if fig1:
        # Si-O coordination numbers
        # Usui 2020 SiO RCN
        #cncol='#E50000'
        cncol='black'
        fs=10
        yshift=15
        # shifted because LDA DFT to match the liquid model from Karki LDA DFT
        #plt.text(3.731+lda_liquid_density_shift,35,'Si-O CN\n3900K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(3.731+lda_liquid_density_shift,23+yshift,'5.1',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.199+lda_liquid_density_shift,50+yshift,'5.6',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        #plt.text(4.199+lda_liquid_density_shift,56,'5200K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.549+lda_liquid_density_shift,77+yshift,'5.6',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        #plt.text(4.549+lda_liquid_density_shift,83,'5400K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.942+lda_liquid_density_shift,108+yshift,'6.3',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        #plt.text(4.942+lda_liquid_density_shift,114,'5800K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(3.25,40,'LDM',fontsize=12)
#        plt.text(3.875,53,'HDM',fontsize=fs)
        plt.text(4.,80,'HDM',fontsize=12)

    if fig1:
        tss = 11
        plt.text(4.,150,'liquid',fontsize=14)
        #plt.text(3.69,45,'LDM',fontsize=tss)
        #plt.text(3.875,53,'HDM',fontsize=tss)
        plt.text(4.58,44,'stv',fontsize=tss)
        plt.text(4.95,90,'CaCl₂',fontsize=tss)
        plt.text(5.45,135,'sft',fontsize=tss)
        plt.text(6.1,215,'pyr',fontsize=tss)
        plt.plot([],[],' ',label='# Si-O coord. number   ')
    if fig9 or fig10:

        tss = 12
        if fig10:
            plt.text(5.45,135,'sft',fontsize=tss)
            plt.text(5,150,'liquid',fontsize=14)
            plt.plot(qtz_startliq_rho,qtz_startliq_P,'kx',label='PT liquidus intersection')

        if fig9:
            plt.text(5.36,120,'sft',fontsize=tss)
            plt.text(4.7,150,'liquid',fontsize=14)
            plt.plot(fs_startliq_rho,fs_startliq_P,'kx',label='PT liquidus intersection')
        plt.text(4.65,41,'stv',fontsize=tss)
        plt.text(4.93,71,'CaCl₂',fontsize=tss)
        #plt.plot([],[],' ',label='# Si-O coord. number')
    if fig13:
        tss = 12
        plt.text(5.37,120,'sft',fontsize=tss)
        plt.text(6.17,230,'pyrite',fontsize=tss)
        plt.text(5.,180,'liquid',fontsize=14)
        plt.text(4.6,37,'stv',fontsize=tss)
        plt.text(5.15,101,'CaCl$_2$',fontsize=tss)
        #plt.text(5.5,265,'melt curve model\npoorly constrained',fontsize=tss)
        #plt.annotate('',xy=(5.5,290), xytext=(5.8,290), arrowprops=dict(arrowstyle='<-'))
    if fig4:
        # Si-O coordination numbers
        # Usui 2020 SiO RCN
        #cncol='#E50000'
        cncol='black'
        fs=12
        # shifted because LDA DFT to match the liquid model from Karki LDA DFT
        plt.text(3.731+lda_liquid_density_shift,35,'Si-O CN\n3900K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(3.731+lda_liquid_density_shift,23,'5.1',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.199+lda_liquid_density_shift,50,'5.6',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.199+lda_liquid_density_shift,56,'5200K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.549+lda_liquid_density_shift,77,'5.6',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.549+lda_liquid_density_shift,83,'5400K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.942+lda_liquid_density_shift,108,'6.3',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        plt.text(4.942+lda_liquid_density_shift,114,'5800K',horizontalalignment='center',     verticalalignment='center',color=cncol,fontsize=fs)
        #plt.text(5.504,145,'6.3',horizontalalignment='center',     verticalalignment='center')

        #plt.plot(4.84,118,'ko') # Shaner 1981
        #plt.text(4.8,115,'7090K',horizontalalignment='right')

    if 0:
        #ind = np.where((driver.P > 108) & (driver.rho < 5.6))[0]
        ind = np.where((driver.P > 100) & (driver.rho < 5.6))[0]
        #plt.plot(driver.rho[ind],driver.P[ind],'c--',label='D10 PbO2 T=4000K')
        plt.plot(driver.rho[ind],driver.P[ind],'c--',label='D10 T=4000K')
        #plt.plot(driver.rho[ind],driver.P[ind]+17,'c--',label='')
        #plt.plot(driver.rho[ind],driver.P[ind]+34,'c--',label='')
        #plt.plot(driver.rho[ind],driver.P[ind]+10,'c--',label='D10 PbO2 4000K isotherm')
        plt.plot([5.1,5.15],[100,100],'c--')
        ind = np.where((driverstv.P < 105) & (driverstv.P >10))[0]
        plt.plot(driverstv.rho[ind],driverstv.P[ind],'c--')#,label='D10 Stv-CaCl2 T=4000K')
        #plt.plot(driverstv.rho[ind],driverstv.P[ind]-15.,'c--')
        #plt.plot(driverstv.rho[ind],driverstv.P[ind]+13.,'c--')
        #plt.plot(driverstv.rho[ind],driverstv.P[ind]+26,'c--')


    #nd = np.where((driverstv.P < 105) & (driverstv.P >10))[0]
    #plt.plot(driverstv.rho[ind],driverstv.P[ind],'c--')#,label='D10 Stv-CaCl2 T=4000K')


    if 0:
        plt.errorbar(eosdat_list[G16mc['id']].rho_l[1:4],eosdat_list[G16mc['id']].P[1:4],fmt='d--',mfc='white',color='#FFA500',label=G16mc['lab']) # pyrite structure
        plt.xlim(3.,6.5)
        plt.ylim(0,350)
        plt.text(6,250,r'pyrite',horizontalalignment='left',color='black')

    # plot on top
    #if fig5:
    #    ind = np.where(M70813.P > Pt2b)[0]
    #    plt.plot(M70813.rho[ind],M70813.P[ind],'g-',linewidth=3,label='Qtz Hugoniot')
    if fig6:
        ind = np.where( (qtz_hug_stv.P < ((Pt2a+Pt2b)/2)) )[0]
        plt.plot(qtz_hug_stv.rho[ind],qtz_hug_stv.P[ind],'g-',linewidth=3)
    #if fig2:
    #    ind = np.where(M70813.P > Pt2b)[0]
    #    plt.plot(M70813.rho[ind],M70813.P[ind],'r.',linewidth=3,label='Omega Hugoniot')


    #plt.plot(bosfsrho,bosfsp,'bd',label='B88 FS')
    

    # ------------------------------------------------------------------------
    #            0        1              2        3            4         5
    #           3100    3300           3500     3700         3900      4000
    isen_col = ['blue','darkturquoise','green','yellowgreen','orange','red']
    if fig6 or fig7:
#    if 1:
        lines='-'
        sval = 3500
        karkipvt.isenval_liquid(sval)
        plt.plot(karkipvt.rhoarr,karkipvt.isenval_p,lines,label='Liq S='+str(sval)+' J/K/kg',color=isen_col[2])
        sval = 3700
        karkipvt.isenval_liquid(sval)
        plt.plot(karkipvt.rhoarr,karkipvt.isenval_p,lines,label='Liq S='+str(sval),color=isen_col[3])
        sval = 3900
        karkipvt.isenval_liquid(sval)
        plt.plot(karkipvt.rhoarr,karkipvt.isenval_p,lines,label='Liq S='+str(sval),color=isen_col[4])
        sval = 4000
        karkipvt.isenval_liquid(sval)
        plt.plot(karkipvt.rhoarr,karkipvt.isenval_p,lines,label='Liq S='+str(sval),color=isen_col[5])
 
    if fig6 or fig7:
        lines=':'
        lw=2
        ind = np.where(stvpvt.rhoarr < 5.2)[0]
        sval=3100
        stvpvt.isenval_stv(sval)
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isenval_p[ind],lines,label='stv S='+str(sval),color=isen_col[0],linewidth=lw)
        sval=3300
        stvpvt.isenval_stv(sval)
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isenval_p[ind],lines,label='stv S='+str(sval),color=isen_col[1],linewidth=lw)
        sval=3500
        stvpvt.isenval_stv(sval)
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isenval_p[ind],lines,label='stv S='+str(sval),color=isen_col[2],linewidth=lw)
        sval=3700
        stvpvt.isenval_stv(sval)
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isenval_p[ind],lines,label='stv S='+str(sval),color=isen_col[3],linewidth=lw)
        sval=3900
        stvpvt.isenval_stv(sval)
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isenval_p[ind],lines,label='stv S='+str(sval),color=isen_col[4],linewidth=lw)
        sval=4000
        stvpvt.isenval_stv(sval)
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isenval_p[ind],lines,label='stv S='+str(sval),color=isen_col[5],linewidth=lw)

    if fig7:
        lines='--'
        ind = np.where(sftpvt.rhoarr > 5.0)[0]
        sval=3100
        sftpvt.isenval_sft(sval)
        plt.plot(sftpvt.rhoarr[ind],sftpvt.isenval_p[ind],lines,label='sft S='+str(sval),color=isen_col[0])
        sval=3300
        sftpvt.isenval_sft(sval)
        plt.plot(sftpvt.rhoarr[ind],sftpvt.isenval_p[ind],lines,label='sft S='+str(sval),color=isen_col[1])
        sval=3500
        sftpvt.isenval_sft(sval)
        plt.plot(sftpvt.rhoarr[ind],sftpvt.isenval_p[ind],lines,label='sft S='+str(sval),color=isen_col[2],linewidth=lw)
        sval=3700
        sftpvt.isenval_sft(sval)
        plt.plot(sftpvt.rhoarr[ind],sftpvt.isenval_p[ind],lines,label='sft S='+str(sval),color=isen_col[3],linewidth=lw)
        sval=3900
        sftpvt.isenval_sft(sval)
        plt.plot(sftpvt.rhoarr[ind],sftpvt.isenval_p[ind],lines,label='sft S='+str(sval),color=isen_col[4],linewidth=lw)
        sval=4000
        sftpvt.isenval_sft(sval)
        plt.plot(sftpvt.rhoarr[ind],sftpvt.isenval_p[ind],lines,label='sft S='+str(sval),color=isen_col[5],linewidth=lw)
        #plt.plot(sftpvt.rhoarr,sftpvt.sp5800,'r-',label='sft test')

    if 0:
        # low melt fit fs
        #[1.60478081 1.41001945]
        uparr = np.arange(4,6,.1)
        usarr = 1.41+1.60478*uparr
        parr = 2.2*uparr*usarr
        rarr = 1/((1/2.2)*(1-uparr/usarr))
        plt.plot(rarr,parr,'b--',label='L83 fit to low-P liquid')

        # low melt Marsh and L83 fs liquid fit =  [1.78054881 0.51652784]
        uparr = np.arange(4,6,.1)
        usarr = 0.5165+1.78*uparr
        parr = 2.2*uparr*usarr
        rarr = 1/((1/2.2)*(1-uparr/usarr))
        plt.plot(rarr,parr,'r--',label='M80+L83 fit to low-P liquid')

        # fit to low Z
        #[1.59150385 1.08393542]
        uparr = np.arange(4,6,.1)
        usarr = 1.084+1.5915*uparr
        parr = 2.2*uparr*usarr
        rarr = 1/((1/2.2)*(1-uparr/usarr))
        plt.plot(rarr,parr,'b-',label='fit to lowest Z liquid')

        # low fs
        #[1.58477234 1.21747131]
        uparr = np.arange(2.5,4,.1)
        usarr = 1.2174+1.58477*uparr
        parr = 2.2*uparr*usarr
        rarr = 1/((1/2.2)*(1-uparr/usarr))
        plt.plot(rarr,parr,'b-',label='W62 fit to FS solid')



    # ------------------------------------------------------------------------
    # ISOTHERMS
    if fig11 or fig14:
        #            0        1              2        3            4         5
        #           3000    4000           5000     5900         7000      xxxx
        isen_col = ['blue','darkturquoise','green','yellowgreen','orange','red']
        #plt.plot(karkipvt.rhoarr,karkipvt.p3000,':',label='K07 L T=3000K')
        #plt.plot(karkipvt.rhoarr,karkipvt.p4000,':',label='K07 L T=4000K')
        ind = np.where((karkipvt.rhoarr < 5.15))[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p3000[ind],':',label='K07 L T=3000K',color=isen_col[1])
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p4000[ind],':',label='K07 L T=4000K',color=isen_col[1])
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p5000[ind],':',label='K07 L T=5000K',color=isen_col[2])
#        plt.plot(karkipvt.rhoarr[ind],karkipvt.p5900[ind],':',label='K07 L T=5900K',color=isen_col[3])
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p6000[ind],':',label='K07 L T=6000K',color=isen_col[3])
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p7000[ind],':',label='K07 L T=7000K',color=isen_col[4])
        if 1:
            ind = np.where((stvpvt.rhoarr > 4.1) & (stvpvt.rhoarr < 5.3))[0]
            plt.plot(stvpvt.rhoarr[ind],stvpvt.p3000[ind],'--',label='stv T=3000K')
            plt.plot(stvpvt.rhoarr[ind],stvpvt.p4000[ind],'--',label='stv T=4000K',color=isen_col[1])
            plt.plot(stvpvt.rhoarr[ind],stvpvt.p5000[ind],'--',label='stv T=5000K',color=isen_col[2])
            plt.plot(stvpvt.rhoarr[ind],stvpvt.p6000[ind],'--',label='stv T=6000K',color=isen_col[3])
            plt.plot(stvpvt.rhoarr[ind],stvpvt.p7000[ind],'--',label='stv T=7000K',color=isen_col[4])
            #plt.plot(sftpvt.rhoarr,sftpvt.p3000,'--',label='sft T=3000K')
        ind = np.where((sftpvt.rhoarr > 5.35) & (sftpvt.rhoarr < 5.6))[0]
        #ind = np.where(sftpvt.p4000>50)[0]
        plt.plot(sftpvt.rhoarr[ind],sftpvt.p3000[ind],'-',label='sft T=3000K')
        plt.plot(sftpvt.rhoarr[ind],sftpvt.p4000[ind],'-',label='sft T=4000K',color=isen_col[1])
        #plt.plot(sftpvt.rhoarr[ind],sftpvt.p4800[ind],'-',label='sft T=4800K')
        plt.plot(sftpvt.rhoarr[ind],sftpvt.p5000[ind],'-',label='sft T=5000K',color=isen_col[2])
        plt.plot(sftpvt.rhoarr[ind],sftpvt.p6000[ind],'-',label='sft T=6000K',color=isen_col[3])
        plt.plot(sftpvt.rhoarr[ind],sftpvt.p7000[ind],'-',label='sft T=7000K',color=isen_col[4])

    if fig12:
        #            0        1              2        3            4         5
        #           3000    4000           5000     5900         7000      xxxx
        isen_col = ['blue','darkturquoise','green','yellowgreen','orange','red']
        #plt.plot(karkipvt.rhoarr,karkipvt.p3000,':',label='K07 L T=3000K')
        #plt.plot(karkipvt.rhoarr,karkipvt.p4000,':',label='K07 L T=4000K')
        ind = np.where((karkipvt.rhoarr < 5.65))[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p5000[ind],':',label='K07 L T=5000K',color=isen_col[2])
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p5900[ind],':',label='K07 L T=5900K',color=isen_col[3])
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p7000[ind],':',label='K07 L T=7000K',color=isen_col[4])
    if fig10:
        #            0        1              2        3            4         5
        #           3000    4000           5000     5900         7000      xxxx
        isen_col = ['blue','darkturquoise','green','yellowgreen','orange','red']
        #plt.plot(karkipvt.rhoarr,karkipvt.p3000,':',label='K07 L T=3000K')
        #plt.plot(karkipvt.rhoarr,karkipvt.p4000,':',label='K07 L T=4000K')
        #ind = np.where((karkipvt.rhoarr < 5.65))[0]
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p5000[ind],':',label='K07 Liquid T=5000K',linewidth=2,color='dodgerblue')
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p5900[ind],':',label='K07 L T=5900K',color=isen_col[3])
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p7000[ind],':',label='K07 L T=7000K',color=isen_col[4])
    if fig9:
        #            0        1              2        3            4         5
        #           3000    4000           5000     5900         7000      xxxx
        isen_col = ['blue','darkturquoise','green','yellowgreen','orange','red']
        #plt.plot(karkipvt.rhoarr,karkipvt.p3000,':',label='K07 L T=3000K')
        #plt.plot(karkipvt.rhoarr,karkipvt.p4000,':',label='K07 L T=4000K')
        #ind = np.where((karkipvt.rhoarr < 5.65))[0]
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p4500[ind],':',label='K07 Liquid T=4500K',linewidth=2,color='deepskyblue')
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p5900[ind],':',label='K07 L T=5900K',color=isen_col[3])
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p7000[ind],':',label='K07 L T=7000K',color=isen_col[4])
        #plt.plot(fs_hug_tra.rho,fs_hug_tra.P,'m-',label='transition')
    if 0: #if fig13:
        #            0        1              2        3            4         5
        #           3000    4000           5000     5900         7000      xxxx
        isen_col = ['blue','darkturquoise','green','yellowgreen','orange','red']
        #plt.plot(karkipvt.rhoarr,karkipvt.p3000,':',label='K07 L T=3000K')
        #plt.plot(karkipvt.rhoarr,karkipvt.p4000,':',label='K07 L T=4000K')
        ind = np.where((karkipvt.rhoarr < 5.65))[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.p5000[ind],':',label='K07 Liquid T=5000K',linewidth=2,color='dodgerblue')
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p5900[ind],':',label='K07 L T=5900K',color=isen_col[3])
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.p7000[ind],':',label='K07 L T=7000K',color=isen_col[4])


    #tval=3000
    #karkipvt.isotherm_liquid(tval)
    #plt.plot(karkipvt.rhoarr,karkipvt.isotherm_p,'k-')
    #plt.plot(qtz_hug_mix.rho,qtz_hug_mix.P,'-',color='orange')
    #plt.plot(qtz_hug_mix2.rho,qtz_hug_mix2.P,'-',color='orange')

    if fig6:
        # isotherms?
        # 1 amu/angstrom^3 = 1.6605391 g/cc
        # [v0, K0, Kp]
        # A^3, GPa, -
        #stv = [23.1559, 318.326, 4.3726]
        coesite = [546.80,97.4,4.3]
        amu = 60.0843 # SiO2 in amu
        # Vinet Coesite
        par = coesite
        rho0 = 16*amu*1.6605391/par[0]
        print('coesite = ',rho0)
        rhoarr = np.arange(rho0,rho0*2,.0001)
        eta = np.power(rho0/rhoarr,0.333)
        parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
        ind = np.where((parr <= 5))[0]
        plt.plot(rhoarr[ind],parr[ind],'-',color=colvinet,linewidth=3)
        fs=8
        plt.text(2.75,2,'coesite',color='black',fontsize=fs)
    if 0:
        #tval=4700
        #karkipvt.isotherm_liquid(tval)
        #ind = np.where((karkipvt.isotherm_p > 0) & (karkipvt.isotherm_p < 160))[0]
        #plt.plot(karkipvt.rhoarr[ind],karkipvt.isotherm_p[ind],'m:')
        tval=4000
        stvpvt.isotherm_stv(tval)
        ind = np.where((stvpvt.isotherm_p > 0) & (stvpvt.isotherm_p < 160))[0]
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isotherm_p[ind],'m:')
        tval=3000
        stvpvt.isotherm_stv(tval)
        ind = np.where((stvpvt.isotherm_p > 0) & (stvpvt.isotherm_p < 160))[0]
        plt.plot(stvpvt.rhoarr[ind],stvpvt.isotherm_p[ind],'m:')
    
        plt.plot(stvdat.rho,stvdat.P,'b--')

    if 0:
        # check Robinson Pyrite 5000 K data points
        # first is 5000 K isotherm
        vrshift=0.125 # to match Hugoniot with PBE functional
        pyr_rho=np.asarray([5.1936,5.4465,5.7446,6.0105,6.3665])+vrshift
        pyr_P=np.asarray([111.578,141.018,180.411,226.62,286.588])
        # single disordered pyrite from liquid at 5000 K; seems to dense for the qtz Hugoniot
        plt.plot(pyr_rho,pyr_P,'x-',label='Pyrite')
        plt.plot(5.0195+vrshift,102.51,'D',label='d-pyr')

    #if fig9 or fig10:
    #    plt.plot(eosdat_list[G24mclrho['id']].rho_liquid,eosdat_list[G24mclrho['id']].P_liquid,'-',color='red',label=G24mclrho['lab']+' (two-phase)')
    # ------------------------------------------------------------------------
    plt.xlabel(r'Density (g/cm$^3$)')
    plt.ylabel('Pressure (GPa)')
    #get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()
    #specify order of items in legend
    order = range(len(handles))

    plt.xlim(3.5,6.)
    plt.ylim(0,225)
    figname = 'B0-SiO2-P-rho-phase-diagram'
    if not(fig9 or fig10):
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=10)
    if fig1:
        # plt.title(r'SiO$_2$ Phase Diagram')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=9)
        figname = 'B1-SiO2-P-rho-phase-diagram-Fig1b'
        plt.text(3.1,230,'B',fontsize=16)
        plt.xlim(3.,6.3)
        plt.ylim(0,250)
        plt.ylabel('\nPressure (GPa)')
    if fig2:
        plt.title(r'SiO$_2$ Phase Diagram with Fused Silica and Quartz Hugoniots')
        figname = 'B2-SiO2-P-rho-phase-diagram-both-omega'
    if fig3:
        plt.title(r'SiO$_2$ Phase Diagram with Fused Silica and Quartz Hugoniots')
        figname = 'B3-SiO2-P-rho-phase-diagram-both-Rusbank'
    if fig4:
        plt.title(r'SiO$_2$ Phase Diagram with Fluid Coordination Numbers')
        figname = 'B4-SiO2-P-rho-phase-diagram-CNs'
    if fig5:
        #plt.title(r'SiO$_2$ Phase Diagram with Fused Silica and Quartz Hugoniots')
        figname = 'B5-SiO2-P-rho-phase-diagram-both-clean-no-amodeo'
    if fig6:
        #plt.title(r'SiO$_2$ Phase Diagram with Isentropes')
        figname = 'B6-SiO2-P-rho-phase-diagram-with-isentropes'
        plt.xlim(2.5,5.5)
        plt.ylim(0,130)
    if fig7:
        #plt.title(r'SiO$_2$ Phase Diagram with Isentropes')
        figname = 'B7-SiO2-P-rho-phase-diagram-with-isentropes-zoom'
        plt.xlim(4.4,5.6)
        plt.ylim(50,160)
    if fig8:
        #plt.title(r'SiO$_2$ Phase Diagram with Isentropes')
        figname = 'B8-SiO2-P-rho-phase-diagram-sft-trans'      
        plt.title('New SiO2 Hugoniot Models (8/30/2023)')

    if fig9:
        plt.plot([],[],' ',label=' ')#'spacer')
        #plt.title(r'SiO$_2$ Phase Diagram')
        figname = 'B9-SiO2-P-rho-phase-diagram-FS-Fig8A'
        #plt.text(4.25,225,'A',fontsize=16)
        plt.xlim(4.2,5.4)
        plt.ylim(0,350)
        handles, labels = plt.gca().get_legend_handles_labels()
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        order = range(len(handles))
        print(len(order))
#        order = [0,1,2,3,6,4,5,7,8,9,10,11,13,12,14]
#        order = [0,1,2,3,7,4,5,6,8,9,10,11,12,14,13,15]
        order = [0,1,2,3,6,5,7,8,9,10,11,12,13,14,15,4]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],ncol=2,loc='upper left',fontsize=10,title='A. Fused silica')
    if fig10:
        plt.plot([],[],' ',label=' ')#'spacer')
        plt.plot([],[],' ',label=' ')#'spacer')
        #plt.title(r'SiO$_2$ Phase Diagram')
        figname = 'B10-SiO2-P-rho-phase-diagram-qtz-Fig8B'
        #plt.text(4.45,225,'B',fontsize=16)
        #plt.plot(cst_hug_new.rho,cst_hug_new.P,'k-')
        plt.xlim(4.4,5.6)
        plt.ylim(0,300)
        handles, labels = plt.gca().get_legend_handles_labels()
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        order = range(len(handles))
        print(len(order))
#        order = [0,1,2,3,6,4,5,8,7,10,9]
#        order = [0,1,2,3,7,4,5,6,8,11,10]
        order = [0,1,2,3,5,4,6,10,9,12,11,7,8]

        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='upper left',fontsize=10,title='B. α-Quartz z-cut',ncol=2)
    if fig11:
        plt.title(r'SiO$_2$ Phase Diagram')
        figname = 'B11-SiO2-P-rho-phase-diagram-check'
        plt.xlim(3.1,5.6)
    if fig12:
        plt.title(r'SiO$_2$ Phase Diagram')
        figname = 'B12-SiO2-P-rho-phase-diagram-check'
        plt.xlim(3.1,5.6)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=10)
        #plt.xlim(4.4,6.5)
        plt.ylim(0,250)
        #plt.legend(loc='upper left',fontsize=12,title='B. α-Quartz z-cut')
    if fig13:
        plt.plot([],[],' ',label=' ')#'spacer')
        plt.plot([],[],' ',label=' ')#'spacer')
        plt.plot([],[],' ',label=' ')#'spacer')
        
        
        figname = 'B13-SiO2-P-rho-phase-diagram-cst-Fig8C'      
        plt.xlim(4.4,6.3)
        plt.ylim(0,370)
        handles, labels = plt.gca().get_legend_handles_labels()
        for i in range(len(labels)):
            print(i,labels[i],handles[i])
        #order = range(len(handles))
        print(len(order))
        order = [0,1,2,3,5,4,6,10,11,12,7,8,9]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='upper left',fontsize=10,title='C. Coesite',ncol=2)
    if fig14:
        plt.title(r'SiO$_2$ Phase Diagram')
        figname = 'B14-SiO2-P-rho-phase-diagram-isotherms'
        plt.xlim(3.1,5.6)        
    if fig15:
        plt.title(r'SiO$_2$ Phase Diagram')
        figname = 'B15-SiO2-P-rho-phase-diagram-offset'
        plt.xlim(3.1,5.6)        
    if 0:
        # check 0 K, 300 K transitions
        plt.plot([0,10],[Pt2_0k,Pt2_0k])
        plt.plot([0,10],[Pt3_0k,Pt3_0k])
        plt.plot([0,10],[Pt4_0k,Pt4_0k])
    if 0:
        print(Pt2,Pt3,Pt4)
        # check melt curve transitions
        #plt.plot([0,10],[Pt2,Pt2])
        plt.plot([0,10],[fs_p_startliq,fs_p_startliq])
        plt.plot([0,10],[qtz_p_startliq,qtz_p_startliq])
        plt.plot([0,10],[fs_wiggle_P1,fs_wiggle_P1],'--')
        plt.plot([0,10],[fs_wiggle_P2,fs_wiggle_P2],'--')
        plt.plot([0,10],[qtz_wiggle_P1,qtz_wiggle_P1],':')
        plt.plot([0,10],[qtz_wiggle_P2,qtz_wiggle_P2],':')
    #plt.plot([0,10],[130,130])
    #plt.plot([0,10],[140,140])
    #plt.plot([0,10],[Pt3,Pt3])
    #plt.plot([0,10],[75.5,75.5])
    #plt.plot([0,10],[70,70])
    #plt.plot([0,10],[65.9,65.9])
    #plt.plot([0,10],[84,84])
    #plt.plot([4.6,4.6],[0,250])
    #plt.plot([5.11,5.11],[0,250])
    # plt.plot([0,10],[120,120]) # remember where cst enters the melt curve

    plt.tight_layout()
    plt.savefig(outputfigpath+figname+'.pdf')
    if fig1 or fig9 or fig10 or fig13:
        if not(fig12):
            plt.savefig(figurefigpath+figname+'.pdf')
    print('Saved figure: '+figname+'.pdf')
    plt.show()

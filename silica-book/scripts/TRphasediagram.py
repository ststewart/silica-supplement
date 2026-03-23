# TRphasediagram.py
def plotomegaTR(hug,label=False,branch=1):
    omegapmin=85
    omegapmax=111
    usmax=9.5
    mss = 5
    col9='olive'

    if not((branch==4) or (branch==5) or (branch==2)):
    
        zzz = np.where((hug.P > omegapmin) & (hug.P < omegapmax))[0]
        if label:
            plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss,label='CaCl$_2$-type branch (Omega)')
        else:
            plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss)#,label='CaCl$_2$-type')

    if (branch==1) or (branch==2):
        col9='black'
        usmin = 9.23
        usmax = 9.5
        zzz = np.where((hug.us > usmin) & (hug.us < usmax))[0]
        if label:
            plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss,label='Shock front transition T↓ (Omega)')
        else:
            plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss)#,label='Transition')
    if (branch==1) or (branch==4) or (branch==5):
        if not(branch==5):
            #mss = 2
            col9='springgreen'
            #zzz = np.where((hug.us >= usmax) & (hug.rho <=5.2))[0]
            zzz = np.where((hug.us >= usmax) & (hug.T <= 6200))[0]
            if label:
                plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss,label='Cool liquid shock front T↑ (Omega)')
            else:
                plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss)#,label='Disordered phases')
        if not(branch==4):
            rhotrans=5.2
            #mss = 2
            col9='green'
            #zzz = np.where(hug.rho > rhotrans)[0]
            zzz = np.where(hug.T > 6200)[0]
            if label:
                plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss,label='Liquid (Omega)')
            else:
                plt.plot(hug.rho[zzz],hug.T[zzz],'.',color=col9,markersize=mss)#,label='Liquid')

#pind = 1
#for pind in [1,2]:

if 1:
    fig1=False
    fig2=False
    fig3=False
    fig4=False
    fig5=False
    fig6=False
    fig7=False
    if pind ==1:
        fig1=True # density-temperature qtz Hugoniot interpretation
    if pind ==2:
        fig1=True # density-temperature qtz Hugoniot interpretation
    if pind ==3:
        fig1=True # density-temperature qtz Hugoniot interpretation
    if pind == 4:
        fig2=True
    if pind == 5:
        fig3=True # overview plot
    if pind == 6:
        fig4=True # solid to metastable melting metastable surface
    if pind == 7:
        fig7=True

        
    widthsol2=2
    plt.figure(figsize=(10, 6))
    plt.rcParams.update({'font.size': 14})

    if 0:
        # Connelly vapor curves
        plt.plot(C16bsrho,C16bstemp,'k-',label='C16 HSC-vWRK vapor curve')
        plt.plot(C16bdrho,C16bdtemp,'k--',label='C16 HSC-CSRK vapor curve')
        #plt.plot(C16rsrho,C16rstemp,'r-',label='C16 SiO(l)-vWRK')
        plt.plot(1.036,6074,'ko')
        plt.plot(.555,5548,'ko')
        #plt.plot(.782,5675,'ro')
        plt.plot(kdome.rho_l,kdome.T,'-',color='grey',label='K12 MANEOS vapor curve')
        plt.plot(np.min(kdome.rho_l),np.max(kdome.T),'o',color='grey')
        plt.plot(kdome.rho_v,kdome.T,'-',color='grey')

    if 1:
        # Karki DFT Isentropes
        ind = np.where(karkipvt.rhoarr < 4.2)[0]
        iii = np.where(karkipvt.stempsarr == 3500)[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.st3500[ind],':',label='Liquid S='+str(np.int_(np.round(karkipvt.svalsarr[iii])))+' J/K/kg')
        ind = np.where(karkipvt.rhoarr < 4.6)[0]
        iii = np.where(karkipvt.stempsarr == 4500)[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.st4500[ind],':',label='S='+str(np.int_(np.round(karkipvt.svalsarr[iii]))))
        iii = np.where(karkipvt.stempsarr == 5000)[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.st5000[ind],':',label='S='+str(np.int_(np.round(karkipvt.svalsarr[iii]))))
        ind = np.where(karkipvt.rhoarr < 5)[0]
        iii = np.where(karkipvt.stempsarr == 5500)[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.st5500[ind],':',label='S='+str(np.int_(np.round(karkipvt.svalsarr[iii]))))
        ind = np.where(karkipvt.rhoarr < 5)[0]
        iii = np.where(karkipvt.stempsarr == 7000)[0]
        plt.plot(karkipvt.rhoarr[ind],karkipvt.st7000[ind],':',label='S='+str(np.int_(np.round(karkipvt.svalsarr[iii]))))
        print('liquid=',karkipvt.svalsarr,karkipvt.stempsarr)
        #liquid= [2873.13445363 3149.4759091  3388.85376607 3600.         3788.87674241 3959.73653905 4115.71931245 4259.20955435 4392.06076792] 
        #    [3000 3500 4000 4500 5000 5500 6000 6500 7000]
        #solid =  [2952.63563485 3186.20612626 3388.53408893 3567.         3726.64319745 3871.05803596 4002.89845408 4124.179882   4236.46894549] 
        #    [3000 3500 4000 4500 5000 5500 6000 6500 7000]

        # Solid stv Isentropes
        print('solid = ',stvpvt.svalsarr,stvpvt.stempsarr)
        rhomin=2.5
        ind = np.where(stvpvt.rhoarr > rhomin)[0]
        iii = np.where(stvpvt.stempsarr == 3500)[0]
        plt.plot(stvpvt.rhoarr[ind],stvpvt.st3500[ind],'--',label='Solid S='+str(np.int_(np.round(stvpvt.svalsarr[iii])))+' J/K/kg')
        ind = np.where(stvpvt.rhoarr > rhomin)[0]
        iii = np.where(stvpvt.stempsarr == 4500)[0]
        plt.plot(stvpvt.rhoarr[ind],stvpvt.st4500[ind],'--',label='S='+str(np.int_(np.round(stvpvt.svalsarr[3]))))
        ind = np.where(stvpvt.rhoarr > rhomin)[0]
        iii = np.where(stvpvt.stempsarr == 5000)[0]
        plt.plot(stvpvt.rhoarr[ind],stvpvt.st5000[ind],'--',label='S='+str(np.int_(np.round(stvpvt.svalsarr[4]))))
        ind = np.where(stvpvt.rhoarr > rhomin)[0]
        iii = np.where(stvpvt.stempsarr == 5500)[0]
        plt.plot(stvpvt.rhoarr[ind],stvpvt.st5500[ind],'--',label='S='+str(np.int_(np.round(stvpvt.svalsarr[5]))))

    # sts melt curve
    #ind = np.where(stsmcsmooth.P < Pt3a)[0]
    #ind = np.where(stsmcsmooth.P < Pt4)[0]
    #plt.plot(stsmcsmooth.rho[ind],stsmcsmooth.T[ind],'m-',label='liquidus',linewidth=3)
    if 1:
        plt.plot(stsmcsmooth.rho,stsmcsmooth.T,'m-',linewidth=3,label='proposed liquidus')
    #ind = np.where((stsmcsmooth.P >= Pt3a) & (stsmcsmooth.P <200))[0]
    #plt.plot(stsmcsmooth.rho[ind],stsmcsmooth.T[ind],'m:',linewidth=3)
    #sind = np.where(stsmcsmooth.P < 65)[0]
    if not(pind==11):
        plt.plot(stsmcsmooth.rho_s,stsmcsmooth.T,'-',label='proposed solidus',color='grey',linewidth=3)
        sind = np.where((stsmcsmooth.P < 110) & (stsmcsmooth.P > 0))[0]
        plt.plot(stsmcsmooth.rho_s[sind],stsmcsmooth.T[sind],'-',color='grey',linewidth=3)
            #plt.plot(stsmcsmooth.S_l,stsmcsmooth.T-240,'m-',label='liquidus')

    # find the intersections in rho-T for the 4500K reference isentrope
    #plt.plot([1.97,3.23],[3100,3425],'kP',label='dT=325 K')
    #plt.plot([1.8,2.2],[3100,3100],':')
    #scross = np.interp(3100,np.flip(kdome.T),np.flip(kdome.ent_l))
    #print('S at the 4500 K ref isentrope = ',scross)

    # Coesite sliver
    cvol1 = 6.85 # cm3/mol
    cvol2 = 7.15 # cm3/mol
    crho1 = sio2mol/cvol1/3
    crho2 = sio2mol/cvol2/3
    plt.fill_between([crho1,crho2],[2000,2000],[3000,3000],color='grey',alpha=0.3)
    qrho2 = sio2mol/8.15/3
    plt.fill_between([qrho2,2.65],[2000,2000],[2723,2723],color='grey',alpha=0.3)
    #plt.text(2.45,2010,'qtz')
    #plt.text(2.8,2010,'cst')
    #plt.text(4.25,2010,'stv')
    #print('coesite = ',crho1/3,crho2/3,sio2mol,qrho2,crho1,crho2)

    if not(pind==8 or pind==10):
        if not(fig4 or fig6):
            plt.text(4.35,2100,'stv')
            plt.text(4.,5300,'liquid')
            plt.text(3.85,4000,'liq+stv')
        else:
            plt.text(4.4,6000,'liquid')
            plt.text(4.4,2100,'stv')
            
        plt.text(4.87,2100,'CaCl$_2$')
        plt.text(5.35,2100,'sft')
    if (pind==8) or (pind==10):
        plt.text(4.4,6000,'liquid')
        plt.text(5.0,4000,'metastable liquid\n shock front')
       

    if fig3 or fig6 or fig7:
        if not(pind==9):
            ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'z')[0]
            plt.errorbar(shockdat_list[B88RRqtz['id']].rho[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='yellow',\
                 yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label=B88RRqtz['lab']+'(z)')
            #plt.plot([shot26.rho],[shot26.T],'gD',mfc='yellow',label='Gas gun qtz(z) new') # MG release

    if 0:
        plt.errorbar(shockdat_list[B88Rfs['id']].rho,shockdat_list[B88Rfs['id']].T,fmt='bo',mfc='lightblue',yerr=shockdat_list[B88Rfs['id']].T_std,label=B88Rfs['lab'])
        ind = np.where(shockdat_list[B88RRqtz['id']].orientation == 'z')[0]
        plt.errorbar(shockdat_list[B88RRqtz['id']].rho[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='yellow',\
             yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label=B88RRqtz['lab']+'(z)')
#plt.errorbar(shockdat_list[B88RRqtz['id']].rho,shockdat_list[B88RRqtz['id']].T,fmt='go',mfc='red',yerr=shockdat_list[B88RRqtz['id']].T_std,label=B88RRqtz['lab'])

    if pind==10:
        ind = np.where((shockdat_list[B88RRqtz['id']].orientation == 'pz') & (shockdat_list[B88RRqtz['id']].P >120))[0]
        plt.errorbar(shockdat_list[B88RRqtz['id']].rho[ind],shockdat_list[B88RRqtz['id']].T[ind],fmt='go',mfc='orange',\
             yerr=shockdat_list[B88RRqtz['id']].T_std[ind],label=B88RRqtz['lab']+'(perp. z)')
        

    if pind<4:
        #plt.errorbar(shockdat_list[B88qtz['id']].rho,shockdat_list[B88qtz['id']].T,fmt='go',mfc='white',yerr=shockdat_list[B88qtz['id']].T_std,label=B88qtz['lab'])
#        plt.errorbar(fs_hug_Hicks06.rho,fs_hug_Hicks06.T,fmt='b+',\
#                     yerr=fs_hug_Hicks06.T_std,label='H06 FS')
        ind = np.where(fs_hug_liq.P > 80)[0]
        plt.plot(fs_hug_liq.rho[ind],fs_hug_liq.T[ind],'b-',label='FS Hug Liquid NEW')
        #plt.errorbar(shockdat_list[L04qtz['id']].rho,shockdat_list[L04qtz['id']].T,fmt='gv',mfc='white',xerr=shockdat_list[L04qtz['id']].rho_std,yerr=shockdat_list[L04qtz['id']].T_std,label=L04qtz['lab']+'(x)')


    if 0:
        # draw solid phase boundaries
        # CaCl2 to sft
        plt.plot(p23bnd.rho_a,p23bnd.T,'c-',label='CaCl2-sft boundary')
        plt.plot(p23bnd.rho_b,p23bnd.T,'c-',)
        # stv to CaCl2
        plt.plot(p12bnd.rho_a,p12bnd.T,'c--',label='stv-CaCl2 boundary')
        #plt.plot(p12bnd.rho_b,p12bnd.T,'c--')
   # if pind==6 or pind==11:
   ##     # CaCl2 to sft
    #    plt.plot(p23bnd.rho_a,p23bnd.T,'c-',label='CaCl2-sft boundary')
    #    plt.plot(p23bnd.rho_b,p23bnd.T,'c-',)
        

    if pind<4:
        # plot model Hugoniots from new fits
        plt.plot(fs_hug_stv.rho,fs_hug_stv.T,'b--',label='FS Hug solid NEW')
        #ind = np.where((qtz_hug_solid.P > 111) & (qtz_hug_sft.P < 140))[0]
        #plt.plot(qtz_hug_sft.rho[ind],qtz_hug_sft.T[ind],'r--',linewidth=3,label='Qtz Hug Sft NEW')
        ind = np.where(qtz_hug_solid.P < 111)[0]
        plt.plot(qtz_hug_solid.rho[ind],qtz_hug_solid.T[ind],'g--',label='Qtz Hug solid NEW')
#        ind = np.where(qtz_hug_stv.P < 70)[0]
#        plt.plot(qtz_hug_stv.rho[ind],qtz_hug_stv.T[ind],'g--')
    if 1:
        ind = np.where(qtz_hug_solid.P < qtz_wiggle_P1)[0]
        plt.plot(qtz_hug_solid.rho[ind],qtz_hug_solid.T[ind],'g-',label='Qtz Hugoniot solid')
        ind = np.where(qtz_hug_liq.P > qtz_wiggle_P2)[0]
        plt.plot(qtz_hug_liq.rho[ind],qtz_hug_liq.T[ind],'g--',label='Qtz Hugoniot liquid')
    if 1:
        ind = np.where(fs_hug_stv.P < fs_wiggle_P1)[0]
        plt.plot(fs_hug_stv.rho[ind],fs_hug_stv.T[ind],'b-',label='FS Hugoniot solid')
        ind = np.where(fs_hug_liq.P > fs_wiggle_P2)[0]
        plt.plot(fs_hug_liq.rho[ind],fs_hug_liq.T[ind],'b--',label='FS Hugoniot liquid')
    if pind==9:
        ind2 = np.where((qtz_hug_liq.us > qtz_usmax) & (qtz_hug_liq.us <qtz_us_startliq))[0]
        # interpolate for the liquidus temperature
        #t_liq = np.interp(qtz_hug_liq.P[ind2],stsmcsmooth.P,stsmcsmooth.T)
        #plt.plot(qtz_hug_liq.rho[ind2],t_liq,'-',color='springgreen',linewidth=3,label='Cool liquid shock front T↑')
        ind = np.where((stsmcsmooth.P > Pt3) & (stsmcsmooth.P <= qtz_p_startliq))[0]
        plt.plot(stsmcsmooth.rho[ind],stsmcsmooth.T[ind],'-',color='springgreen',linewidth=3,label='Cool liquid shock front T↑')
        plt.text(4.9,6100,'~4700 K',fontsize=10)
        plt.annotate('~5500 K',[5.05,6125],fontsize=10)

    if 0:
        if fig3:
            plotomegaTR(M70811,label=1)
            plotomegaTR(M70813)
            #plotomegaTR(M70817)
            #plotomegaTR(M70820)
        if fig4 and not(fig5): # plot solid and amorphization branches
            plotomegaTR(M70811,label=1,branch=2)
            plotomegaTR(M70813,branch=2)
            #plt.plot([np.interp(Pt2,stsmcsmooth.P,stsmcsmooth.rho),5.15],[4900,4900],'m--',label='Metastable stv liquidus')
        if fig5: # plot solid and amorphization branches
            plotomegaTR(M70811,label=1,branch=3)
            plotomegaTR(M70813,branch=3)
        if fig6 and (pind==8 or pind==10): # plot glass heating and liquid branches
            plotomegaTR(M70811,label=1,branch=4)
            plotomegaTR(M70813,branch=4)
            #plt.plot([np.interp(Pt2,stsmcsmooth.P,stsmcsmooth.rho),5.15],[4900,4900],'m--',label='Metastable stv liquidus')
        if fig7: # plot glass heating and liquid branches
            plotomegaTR(M70811,label=1,branch=4)
            plotomegaTR(M70813,branch=4)
            #plt.plot([np.interp(Pt2,stsmcsmooth.P,stsmcsmooth.rho),5.15],[4900,4900],'m--',label='Metastable stv liquidus')
            # calculate metastable liquidus
            
        if pind==9: # plot glass heating and liquid branches
            plotomegaTR(M70811,label=1,branch=5)
            plotomegaTR(M70813,branch=5)
  
    if pind==11:

        ind2 = np.where(qtz_hug_glass.T <=6200)[0]
        ind2 = np.where((qtz_hug_glass.us >9.23) & (qtz_hug_glass.us <10.6))[0]
        # find the liquid densities that correspond to the P and T of the Hugoniot data
        rnew = []
        tnew = []
        pnew = []
        for iii in range(len(ind2)):
            karkipvt.isotherm_liquid(qtz_hug_glass.T[ind2[iii]])
            rnew.append(np.interp(qtz_hug_glass.P[ind2[iii]],karkipvt.isotherm_p,karkipvt.rhoarr))
            #pnew.append(np.interp(qtz_hug_glass.P[ind2[iii]],karkipvt.isotherm_p,karkipvt.rhoarr))
            tnew.append(qtz_hug_glass.T[ind2[iii]])


    plt.plot(shockdat_list[A23fs['id']].rho,shockdat_list[A23fs['id']].T,'bD',mfc='white')
    plt.plot(shockdat_list[A23qtz['id']].rho,shockdat_list[A23qtz['id']].T,'gD',mfc='white')
   
    
    plt.xlabel(r'Density (g/cm$^3$)')
    plt.ylabel('Temperature (K)')
    #plt.title('Model Vapor Curves & DFT Isentropes')
    plt.ylim(2000,10000)
    #plt.xlim(0,5.5)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    #plt.xlim(2,6.)
    #plt.legend()

    if fig1:
        figname='C1-SiO2-T-rho-phase-diagram-hugs'
        plt.ylim(2000,9000)
        plt.xlim(2,5.5)    

    if pind==2:
        figname='C1b-SiO2-T-rho-phase-diagram-hugs'
        plt.ylim(2000,20000)
        plt.xlim(2,6)    
    if pind==3:
        figname='C3-SiO2-T-rho-phase-diagram-hugs-new'
        plt.title('New SiO2 Hugoniot Models (8/31/2023)')
        plt.ylim(2000,9000)
        plt.xlim(2,5.5)      
    if fig2:
        figname='C4-Phase-check'
        plt.ylim(2000,7000)
        plt.xlim(2,5.5)       
    if fig3:
        figname='C5-TR-'
        plt.ylim(2000,7000)
        plt.xlim(2,5.5)
    if fig4 and not(fig5):
        figname='C6-TR-'
        plt.ylim(2000,7000)
        plt.xlim(2.,5.5)
        plt.title('Shock front transition from xtl to melt emission')
        #plt.text(4.2,6700,'shock front states with decreasing crystallization (black)',fontsize=10)
    if fig5:
        figname='C7-TR-'
        plt.ylim(2000,7000)
        plt.xlim(2.,5.5)
        plt.title('Shock front transition from xtl to melt emission')
        #plt.text(4.2,6700,'equilibrium states behind partially liquid shock front (black)',fontsize=10)
    if fig6 and (pind==8):
        figname='C8-TR-'
        plt.ylim(2000,7000)
        plt.xlim(2.,5.5)
        plt.title('Hugoniot branch with cool liquid shock front')
        plt.text(4.2,6700,'cool liquid states in shock front (lime)',fontsize=10)
    if fig6 and (pind==9):
        figname='C9-TR-'
        plt.ylim(2000,7000)
        plt.xlim(2.,5.5)
        plt.title('Hugoniot branch with cool liquid shock front')
        plt.text(4.2,6700,'shock states behind cool liquid shock front (lime)',fontsize=10)
    if fig6 and (pind==10):
        figname='C10-TR-'
        plt.ylim(2000,7000)
        plt.xlim(2.,5.5)
        plt.title('Orientation effects')
        #plt.text(4.3,6700,'cool liquid states in shock front (lime)',fontsize=10)
    if fig7:
        figname='C7-TR-'
        plt.ylim(2000,8000)
        plt.xlim(2.25,5.5)
        #plt.title('Hugoniot branch with cool liquid shock front')
        #plt.text(4.2,6700,'equilibrium states behind partially opaque shock front',fontsize=10)
    plt.tight_layout()
    plt.savefig(outputfigpath+figname+'.pdf')
    print('Saved figure: '+figname+'.pdf')
    plt.show()

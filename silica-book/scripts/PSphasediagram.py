# PSphasediagram.py

#pind = 1
#for pind in [1,2]:

if 1:
    fig1=False
    fig2=False
    if pind ==1:
        fig1=True 

# PS-plot

plt.figure(figsize=(10, 6))
plt.rcParams.update({'font.size': 14})

# plot liquidus and solidus
plt.plot(stsmcsmooth.S_l,stsmcsmooth.P,'m-',label='liquidus')
ind = np.where(stsmcsmooth.P <= Pt4)[0]
plt.plot(stsmcsmooth.S_s[ind],stsmcsmooth.P[ind],'-',color='grey',label='solidus')

if not(fig1):
    # ALPHA QUARTZ HUGONIOT BRANCHES
    # ALPHA QUARTZ HUGONIOT BRANCHES
    #ind = np.where(qtz_hug_liq.P>qtz_p_startliq)[0]
    ind = np.where(qtz_hug_liq.P>146)[0]
    plt.plot(qtz_hug_liq.S[ind],qtz_hug_liq.P[ind],'g-',linewidth=1,label='Qtz Hugoniot liquid')
    ind = np.where(qtz_hug_sft.P > 113)[0]
    plt.plot(qtz_hug_stv.S,qtz_hug_stv.P,'g--',label='Qtz Hugoniot solid')
    # CaCl2 region
    ind = np.where((qtz_hug_cacl2.P < 112)& (qtz_hug_cacl2.P > 63))[0]
    plt.plot(qtz_hug_cacl2.S[ind],qtz_hug_cacl2.P[ind],'g--')#,label='qtz Hugoniot-CaCl2 region NEW')
    # Disorder region
    ind = np.where((qtz_hug_glass.P > qtz_wiggle_P1) & (qtz_hug_glass.P < qtz_p_startliq))[0]
    #ind = np.where((qtz_hug_sft.P > 113) & (qtz_hug_sft.P < 145))[0]
    plt.plot(qtz_hug_glass.S[ind],qtz_hug_glass.P[ind],'-',linewidth=2,color='lime',label='Qtz Hugoniot metastable')


if 0:
    # Plot original Kraus et al. 2012 Hugoniot
    kqtzparr = np.arange(100,500) # GPa
    kqtzsarr = 1793+22.4*kqtzparr-0.0426*np.power(kqtzparr,2)+4.319e-05*np.power(kqtzparr,3)-1.742e-8*np.power(kqtzparr,4)
    kqtzserrarr = 127.3+0.328*kqtzparr-0.000193*np.power(kqtzparr,2)
    ind = np.where(kqtzparr>110)[0]
    plt.plot(kqtzsarr[ind],kqtzparr[ind],'g:',linewidth=1,label='K12 qtz Hugoniot-liq region')

if not(fig1):
    # FUSED SILICA HUGONIOT BRANCHES
    ind = np.where(fs_hug_liq.P > 78)[0]
    plt.plot(fs_hug_liq.S[ind],fs_hug_liq.P[ind],'b-',linewidth=1,label='FS Hugoniot liquid')
    plt.plot(fs_hug_mix1.S,fs_hug_mix1.P,'-',color='cyan',linewidth=2,label='FS Hugoniot metastable')
    ind = np.where((fs_hug_stv.P < 66) & (fs_hug_stv.P > 30))[0]
    plt.plot(fs_hug_stv.S[ind],fs_hug_stv.P[ind],'b--',label='FS Hugoniot solid')
    #plt.plot(fs_hug_mix2.S,fs_hug_mix2.P,'r-',linewidth=4,label='Calculated mixed FS hug NEW')

if 0:
    # only 70 to 110 GPa
    kfsparrx = np.arange(70,110)
    kfssarrx = 2284+26.82*kfsparrx-0.0723*kfsparrx**2
    plt.plot(kfssarrx,kfsparrx,'b:',label='K12 fs Hugoniot-liq region')
    
    # plot tie points
    plt.plot([3567,3820,3820],[70,70,113],'X',color='cyan',label='K12 S tie pts')
    plt.plot(karkipvt.refpt[3],karkipvt.refpt[1],'kP',label='Liquid S ref pt')
    plt.plot(stvpvt.refpt[3],stvpvt.refpt[1],'kX',label='Stv S ref pt')
    plt.plot(sftpvt.refpt[3],sftpvt.refpt[1],'kx',label='Sft S ref pt')

if not(fig1):
    # denote the pressure range for temperature wiggles in omega data
    plt.plot([0,10000],[qtz_wiggle_P1,qtz_wiggle_P1],':',color='grey',label='Temp inversions')
    plt.plot([0,10000],[qtz_wiggle_P2,qtz_wiggle_P2],':',color='grey')
    plt.plot([0,10000],[fs_wiggle_P1,fs_wiggle_P1],':',color='grey')
    plt.plot([0,10000],[fs_wiggle_P2,fs_wiggle_P2],':',color='grey')

plt.xlabel('Sp. Entropy (J/K/kg)')
plt.ylabel('Pressure (GPa)')

xloc=2600
plt.text(xloc,130,r'seifertite')
plt.text(xloc,120,r'$\alpha$-PbO$_2$-type')
plt.text(xloc,90,'CaCl$_2$-type')
plt.text(xloc,25,'stishovite')
plt.text(4200,55,'liquid')
plt.xlim([2500,4900])    
plt.ylim([0,Pt4])    

#plt.ylim([100,165])    


if 0:
    tval=6000
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 30) & (stvpvt.isotherm_p < 115))[0]
    plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'c:')
    
    tval=5900
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 30) & (stvpvt.isotherm_p < 115))[0]
    plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'c:')
    
    tval=5000
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 30) & (stvpvt.isotherm_p < 115))[0]
    plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'c:')
    tval=4500
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 30) & (stvpvt.isotherm_p < 115))[0]
    plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'c:')
    tval=4000
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 30) & (stvpvt.isotherm_p < 115))[0]
    plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'c:')
    
    tval=6000
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 100) & (sftpvt.isotherm_p < 160))[0]
    plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'c:')
    tval=5900
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 100) & (sftpvt.isotherm_p < 160))[0]
    plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'c:')
    
    tval=5000
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 100) & (sftpvt.isotherm_p < 160))[0]
    plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'c:')
    tval=4500
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 100) & (sftpvt.isotherm_p < 160))[0]
    plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'c:')
    tval=4000
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 100) & (sftpvt.isotherm_p < 160))[0]
    plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'c:')
    
    
    tval=4500
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 0) & (karkipvt.isotherm_p < 160))[0]
    plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m:')
    tval=5000
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 0) & (karkipvt.isotherm_p < 160))[0]
    plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m:')
    tval=5900
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 0) & (karkipvt.isotherm_p < 160))[0]
    plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m:')
    tval=6000
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 0) & (karkipvt.isotherm_p < 160))[0]
    plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m:')

if 0:
    tval=4900
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 10) & (karkipvt.isotherm_p < 160))[0]
    #plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m--',label='4900 K liq')
    p_start = 112
    s_start = np.interp(p_start,karkipvt.isotherm_p[ind],karkipvt.isotherm_s[ind])
    print(p_start, s_start)
    
    
    tval=6150
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 10) & (karkipvt.isotherm_p < 160))[0]
    #plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m--',label='6150 K liq')
    p_end = 146
    s_end = np.interp(p_end,karkipvt.isotherm_p[ind],karkipvt.isotherm_s[ind])
    print(p_end, s_end)
    plt.plot([s_start,s_end],[p_start,p_end],'m-.',label='Metastable stv liquidus')
if 1:    
    tval=4000
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 80) & (sftpvt.isotherm_p < 160))[0]
    plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'c:',label='4000 K sft')
    tval=5000
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 80) & (sftpvt.isotherm_p < 160))[0]
    plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'k:',label='5000 K sft')
    tval=5900
    sftpvt.isotherm_sft(tval)
    ind = np.where((sftpvt.isotherm_p > 80) & (sftpvt.isotherm_p < 160))[0]
    #plt.plot(sftpvt.isotherm_s[ind],sftpvt.isotherm_p[ind],'c--',label='6150 K sft')
    
    tval=4000
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 60) & (stvpvt.isotherm_p < 120))[0]
    plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'c--',label='4000 K CaCl2')
    tval=5000
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 60) & (stvpvt.isotherm_p < 120))[0]
    plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'k--',label='5000 K CaCl2')
    tval=5900
    stvpvt.isotherm_stv(tval)
    ind = np.where((stvpvt.isotherm_p > 60) & (stvpvt.isotherm_p < 120))[0]
    #plt.plot(stvpvt.isotherm_s[ind],stvpvt.isotherm_p[ind],'--',label='6150 K stv')
    
if 0:
    tval=6100
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > Pt2) & (karkipvt.isotherm_p < 140))[0]
    plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m--',label='T=6100 K liq')
    # Karki model not working well in this density-temperature range; need to look at this range some more
    tval=6300
    karkipvt.isotherm_liquid(tval)
    ind = np.where((karkipvt.isotherm_p > 100) & (karkipvt.isotherm_p < 160))[0]
    plt.plot(karkipvt.isotherm_s[ind],karkipvt.isotherm_p[ind],'m-.',label='T=6300 K liq')

#plt.plot(Strans+sshift,Ptrans,'o')
#plt.text(Strans+20+sshift,Ptrans,str(np.round(Ttrans))+'K',size=10)
#plt.plot(Strans,Ptrans,'o')
#plt.text(np.interp(113,qtz_hug_glass.P,qtz_hug_glass.S)+20,113,'4900 K',size=10)
#plt.plot(np.interp(113,qtz_hug_glass.P,qtz_hug_glass.S),113,'o')
#plt.text(np.interp(113,qtz_hug_glass.P,qtz_hug_glass.S)+20,prel,str(np.round(trel))+'K',size=10)

#trel = np.interp(113,qtz_hug_liq.P,qtz_hug_liq.T)
#prel = np.interp(np.interp(113,qtz_hug_glass.P,qtz_hug_glass.S),qtz_hug_liq.S,qtz_hug_liq.P)
#print(prel)
#plt.text(np.interp(113,qtz_hug_glass.P,qtz_hug_glass.S)+20,113,'4900 K',size=10)
#plt.plot(np.interp(113,qtz_hug_glass.P,qtz_hug_glass.S),113,'o')
#plt.plot(np.interp(113,qtz_hug_glass.P,qtz_hug_glass.S),prel,'o')

plt.plot(p23bnd.S_a,p23bnd.P,'-',linewidth=2,color='cadetblue',label='CaCl2-Sft')
plt.plot(p23bnd.S_b,p23bnd.P,'-',linewidth=2,color='cadetblue')
plt.plot(p12bnd.S_a,p12bnd.P,'--',linewidth=2,color='cadetblue',label='Stv-CaCl2')

#plt.plot([0,5000],[140,140])
#plt.plot([0,5000],[150,150])
#plt.plot([3820,3820],[0,200])
#plt.plot([3600,3600],[0,200])
#plt.plot([4000,4000],[0,200])

if not(fig1):
    plt.title('Disorder Region on Quartz Hugoniot')


plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=10)
figname='Disorder-P-S-plot'
plt.tight_layout()
plt.savefig(outputfigpath+figname+'.pdf')
plt.show()

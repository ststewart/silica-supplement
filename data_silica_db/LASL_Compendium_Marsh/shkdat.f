      program shkdat
c-----------------------------------------------------------------------
c
c                              g. i. kerley
c               division 1533, sandia national laboratories
c                            19 december 1984
c
c  a program to make a file of hugoniot data for plotting.  data are
c  either entered at the tty or read from a file named by the user.  to
c  enter data from the tty, enter tty when prompted.  then enter values
c  of initial density rzro, shock velocity us, and particle velocity up,
c  terminating input by entering zero for rzro.  to read data from a
c  file, first enter the file name when prompted.  then, when prompted
c  for the material description, enter the first 30 characters of the
c  material name as it is given in la-7887-m.  the program will make an
c  output file with a name entered by the user.  for each material, two
c  comment cards (! in column 1) will be printed. the first card will
c  contain the the word start, along with a label of up to 65 characters
c  entered by the user.  the second card will give headers for the col-
c  umns of data that follow, namely the particle velocity up, shock vel-
c  ocity us, density rho, and pressure p.
c
c  the following files will be of interest to users.
c
c    ud:[gkerley.user]hugla.dat   -  los alamos compendium
c    ud:[gkerley.user]hugl3.dat   -  lawrence livermore compendium
c    ud:[gkerley.user]index.dat   -  index of data on hugla and hugl3
c    ud:[gkerley.user]shkdat.exe  -  this code (executable file)
c    ud:[gkerley.user]shkdat.for  -  this code (fortran source file)
c
c-----------------------------------------------------------------------
      character fname*30,fsv1*30,fsv2*30,matid*30,key*30,labl*70,last*18
      dimension up(500),us(500),rho(500),p(500),rzro(500)
      data last,fsv1,fsv2/'              9999',' ',' '/
 1    continue
      write(6,2)
 2    format(38h enter name of input file, tty, or end)
      read(5,3)fname
 3    format(a)
      if(fname.eq.'end') stop
      if(fname.eq.'tty') go to 30
      write(6,4)
 4    format(34h enter material description or end)
      read(5,3)matid
      if(matid.eq.'end') stop
      if(fname.ne.fsv1) then
        if(fsv1.ne.' ') close(unit=1)
        open(unit=1,file=fname,status='old',readonly)
        fsv1 = fname
      else
        rewind 1
        endif
 5    read(1,3)key
      if(key.eq.last) go to 20
      if(key.ne.matid) go to 5
      read(1,3)key
      read(1,3)key
      write(6,3)key
      read(1,3)key
      write(6,3)key
      n = 1
      rsum = 0.
 10   read(1,11)iend,rzro(n),us(n),up(n)
 11   format(14x,i2,3x,e10.0,10x,2e6.0)
cdebug      write(6,12)rzro(n),us(n),up(n)
cdebug 12   format(1x,3(1pe15.5))
      if(iend.eq.99) go to 50
      rho(n) = rzro(n)/(1.-up(n)/us(n))
      p(n) = rzro(n)*us(n)*up(n)
      rsum = rsum+rzro(n)
      n = n+1
      go to 10
 20   write(6,21)
 21   format(24h material is not on file)
      go to 1
 30   write(6,31)
 31   format(47h enter rzro, us, and up, terminating with zeros)
      n = 1
      rsum = 0.
 40   read(5,*)rzro(n),us(n),up(n)
      if(rzro(n).eq.0) go to 50
      rho(n) = rzro(n)/(1.-up(n)/us(n))
      p(n) = rzro(n)*us(n)*up(n)
      rsum = rsum+rzro(n)
      n = n+1
      go to 40
 50   n = n-1
      if(n.eq.0) go to 60
      fn = n
      rsum = rsum/fn
      write(6,51)rsum
 51   format(26h average initial density =,1pe12.5)
      write(6,55)
 55   format(26h enter name of output file)
      read(5,3)fname
      if(fname.ne.fsv2) then
        if(fsv2.ne.' ') close(unit=2)
        open(unit=2,file=fname,status='new')
        fsv2 = fname
        endif
      write(6,52)
 52   format(26h enter label for plot file)
      read(5,3)labl
      if(labl.eq.'end') stop
      write(2,53)labl
 53   format(3h!  ,a/4h!   ,8hup(km/s),5x,8hus(km/s),4x,
     1  9hrho(g/cc),5x,6hp(gpa))
      write(2,54)(up(i),us(i),rho(i),p(i),i=1,n)
 54   format(4(1pe13.5))
c  replace the two preceding cards with the following, in order
c  to print out rzro, us, up, pzro, and ezro for recentering.
c$$$      pzro = 0.
c$$$      ezro = 0.
c$$$      write(2,54)(rzro(i),us(i),up(i),pzro,ezro,i=1,n)
c$$$ 54   format(5(1pe13.5))
      go to 1
 60   write(6,61)
 61   format(33h error - no points were generated)
      go to 1
      end

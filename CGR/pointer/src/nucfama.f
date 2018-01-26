c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::::::: 
      subroutine fmcx()
      implicit double precision (o-z, a-h)
      character a*1, icc*1, d*1, fc*1
      character ibb*6, icdc*6
      character ifmt*108
      dimension d(45), fc(12)
      common /convr1/ in(25), inf(60), isel(12), kid(3, 12), ntrn, itrn(
     &7, 60), itr(60), cc(500), jp, inv(60), ip, iq, itot
      common /convr2/ ifmt(2), ibb(100), icdc(500), icc(501)
      common /inter1/ index, icol, i1, i2, nchar
      common /inter2/ a(800)
c     NIN = COUNTS NUMBER OF INPUT FIELDS (S,R,F,A)                     
c         
c     NP  = COUNTS NUMBER OF INPUT FIELDS (F,A)                         
c         
c     IFMT(108,1) = STORES A-FORMAT FOR INPUT                           
c         
c     IFMT(108,2) = STORES F-FORMAT FOR INPUT                           
c         
c     NSID = COUNTS NUMBER OF SELECTION/REJECTION FIELDS                
c         
c     ISEL = FIELD THAT EACH SELECTION/REJECTION FIELD IS IN.           
c         
c                                                                       
c         
      data d / '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', ' ', 
     &',', '.', '-', '(', ')', '=', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 
     &'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
     &'U', 'V', 'W', 'X', 'Y', 'Z', '*', ' ' /
      data fc / 'F', 'A', 2*' ', 'S', 'R', 'X', ',', ')', 3*' ' /
      nsid = 0
      nin = 0
      jp = 0
      np = 0
      do 25 i = 1, 2
      do 25 j = 1, nchar
c                                                                       
c         
c     INTERPRET FIRST SET OF PARENTHESES -  RECORD LENGTH IN WORDS      
c         
c   W A R N I N G -- FM HAS WORDS BUT IN(12) STORES NO. OF CHARACTERS IN
c REC. -- 
c                                                                       
c         
   25 ifmt(i)(j:j) = a(j)
      call parn
      if (index) 40, 340, 40
   40 in(2) = xnon(icol)
      in(2) = in(2)
c                                                                       
c         
c      INTERPRET 2ND SET OF (  )  -  INPUT FORMAT                       
c         
c     STORE COLUMN WITH ( IN L.                                         
c         
c                                                                       
c         
      icol = i2
      call parn
      if (index) 95, 340, 95
   95 l = i1
      goto 105
c                                                                       
c         
c     DETERMINE REPETITION FACTOR                                       
c         
c                                                                       
c         
  100 call nspacen
  105 irep = xno2n(icol)
      if (irep) 115, 110, 115
c                                                                       
c         
c     DETERMINE WHAT FIELD IT IS                                        
c         
c                                                                       
c         
  110 irep = 1
  115 do 120 i = 1, 9
      if (a(icol) .eq. fc(i)) goto 125
  120 continue
c-----I = 1  F FIELD                                                    
c         
c     I = 2  A FIELD                                                    
c         
c     I = 3  NOT USED                                                   
c         
c     I = 4  NOT IN USE                                                 
c         
c     I = 5  S (SELECTION)                                              
c         
c     I = 6  R (REJECTION)                                              
c         
c     I = 7  X (SKIP)                                                   
c         
c     I = 8  COMMA                                                      
c         
c     I = 9  RIGHT RARENTHESIS                                          
c         
      goto 340
c                                                                       
c         
c     X - SPACING                                                       
c         
c                                                                       
c         
  125 goto (200, 200, 340, 340, 280, 280, 150, 100, 290), i
  150 call nspacen
  155 if (a(icol) .ne. d(12)) goto 105
c                                                                       
c         
c     A, F FIELDS                                                       
c         
c                                                                       
c         
      goto 100
  200 i3 = icol
      j = xno1n(icol)
      if (j) 205, 340, 205
  205 do 215 j = 1, irep
      nin = nin + 1
      np = np + 1
c                                                                       
c         
c     JP COUNTS NO. OF FLOATING PT. FIELDS.  INV HAS A SPECIAL PURPOSE  
c         
c     IN OVERLAY 2.                                                     
c         
c                                                                       
c         
      if (i .eq. 2) goto 215
  210 jp = jp + 1
      inv(jp) = nin
      inf(jp) = np
  215 continue
      if (i .eq. 2) then
      do 216 j = i3, icol - 2
  216 ifmt(2)(j:j) = ifmt(2)(j + 1:j + 1)
      ifmt(2)(icol - 1:icol - 1) = 'X'
      if (irep .gt. 1) then
      write(unit=62, fmt='(X,A45)') 
     &'NO REPETITION COUNT ALLOWED FOR A FIELDS'
      goto 340
      end if
      goto 155
c                                                                       
c         
c     BLANK OUT DECIMAL IN F - FIELDS FOR A - FORMAT                    
c         
c                                                                       
c         
      end if
  220 icol = i3
      ifmt(1)(icol:icol) = d(18)
  225 icol = icol + 1
      if (a(icol) .ne. d(13)) goto 225
  230 ifmt(1)(icol:icol) = d(11)
      icol = icol + 1
      if ((a(icol) .eq. d(12)) .or. (a(icol) .eq. d(16))) goto 115
c                                                                       
c         
c     SELECTION AND REJECTION FIELDS                                    
c         
c                                                                       
c         
      goto 230
  280 nsid = nsid + 1
      in(i + 3) = in(i + 3) + 1
      isel(nsid) = nin + 1
      j = 0
      if (a(icol + 1) .eq. d(1)) j = 1
      if (a(icol + 1) .eq. d(2)) j = 2
      if (a(icol + 1) .eq. d(3)) j = 3
      if (a(icol + 1) .eq. d(4)) j = 4
      if (a(icol + 1) .eq. d(5)) j = 5
      if (a(icol + 1) .eq. d(6)) j = 6
  277 if (j .eq. 0) goto 340
  278 nin = nin + 1
      ifmt(1)(icol:icol) = d(18)
      ifmt(2)(icol:icol) = d(j)
      ifmt(1)(icol + 1:icol + 1) = d(j)
      ifmt(2)(icol + 1:icol + 1) = d(41)
      icol = icol + 2
c                                                                       
c         
c     STORE NO. OF FIELDS AND NO. OF VARIABLE FIELDS.                   
c         
c                                                                       
c         
      goto 155
  290 in(17) = nin
c                                                                       
c         
c     PACK A AND F FORMATS TO LEFT                                      
c         
c                                    
      in(19) = np
  300 ifmt(1)(1:1) = '('
      ifmt(2)(1:1) = '('
      do 333 j = 1, 2
      k = 2
      do 320 i = l, icol
      if (ifmt(j)(i:i) .eq. d(11)) goto 320
      ifmt(j)(k:k) = ifmt(j)(i:i)
      k = k + 1
  320 continue
      do 325 i = k, nchar
  325 ifmt(j)(i:i) = d(11)
  333 continue
      index = 1
      return 
  340 index = 0
      end
      subroutine nspacen()
      implicit double precision (o-z, a-h)
      character a*1
      common /inter1/ index, icol, i1, i2, nchar
      common /inter2/ a(800)
c                                                                       
c         
c     FINDS NEXT NON-BLANK COLUMN                                       
c         
c                                                                       
c         
      icol = icol + 1
      entry nsp1n()
    5 index = 1
    6 if (a(icol) .eq. ' ') then
      icol = icol + 1
      if (icol .gt. nchar) goto 55
      else
      goto 7
      end if
      goto 6
    7 return 
c                                                                       
c         
c     CHECKS FOR RIGHT AND LEFT PARENTHESIS STARTING WITH ICOL          
c         
c     I1 ON RETURN = COL. AFTER (   ;   I2 = COL. BEFORE )              
c         
c     ICOL ON RETURN = 1ST NON-BLANK COLUMN AFTER (                     
c         
c                                                                       
c         
      entry parn()
    8 if (a(icol) .ne. '(') then
      icol = icol + 1
      if (icol .gt. nchar) goto 55
      else
      goto 9
      end if
      goto 8
    9 icol = icol + 1
      i1 = icol
   10 if (a(icol) .ne. ')') then
      icol = icol + 1
      if ((icol .gt. nchar) .or. (a(icol) .eq. '(')) goto 55
      else
      goto 11
      end if
      goto 10
   11 i2 = icol - 1
      icol = i1
      goto 5
   55 index = 0
      end
      
c     program nuca77
      subroutine nucfama(datfile,jobfile,profile,terfile)
      implicit double precision (o-z, a-h)
      character title*12
      character ctime*24
      character filename*45, datfile*30, jobfile*30, profile*30,
     &   terfile*30
      integer idebug
      dimension title(7)     
      common /pets/ cat(12), dog, iguana(27), idebug
      common /connet/ lindex
c     call getarg(1,datfile)
c     call getarg(2,jobfile)
c     call getarg(3,profile)
c     call getarg(4,terfile)
      open(unit=48,file='INTERMED.1',status='new',
     & form='unformatted') 
      open(unit=52,file='INTERMED.ALP',status='new',
     & form='FORMATTED') 
      open(unit=49,file=datfile,status='old',form='FORMATTED') 
      open(unit=60,file=jobfile,status='old',form='FORMATTED') 
      open(unit=61,file=profile,status='new',form='FORMATTED') 
      open(unit=62,file=terfile,status='new',form='FORMATTED') 
      read(unit=60, fmt=1) title
      if (title(1)(1:5) .eq. 'DEBUG') then
      idebug = 1
      else
      backspace(unit=60, err=20051)
      idebug = 0
      end if
20051 continue
      write(unit=52, fmt=5) 0, idebug
    5 format(2i6)
      read(unit=60, fmt=1) title
      if (title(1)(1:5) .eq. 'TITLE') then
      write(unit=52, fmt=4) title
    4 format(a12)
      write(unit=61, fmt=2) title
      write(unit=62, fmt=2) title
      else
      backspace(unit=60, err=20052)
20052 continue 
      title(1) = 'TITLE=?     '
      do 3 i = 2, 7
    3 title(i) = '            '
      write(unit=61, fmt=2) title
      write(unit=62, fmt=2) title
      write(unit=52, fmt=4) title
      end if
    1 format(7a12)
    2 format(//6h FILE ,7a12)
      call fdate(ctime)
      write(unit=62,fmt=1361) ctime
      write(unit=61,fmt=1362) ctime
 1361 format(//1x,a24,9x,42hNUCFAMA77 CREATED 08-01-84.  TERSE OUTPUT.)
 1362 format(//1x,a24,9x,43hNUCFAMA77 CREATED 08-01-84.  PROLIX OUTPUT.)
      write(unit=62, fmt=1363) 
      write(unit=61, fmt=1363) 
 1363 format(//1x,130(1h-))
      inquire(unit=49, name=filename) 
      write(unit=61, fmt='(''0DATA FILE NAME = '',A45)') filename
      write(unit=62, fmt='(''0DATA FILE NAME = '',A45)') filename
      inquire(unit=60, name=filename) 
      write(unit=61, fmt='(''0JOB FILE NAME = '',A45)') filename
      write(unit=62, fmt='(''0JOB FILE NAME = '',A45)') filename
      inquire(unit=61, name=filename) 
      write(unit=61, fmt='(''0PROLIX FILE NAME = '',A45)') filename
      write(unit=62, fmt='(''0PROLIX FILE NAME = '',A45)') filename
      inquire(unit=62, name=filename) 
      write(unit=61, fmt='(''0TERSE FILE NAME = '',A45)') filename
      write(unit=62, fmt='(''0TERSE FILE NAME = '',A45)') filename
      lindex = 0
      iguana(1) = 86
      call over1
      if (iguana(1) .ne. 2) goto 100
      call over2
  100 continue
      if (lindex .eq. 2) write(unit=*, fmt=*) 
     &'ERROR IN NUCFAMA PROCESSING INPUT'
      write(unit=62, fmt=1363) 
      write(unit=61, fmt=1363) 
      endfile(unit=48) 
      endfile(unit=52) 
      close(unit=48) 
      close(unit=52) 
      close(unit=61) 
      close(unit=62) 
      return
      end
      subroutine over1()
      implicit double precision (o-z, a-h)
cPOINTER ONLY                    
      character mpc*1, xxc*1
cPOINTER ONLY                    
      character iddy*2
      character a*1, mc*1, icc*1
      character ibb*6, icdc*6, ib*6
      character ncode*62
      character ifmt*108
cPOINTER ONLY                    
      dimension mpc(3), xxc(21)
      dimension mc(24)
      dimension ncode(28), ib(11)
      common /pets/ cat(12), dog, iguana(27), idebug
      common /convr1/ in(25), inf(60), isel(12), kid(3, 12), ntrn, itrn(
     &7, 60), itr(60), cc(500), jp, inv(60), ip, iq, itot
      common /convr2/ ifmt(2), ibb(100), icdc(500), icc(501)
      common /inter1/ index, icol, i1, i2, nchar
      common /inter2/ a(800)
      common /poin/ yi(10), zi(10), pre(7), id(7), jaf(3), jsex(2), 
     &jpro(4), npp, iccc, jli
      common /oddsn/ pie(10), fi(10), zeey, npie, izy, nfi
      common /connet/ lindex
cPOINTER                         
      data ncode / 'C(1)*X(I)', 'C(1)*X(I)*X(J)', 'C(1)/X(I)', 
     &'C(1)*X(J)/X(I)', 'X(I)**C(1)', 'C(1)**(C(2)*X(I))', 'EXP(X(I))', 
     &'X(I)+C(1)*X(J)', 'C(1)+C(2)*X(I)', 'LN(X(I))', 'LOG10(X(I))', 
     &'SIN(C(1)*X(I))', 'COS(C(1)*X(I))', 'ARCSIN(C(1)*X(I))', 6*' ', 
     &'C(1) IF X(I) <= A(1), C(2) IF X(I) > A(1)', 
     &'C(1) IF X(I)<A(1), C(2) IF A(1)<=X(I)<A(2), C(3) IF X(I)>=A(2)', 
     &'C(1) IF X(I) = ANY A(K), C(2) OTHERWISE', 
     &'C(1) IF X(I) >< ANY A(K), C(2) OTHERWISE', 
     &'C(1) IF X(I) = ANY A(K), X(I) OTHERWISE', 
     &'C(1) IF X(I) >< ANY A(K), X(I) OTHERWISE', 'C(M) IF X(I) = A(M)'
     &, 'EQUAL INTERVALS FROM A(1) TO A(2) STEP A(3)' /
      data mc / 'C', 'F', 'T', 2*'S', 'A', 'P', 'S', 'P', 'L', 'P', 'C'
     &, 'M', 'R', 'I', 'D', 'F', 'O', 'X', 'R', 'I', 'T', 2*' ' /
      data mpc / 'A', 'N', 'U' /
      data xxc / 'L', 2*'P', 'I', 'A', 'Q', 'P', 'G', 'S', 'G', 'I', 'T'
     &, 'O', 'D', 'F', 'U', 'R', 'N', 'X', 'R', ' ' /
      jndex = -7
      write(unit=62, fmt=4) 
    4 format(//1x,130(1h=))
      write(unit=61, fmt=4) 
      ntrn = 0
      do 5 i = 1, 25
    5 in(i) = 0
      do 6 i = 1, 60
      itr(i) = 0
      inv(i) = 0
      inf(i) = 0
      do 6 j = 1, 7
    6 itrn(j,i) = 0
      do 7 i = 1, 12
      isel(i) = 0
      do 7 j = 1, 3
    7 kid(j,i) = 0
      do 8 i = 1, 100
cSELECTION CRITERIA STORED SEQUENTIA
    8 ibb(i) = ' '
      do 9 i = 1, 2
    9 ifmt(i) = ' '
      do 10 i = 1, 500
      cc(i) = 0.
      icdc(i) = ' '
   10 icc(i) = ' '
      npie = 0
      iccc = 0
      nfi = 0
      in(15) = 1
      in(16) = 2
      in(21) = 4
      in(22) = 5
      in(23) = 6
      jaf(1) = 1
      jaf(2) = 0
      jaf(3) = 9
      jsex(1) = 1
      jsex(2) = 2
      jli = 0
      niaz = 0
      ip = 0
      iq = 0
c---- READ CONTROL CARD                                                 
c         
c---- UP TO 9 CONTINUATION CARDS ARE ALLOWED                            
c         
c---- AN * IN LAST NON-BLANK COLUMN IS A CONTINUATION                   
c         
      itot = 0
   20 i1 = 1
      nchar = 80
   25 read(unit=60, fmt=26, end=148) (a(i),i = i1, nchar)
   26 format(80a1)
 8690 if (nchar .ge. 800) goto 40
      do 404 i = nchar, i1, -1
      if (a(i) .eq. ' ') goto 404
      if (a(i) .eq. '*') goto 35
      goto 40
  404 continue
      goto 55
   35 a(i) = ' '
      i1 = nchar + 1
      nchar = nchar + 80
c---- PRINT CONTROL CARD. GO TO FIRST NON-BLANK COLUMN ON CARD          
c         
      goto 25
   40 icol = 1
      call nsp1n
      i1 = icol
      icol = icol + 1
      if ((a(i1) .eq. 'P') .and. (a(icol) .eq. 'A')) goto 149
      write(unit=62, fmt=41) (a(i),i = 1, nchar)
      write(unit=61, fmt=41) (a(i),i = 1, nchar)
c---- IDENTIFY CONTROL CARD                                             
c         
c---- I = 1 CC  (INDICATES END OF JOB)                                  
c         
c----   = 2 FM                                                          
c         
c----   = 3 TR                                                          
c         
c----   = 4 SI                                                          
c         
c----   = 5 SD                                                          
c         
c----   = 6 AF  (NOT NEEDED)                                            
c         
c----   = 7 PO  (NOT NEEDED)                                            
c         
c----   = 8 SX  (NOT NEEDED)                                            
c         
c----   = 9 PR  (PROBAND CONTROL CARD, NOT IN USE YET)                  
c         
c----   = 10 LI                                                         
c         
c----   = 11 PT (FIRST CONTROL CARD)                                    
c         
c----   = 12 LZ                                                         
c         
c----   = 13 PI                                                         
c         
c----   = 14 FI                                                         
c         
c----   = 15 ZY                                                         
c         
c                                                                       
c         
   41 format(1x,80a1)
      do 50 i = 1, 11
      if ((a(i1) .ne. mc(i)) .or. (a(icol) .ne. mc(i + 11))) goto 50
      goto (60, 120, 130, 140, 140, 730, 3010, 740, 20, 759, 770), i
   50 continue
      if ((a(i1) .eq. 'L') .and. (a(icol) .eq. 'Z')) goto 758
      if ((a(i1) .eq. 'P') .and. (a(icol) .eq. 'I')) goto 70
      if ((a(i1) .eq. 'F') .and. (a(icol) .eq. 'I')) goto 80
      if ((a(i1) .eq. 'Z') .and. (a(icol) .eq. 'Y')) goto 880
   55 write(unit=62, fmt=56) i1, i2, icol, (a(i),i = 1, nchar)
      write(unit=61, fmt=56) i1, i2, icol, (a(i),i = 1, nchar)
   56 format(//1x,18hERROR  I1,I2,ICOL ,3i4/(6x,80a1))
   57 read(unit=60, fmt=26, end=475) (a(i),i = 1, 80)
      write(unit=62, fmt=41) (a(i),i = 1, 80)
      write(unit=61, fmt=41) (a(i),i = 1, 80)
c                                                                       
c         
c     CC    END OF JOB CARD                                             
c         
c                                                                       
c         
      goto 57
c                                                                       
c         
c     AF(A=1,N=0,U=2)                                                   
c         
c                                                                       
c         
   60 goto 148
  730 call parn
      if (index) 751, 55, 751
  751 do 752 i = 1, 3
      if (a(icol) .eq. mpc(i)) goto 753
  752 continue
      goto 55
  753 call nspacen
      jaf(i) = xno1n(icol)
      if (a(icol) .eq. ')') goto 20
      call nspacen
c                                                                       
c         
c     SX(M=1,F=2)                                                       
c         
c                                                                       
c         
      goto 751
  740 call parn
      if (index) 741, 55, 741
  741 if (a(icol) .ne. 'M') goto 742
      inter = 1
      goto 745
  742 if (a(icol) .ne. 'F') goto 55
      inter = 2
  745 call nspacen
      jsex(inter) = xno1n(icol)
      if (a(icol) .eq. ')') goto 20
      call nspacen
c                                                                       
c         
c     LI(  ,  ,  ,  ,  ,  ,ETC.)                                        
c         
c     LZ(  ,  ,  ,  ,  ,  ,ETC.)                                        
c         
c                                                                       
c         
      goto 741
  758 jli = jli + 2
      do 760 i = 1, 10
  760 zi(i) = 0.
      niaz = 0
      call parn
      if (index) 765, 55, 765
  765 niaz = niaz + 1
      if (niaz .gt. 10) goto 55
      zi(niaz) = xnon(icol)
      if (a(icol) .eq. ')') goto 20
      call nspacen
      goto 765
  759 jli = 1
      do 761 i = 1, 10
  761 yi(i) = 0.
      nia = 0
      call parn
      if (index) 762, 55, 762
  762 nia = nia + 1
      if (nia .gt. 10) goto 55
      yi(nia) = xnon(icol)
      if (a(icol) .eq. ')') goto 763
      call nspacen
      goto 762
  763 in(24) = nia
c                                                                       
c         
c     FI(  ,  ,  ,  ,  ,ETC.)                                           
c         
c                                                                       
c         
      goto 20
   80 do 81 i = 1, 10
   81 fi(i) = 0.
      nfi = 0
      call parn
      if (index) 82, 55, 82
   82 nfi = nfi + 1
      if (nfi .gt. 10) goto 55
      fi(nfi) = xnon(icol)
      if (a(icol) .eq. ')') goto 20
      call nspacen
c                                                                       
c         
c     PO(FA=1,MO=2,FP=4,MP=5,CP=6)                                      
c         
c                                                                       
c         
      goto 82
 3010 call parn
      if (index) 3011, 55, 3011
 3011 j = 0
      iddy = ' '
 3012 j = j + 1
      iddy(j:j) = a(icol)
      call nspacen
      if (a(icol) .ne. '=') goto 3012
 3013 if (j - 2) 3014, 3014, 55
 3014 call nspacen
      inter = xnon(icol)
      if (iddy .ne. 'FA') goto 3015
      in(15) = inter
      goto 3018
 3015 if (iddy .ne. 'MO') goto 3016
      in(16) = inter
      goto 3018
 3016 if (iddy .ne. 'FP') goto 3017
      in(21) = inter
      goto 3018
 3017 if (iddy .ne. 'MP') goto 3009
      in(22) = inter
      goto 3018
 3009 if (iddy .ne. 'CP') goto 55
      in(23) = inter
 3018 if (a(icol) .ne. ',') goto 55
 3019 call nspacen
      goto 3011
 3008 if (a(icol) .eq. ')') goto 20
c                                                                       
c         
c FMCX,  TRCX,  SIDCX                                                   
c         
c                                                                       
c         
      goto 55
  120 call fmcx
  125 if (index) 20, 55, 20
  130 call trcx
      goto 125
  140 call sidcx
c                                                                       
c         
c     PI(  ,  ,  ,  ,  ,  ETC.)                                         
c         
c                                                                       
c         
      goto 125
   70 npie = 0
      call parn
      if (index) 71, 55, 71
   71 npie = npie + 1
      if (npie .gt. 10) goto 55
      pie(npie) = xnon(icol)
      if (a(icol) .eq. ')') goto 20
      call nspacen
c                                                                       
c         
c     ZY(VALUE)                                                         
c         
c                                                                       
c         
      goto 71
  880 izy = 1
      call parn
      if (index) 881, 55, 881
  881 zeey = xnon(icol)
c                                                                       
c         
c     PT(LI= ,PT= ,PO= ,ID= ,AF= ,QU= ,PR= ,GN= ,SX= ,GR= ) (J OR C)    
c         
c                                                       2ND SET OF () OP
cTIONAL   
c                                                                       
c         
      goto 20
  770 jndex = 0
      izy = 0
      call parn
      iorder = 0
      if (index) 771, 55, 771
  771 do 772 i = 1, 10
      if ((a(icol) .eq. xxc(i)) .and. (a(icol + 1) .eq. xxc(i + 10))) 
     &goto 7700
  772 continue
      goto 55
 7700 icol = icol + 1
      call nspacen
      iorder = xno1n(icol)
  773 if (i - 2) 774, 775, 776
  774 in(1) = iorder
      goto 780
  775 in(4) = iorder
      goto 780
  776 if (i .eq. 3) goto 777
      if (i .eq. 9) goto 778
      if (i .eq. 10) goto 779
      in(i + 6) = iorder
      goto 780
  777 in(5) = iorder
      goto 780
  778 in(20) = iorder
      goto 780
  779 in(6) = iorder
  780 if (a(icol) .eq. ')') goto 781
      call nspacen
      goto 771
  781 call parn
      if (index) 782, 20, 782
  782 if (a(icol) .eq. 'J') iccc = -1
c---- PRINTOUT FOR INPUT AND OUTPUT                                     
c         
c                                                                       
c         
c---- INPUT SPECIFICATIONS:     IN(2) WORDS PER RECORD                  
c         
c                               NO. OF FIELDS = IN(17)                  
c         
c                                                                       
c         
c         A-FORMAT  IFMT(I,1)                                           
c         
c             SEL./REJ.           NSID    FIELDS    CRITERIA            
c         
c                        SEL.                                           
c         
c               KID(2,I)         FIELD     ISEL       IBB               
c         
c                        REJ.                                           
c         
c                                                                       
c         
c         F-FORMAT  IFMT(I,2)                                           
c         
c             VARIABLE NAMES      FIELD                                 
c         
c                                                                       
c         
      goto 20
  148 if (jndex .eq. (-7)) goto 475
      lindex = 2
      goto 150
  149 if (jndex .eq. (-7)) goto 475
      lindex = 1
      irecs = nchar / 80
      do 1479 i = 1, irecs
 1479 backspace(unit=60, err=150) 
  150 if (jli .eq. 2) in(24) = niaz
      write(unit=61, fmt=155) in(2), in(17), ifmt(1)
  155 format(/1x,21hINPUT SPECIFICATIONS:,5x,i3,18h CHARS. PER RECORD
     &/1x,26x,16hNO. OF FIELDS = ,i3//10x,8hA-FORMAT,5x,a108/)
      if (in(6) .ne. 0) write(unit=61, fmt=157) in(6)
  157 format(22x,15hPEDIGREE  FIELD,i5)
      write(unit=61, fmt=800) in(10), in(5), in(14), in(13)
  800 format(22x,15hFAMILY    FIELD,i5/22x,15hPOSITION  FIELD,i5/22x,
     &15hGEN.      FIELD,i5/22x,15hPROBAND   FIELD,i5)
      if (in(12) .ne. 0) write(unit=61, fmt=656) in(12)
  656 format(22x,15hQUANT.    FIELD,i5)
      if (in(11) .ne. 0) write(unit=61, fmt=659) in(11), jaf
  659 format(22x,15hAFFECT.   FIELD,i5,10x,3hA =,i3,7x,3hN =,i3,7x,3hU =
     &,i3)
      if (jli .eq. 0) goto 158
      if (jli .ne. 2) goto 159
      write(unit=61, fmt=154) in(1), (zi(i),i = 1, in(24))
  154 format(/22x,15hLIABILITY FIELD,i5/27x,11hLZ VALUES: ,10f9.5/)
      goto 158
  159 write(unit=61, fmt=651) in(1), (yi(i),i = 1, in(24))
  651 format(/22x,15hLIABILITY FIELD,i5/27x,11hLI VALUES: ,10f9.5/)
      if (jli .eq. 1) goto 158
      write(unit=61, fmt=156) zi(1)
  156 format(27x,11hLZ VALUE = ,f10.5/)
  158 if (in(20)) 652, 660, 652
  652 write(unit=61, fmt=653) in(20), jsex
  653 format(/22x,12hSEX  FIELD  ,i5,5x,7hMALE = ,i5,5x,9hFEMALE = ,i5/)
  660 write(unit=61, fmt=654) in(15), in(16), in(21), in(22), in(23)
  654 format(/27x,11hPO VALUES: ,5x,9hFATHER = ,i5,10x,9hMOTHER = 
     &,i5/43x,17hFATHER POINTER = ,i5,10x,17hMOTHER POINTER = ,i5,10x,
     &19hCHILDREN POINTER = ,i5/43x,21hCHILD = ANYTHING ELSE/)
      if (nfi .ne. 0) write(unit=61, fmt=691) (fi(i),i = 1, nfi)
  691 format(27x,5hFI = ,10f9.5/)
      if (npie .eq. 0) goto 670
      write(unit=61, fmt=669) (pie(i),i = 1, npie)
  669 format(27x,5hPI = ,10f9.5/)
      goto 6720
  670 write(unit=61, fmt=671) 
  671 format(27x,15hNO PI SPECIFIED/)
 6720 if (izy) 6721, 672, 6721
 6721 write(unit=61, fmt=6723) zeey
 6723 format(27x,5hZY = ,f9.5/)
  672 nsid = in(8) + in(9)
      if (nsid) 165, 200, 165
  165 write(unit=61, fmt=170) nsid
  170 format(20x,9hSEL./REJ.,8x,i5,5x,6hFIELDS,4x,8hCRITERIA)
      nc = 0
      ib(1) = 'SEL.  '
      ib(2) = 'REJ.  '
      do 190 j = 1, nsid
      i1 = nc + 1
      i2 = kid(2,j)
      i3 = isel(i2)
      nc = nc + kid(3,j)
      n = kid(1,j)
  190 write(unit=61, fmt=195) i2, ib(n), i3, (ibb(i),i = i1, nc)
  195 format(23x,i2,3x,a6,3x,5hFIELD,3x,i5,5x,10(a6,2x)/(55x,10(a6,2x)))
  200 write(unit=61, fmt=205) ifmt(2)
  205 format(//10x,8hF-FORMAT,5x,a108/)
  310 write(unit=61, fmt=315) ntrn
  315 format(//18h TRANSFORMATION - ,i3//)
      if (ntrn) 320, 420, 320
  320 ib(1) = '*R.O.*'
      ib(2) = '      '
      ib(5) = 'KIJL  '
      nc = 0
      do 415 i = 1, ntrn
      n = 2
      if (itr(i) .eq. 1) n = 1
      m = itrn(1,i)
      if (m .ne. 41) goto 339
  343 write(unit=61, fmt=344) i
  344 format(//9x,i2,5x,10hTR CODE 41,5x,25hCHANGES LONGITUDE AND LAT,
     &19hITUDE TO KILOMETERS)
      nc = nc + 1
      k = cc(nc)
      goto 353
  339 write(unit=61, fmt=340) i, m, ncode(m), ib(n)
  340 format(//1x,8x,i2,5x,7hTR CODE,i3,3x,5hX(K)=,a62,3x,a6)
  353 do 365 j = 2, 4
      if (itrn(j,i)) 355, 365, 355
  355 write(unit=61, fmt=360) ib(5)(j - 1:j - 1), itrn(j,i)
  360 format(29x,2hX(,a1,2h)=,i3)
      if ((m .eq. 41) .and. (j .eq. 2)) write(unit=61, fmt=360) ib(5)(4:
     &4), k
  365 continue
      k = itrn(5,i)
      if (k) 373, 385, 373
  373 if (m .eq. 41) k = k - 1
      i1 = nc + 1
      nc = nc + k
      write(unit=61, fmt=380) k, (cc(j),j = i1, nc)
  380 format(20x,i6,3x,7hC-VALUE,3x,6hC(M) =,6f12.6/(45x,6f12.6))
  385 if (itrn(6,i)) 390, 415, 390
  390 i1 = nc + 1
      nc = nc + itrn(6,i)
      if (itrn(7,i)) 395, 405, 395
  395 write(unit=61, fmt=400) itrn(6,i), (icdc(j),j = i1, nc)
  400 format(20x,i6,21h   A-VALUE   A(M) =  ,10(a6,2x)/(47x,10(a6,2x)))
      goto 415
  405 write(unit=61, fmt=410) itrn(6,i), (cc(j), icc(j),j = i1, nc)
  410 format(20x,i6,21h   A-VALUE   A(M) =  
     &,6(f10.4,1x,a1)/(47x,6(f10.4,1x,a1)))
  415 continue
  420 ierror = 0
      if ((in(10) .ne. 0) .and. (in(5) .ne. 0)) goto 602
      write(unit=61, fmt=601) 
  601 format(1x,21hNO  ID  OR  PO  FIELD)
      ierror = 1
  602 if (in(1) .eq. 0) goto 604
      if (in(24) .ne. 0) goto 604
      write(unit=61, fmt=603) 
  603 format(1x,21hNO LI/LZ CONTROL CARD)
      ierror = 1
  604 if ((in(11) .ne. 0) .or. (in(12) .ne. 0)) goto 608
      write(unit=61, fmt=607) 
  607 format(1x,26hMISSING AF AND/OR QU FIELD)
      ierror = 1
  608 if (ierror) 475, 4420, 475
 4420 iguana(1) = 2
      goto 480
  475 iguana(1) = 0
      lindex = 2
  480 return 
      end
      subroutine over2()
c     PROGRAM SAMENAME                                                  
c         
c     ITEMP = 18 + 23 *6 :  23 = MAX. NO. OF RECORDS PER FAMILY         
c         
c                                                                       
c         
      implicit double precision (o-z, a-h)
c                DEBUG=1 FOR DEBUGGING PRINTOUT, ELSE=0                 
c         
c                NULLS = HIGHEST POSSIBLE INTEGER                       
c         
c                                                                       
c         
      parameter (nulls = -8388607)
      character icc*1, ia*1
      character iprel*3, ixi*3
      character icdc*6, ix*6, ibb*6, cfam*6, ccfam*6, cgr*6
      character ifmt*108
      character iai*512
      real ifam, icfam
      integer idebug
      dimension itab1(2, 3), itab2(9, 4, 2), ny(19)
      dimension itemp(156), temp(78), ptemp(9), iptemp(18)
      dimension iporel(3), igepo(3), iprel(35)
      dimension ispace(120), ia(512)
      dimension ixi(120), ix(60), b(3751)
      equivalence (ixi, ix)
      equivalence (itemp, temp), (iptemp, ptemp)
      equivalence (nr, itemp(2)), (nc, itemp(3)), (mat, itemp(4)), (
     &iporel, itemp(5)), (iquad, itemp(8)), (igec, itemp(9)), (igepo, 
     &itemp(10)), (isell, itemp(14))
      equivalence (ili, in(1)), (ichar, in(2)), (ip, in(4)), (mbr, in(5)
     &), (igr, in(6)), (iden, in(10)), (iaf, in(11)), (iqt, in(12)), (
     &ipro, in(13)), (igen, in(14)), (ifa, in(15)), (imo, in(16)), (nin
     &, in(17)), (np, in(19)), (isex, in(20)), (ifap, in(21)), (imop, in
     &(22)), (kidp, in(23)), (nia, in(24))
      equivalence (b, ib), (iai, ia)
      common /pets/ cat(12), dog, iguana(27), idebug
      common /convr1/ in(25), inf(60), isel(12), kid(3, 12), ntrn, itrn(
     &7, 60), itr(60), cc(500), jp, inv(60), ipp, iq, itot
      common /convr2/ ifmt(2), ibb(100), icdc(500), icc(501)
      common /poin/ yi(10), zi(10), pre(7), id(7), jaf(3), jsex(2), 
     &jpro(4), npp, iccc, jli
      common /oddsn/ pie(10), fi(10), zeey, npie, izy, nfi
      common /connet/ lindex
      common /add/ ib(7502), x(60), index
      common /xes/ ix
      data iprel / '1A', '1B', '1C', '1D', '1E', '1F', '1G', '2A', '2B'
     &, '2C', '2D', '2E', '2F', '2G', '3A', '3B', '3C', '3D', '3E', '3F'
     &, '3G', '4A', '4B', '4C', '4D', '4E', '4F', '4G', '5A', '5B', '5C'
     &, '5D', '5E', '5F', '5G' /
      inp = 49
cc JH Zhao 6/6/2001
      data ny/19*0/,itab1/6*0/,icfam/0./,icmem/0/,istar2/0/,istart/0/,
     + jktemp/0/,item1m/0/,item2m/0/,istapt/0/,istap2/0/,iccc/0/
      
  610 rewind(unit=49) 
      write(unit=62, fmt=880) 
  880 format(///)
      do 607 j = 1, 3
  607 itab1(1,j) = 0
      do 608 i = 1, 9
      do 608 j = 1, 4
      itab2(i,j,1) = 0
  608 itab2(i,j,2) = 0
      nread = 0
      iguana(2) = 0
      x9 = -9999.
      do 22 i = 1, 7502
   22 ib(i) = nulls
cNO GROUP FIELD   
      itypgr = 0
cID FIELD IS ALPHA
      itypaf = 1
cGROUP FIELD  IS A
      if (igr .ne. 0) itypgr = 1
      do 40 i = 1, jp
cID FIELD IS NUMER
      if (iden .eq. inf(i)) itypaf = 2
cGROUP FIELD IS AL
      if ((igr .ne. 0) .and. (igr .eq. inf(i))) itypgr = 2
   40 continue
      ifam = -10.
      cfam = '\\\\\\\\\\\\'
      cgr = '      '
      ilast = 0
      ilast1 = 0
      nout = 0
      nfam = 0
      nfamb = 0
      nsid = in(8) + in(9)
      icount = 0
      nrej = 0
      ntr = 0
c     INITIAL TAPE BLOCK                                                
c         
      itab = 0
   20 read(unit=inp, fmt='(512A1)', end=350) (ia(i),i = 1, ichar)
      goto 110
   37 iguana(2) = 1
      lindex = 2
      goto 500
  110 read(unit=iai, fmt=ifmt(1)) (ix(i),i = 1, nin)
      if (idebug .ne. 1) goto 130
  335 write(unit=61, fmt=336) (ia(i),i = 1, ichar)
c     COUNT TOTAL NO. OF INPUT RECORDS                                  
c         
  336 format(1x,6hINPUT:,120a1)
c     SELECTION / REJECTION                                             
c         
  130 icount = icount + 1
      if (nsid) 140, 195, 140
  140 itime = 1
      l = 0
      do 165 i = 1, nsid
      l = l + kid(3,i)
      k = kid(2,i)
      k = isel(k)
      do 150 j = itime, l
      if (ix(k) .eq. ibb(j)) goto 160
  150 continue
      if (kid(1,i) .ne. 1) goto 165
  155 nrej = nrej + 1
      goto 777
  160 if (kid(1,i) .ne. 1) goto 155
  165 itime = l + 1
  195 iblank = 0
      do 220 i = 1, jp
      k = inv(i)
      if (ix(k) .ne. '      ') goto 220
  215 iblank = iblank + 1
      ispace(iblank) = inf(i)
  220 continue
      read(unit=iai, fmt=ifmt(2)) (x(i),i = 1, jp)
      do 230 i = jp, 1, -1
      k = inv(i)
  230 x(k) = x(i)
      do 260 i = 1, iblank
      j = ispace(i)
c                                                                       
c         
c     CALL TRANSFORMATION                     PROCESSS A RECORD         
c         
c                                                                       
c         
  260 x(j) = x9
      if (ntrn) 240, 250, 240
  240 call tr
      if (index) 250, 245, 250
  245 ntr = ntr + 1
      goto 777
c                                                                       
c         
c     CHECK FOR CODING ERRORS IN PO,AF,SEX,AND PROBAND FIELDS.          
c         
c                                                                       
c         
  250 itab = itab + 1
      if (iaf) 1001, 1003, 1001
 1001 inter = x(iaf)
      if (((inter .eq. jaf(1)) .or. (inter .eq. jaf(2))) .or. (inter
     & .eq. jaf(3))) goto 1003
      write(unit=61, fmt=1002) 
 1002 format(1x,14hWRONG AF VALUE)
 1005 write(unit=61, fmt=2005) icount
 2005 format(1x,13hERROR RECORD:,i5)
      write(unit=61, fmt=336) (ia(i),i = 1, ichar)
      goto 37
 1003 if (ili) 1006, 1006, 1013
 1013 inter = x(ili)
      if ((inter .gt. 0) .and. (inter .le. in(24))) goto 1006
      write(unit=61, fmt=1004) 
 1004 format(1x,
     &57hWRONG LI VALUE GIVEN THE VALUES SPECIFIED IN THE LI CARD.)
      goto 1005
 1006 if (in(20)) 1000, 1000, 1007
 1007 inter = x(isex)
      if ((inter .eq. jsex(1)) .or. (inter .eq. jsex(2))) goto 1000
      write(unit=61, fmt=1008) 
 1008 format(1x,15hWRONG SEX VALUE)
c     THIS SECTION HAS BEEN ALTERED TO MAKE IT MORE COMPATIBLE          
c         
c     ALPHA FIELDS ARE NOW HANDLED SEPARATELY                           
c         
c                                                                       
c         
      goto 1005
 1000 icmem = x(mbr)
      if (itypaf .eq. 2) then
      icfam = x(iden)
      if (ifam .eq. (-10)) goto 678
      else
      ccfam = ix(iden)
      if (cfam .eq. '\\\\\\\\\\\\') goto 678
      end if
  651 if (itypaf .eq. 2) then
      if (ifam .ne. icfam) goto 675
      else
      if (cfam .ne. ccfam) goto 675
      end if
  652 temp(8) = icfam
      if (itypgr .eq. 2) then
      temp(9) = x(igr)
      else
      if (itypgr .eq. 1) cgr = ix(igr)
c     END OF CHANGES FOR COMPATIBILITY                                  
c         
c                                                                       
c         
      end if
c FATHER'S RECORD                                                       
c         
      nr = nr + 1
      if (icmem .ne. ifa) goto 660
      mat = mat + 1
  655 if (iqt .eq. 0) goto 658
      if (x(iqt) .eq. x9) goto 656
      if (izy .eq. 1) goto 6658
c                                                                       
c         
c     IF AF & QT FIELDS, AND QT .NE. BLANK       04-06-81.              
c         
c                                                                       
c         
 6656 if (iaf) 6657, 6655, 6657
 6657 temp(istar2 + 1) = x(iqt)
      isytem = x(iaf)
      if (isytem .eq. jaf(1)) itemp(istart + 3) = 5
      if (isytem .eq. jaf(2)) itemp(istart + 3) = 4
      if (isytem .eq. jaf(3)) itemp(istart + 3) = 3
      if (isytem .eq. jaf(1)) jktemp = 1
      if (isytem .eq. jaf(2)) jktemp = 0
      if (isytem .eq. jaf(3)) jktemp = 2
      if (icmem .eq. ifa) item1m = jktemp
      if (icmem .eq. imo) item2m = jktemp
      goto 6593
 6658 if (x(iqt) .le. zeey) goto 6656
      itemp(istart + 3) = 1
      goto 659
 6655 temp(istar2 + 1) = x(iqt)
      itemp(istart + 3) = 3
      goto 659
  656 if (iaf) 658, 654, 658
  654 itemp(istart + 3) = 2
      goto 659
  658 isytem = x(iaf)
      if (isytem .eq. jaf(1)) itemp(istart + 3) = 1
      if (isytem .eq. jaf(2)) itemp(istart + 3) = 0
      if (isytem .eq. jaf(3)) itemp(istart + 3) = 2
  659 if (icmem .eq. ifa) item1m = itemp(istart + 3)
      if (icmem .eq. imo) item2m = itemp(istart + 3)
 6593 if (ili) 6591, 6591, 6590
 6591 itemp(istart + 4) = 0
      goto 6592
 6590 itemp(istart + 4) = x(ili)
 6592 istar2 = istar2 + 3
      istart = istart + 6
c MOTHER'S RECORD                                                       
c         
      goto 777
  660 if (icmem .ne. imo) goto 664
      mat = mat + 2
c FATHER'S POINTER RECORD                                               
c         
      goto 655
  664 if (icmem .ne. ifap) goto 670
      ike = 1
  674 if (ip) 6674, 6673, 6674
 6673 write(unit=61, fmt=676) 
  676 format(1x,29hNO POINTEE RELATIONSHIP FIELD)
      goto 37
 6674 do 6675 i = 1, 35
      if (ixi((2 * ip) - 1) .eq. iprel(i)) goto 6678
 6675 continue
      write(unit=61, fmt=6671) 
 6671 format(1x,16hILLEGAL IP VALUE)
      goto 37
 6678 if (ike .ne. 3) goto 6672
      i = i + 35
 6672 iporel(ike) = i
      if (igen) 6402, 6400, 6402
 6400 igepo(ike) = 1
      if (ike .eq. 3) igepo(ike) = 2
      goto 6401
 6402 igepo(ike) = x(igen)
 6401 if (iqt .eq. 0) goto 758
      if (x(iqt) .eq. x9) goto 756
      if (izy .ne. 1) goto 6416
      if (x(iqt) .le. zeey) goto 6416
      iptemp(istapt + 3) = 1
      goto 759
 6416 if (iaf) 7657, 7655, 7657
 7657 ptemp(istap2 + 1) = x(iqt)
      isytem = x(iaf)
      if (isytem .eq. jaf(1)) iptemp(istapt + 3) = 5
      if (isytem .eq. jaf(2)) iptemp(istapt + 3) = 4
      if (isytem .eq. jaf(3)) iptemp(istapt + 3) = 3
      goto 759
  756 if (iaf) 758, 757, 758
  757 iptemp(istapt + 3) = 2
      goto 759
 7655 ptemp(istap2 + 1) = x(iqt)
      iptemp(istapt + 3) = 3
      goto 759
  758 isytem = x(iaf)
      if (isytem .eq. jaf(1)) iptemp(istapt + 3) = 1
      if (isytem .eq. jaf(2)) iptemp(istapt + 3) = 0
      if (isytem .eq. jaf(3)) iptemp(istapt + 3) = 2
  759 if (ili) 7591, 7591, 7590
 7591 iptemp(istapt + 4) = 0
      goto 7592
 7590 iptemp(istapt + 4) = x(ili)
 7592 istap2 = istap2 + 3
      istapt = istapt + 6
c MOTHER'S POINTER RECORD                                               
c         
      goto 777
  670 if (icmem .ne. imop) goto 657
      ike = 2
c CHIDREN'S POINTER RECORD                                              
c         
      goto 674
  657 if (icmem .ne. kidp) goto 700
      ike = 3
c CHILD'S RECORD                                                        
c         
      goto 674
  700 nc = nc + 1
      if (nc .le. 20) goto 720
      write(unit=61, fmt=721) 
  721 format(1x,30hA FAMILY HAS OVER 20 CHILDREN.)
      goto 37
  720 if (igen) 6404, 6403, 6404
 6403 igec = 2
      goto 6405
c     ISELL=0 COMPLETE SELECTION                                        
c         
c     ISELL>1 INCOMPLETE SELECTION                                      
c         
 6404 igec = x(igen)
 6405 if (ipro) 655, 655, 6406
 6406 if (x(ipro) .lt. 1.) goto 655
      if (x(iaf) .eq. jaf(1)) goto 1035
      write(unit=61, fmt=1034) 
 1034 format(1x,22hNORMAL KID W/ PROBAND.)
      goto 1005
 1035 isell = x(ipro)
      if (isell .le. npie) goto 655
      write(unit=61, fmt=1030) 
 1030 format(1x,36hPROBAND VALUE > THAN THE NO. OF PIS.)
c                                                                       
c         
c     PROCESS A FAMILY                                                  
c         
c                                                                       
c         
      goto 1005
  675 lend = 18 + (6 * nr)
c     NEW FILE FOR COMPATIBILITY                                        
c         
      if ((nout + lend) .gt. 7500) goto 680
      if ((itypaf .eq. 1) .or. (itypgr .eq. 1)) write(unit=52, fmt=529) 
     &cfam, cgr
c                                                                       
c         
  529 format(2a6)
      nfam = nfam + 1
      nfamb = nfamb + 1
c                                                                       
c         
c     CALCULATE ITAB1,ITAB2                                             
c         
c                                                                       
c         
      ny(nc + 1) = ny(nc + 1) + 1
      isub = itemp(4)
      if (itemp(4) .eq. 0) isub = 1
      if (itemp(4) .eq. 1) isub = 2
      itab1(1,isub) = itab1(1,isub) + 1
      itab1(2,isub) = itab1(2,isub) + nc
      isub1 = 1
      if ((item1m .gt. 2) .or. (item2m .gt. 2)) goto 2900
      isub1 = 4
      goto 3333
 2900 if ((item1m .gt. 2) .and. (item2m .gt. 2)) goto 3333
      if (item1m .gt. 2) isub1 = 2
      if (item2m .gt. 2) isub1 = 3
 3333 icross = 9
      if ((item1m .eq. 2) .and. (item2m .eq. 2)) goto 3335
      if (item1m .ne. 2) goto 2901
      if ((item2m .eq. 0) .or. (item2m .eq. 4)) icross = 7
      if ((item2m .eq. 1) .or. (item2m .eq. 5)) icross = 8
      goto 3335
 2901 if (item2m .ne. 2) goto 2902
      if ((item1m .eq. 0) .or. (item1m .eq. 4)) icross = 3
      if ((item1m .eq. 1) .or. (item1m .eq. 5)) icross = 6
      goto 3335
 2902 if ((item1m .eq. 1) .or. (item1m .eq. 5)) goto 2903
      if ((item2m .eq. 0) .or. (item2m .eq. 4)) icross = 1
      if ((item2m .eq. 1) .or. (item2m .eq. 5)) icross = 2
      goto 3335
 2903 if ((item2m .eq. 0) .or. (item2m .eq. 4)) icross = 4
      if ((item2m .eq. 1) .or. (item2m .eq. 5)) icross = 5
 3335 isub = 2
      if (isell .eq. 0) isub = 1
c                                                                       
c         
c     END OF ITAB1,ITAB2 SECTION                                        
c         
c                                                                       
c         
      itab2(icross,isub1,isub) = itab2(icross,isub1,isub) + 1
      ia1 = 0
      ia2 = 0
      if (iporel(1) .ne. 0) goto 2675
      if ((mat .eq. 0) .or. (mat .eq. 2)) goto 1675
      if (itemp(21) .eq. 2) goto 1675
 2675 ia1 = 1
 1675 if (iporel(2) .ne. 0) goto 2676
      if (mat - 2) 1676, 2680, 2681
 2680 if (itemp(21) .eq. 2) goto 1676
      goto 2676
 2681 if (itemp(27) .eq. 2) goto 1676
 2676 ia2 = 1
 1676 if ((ia1 + ia2) - 1) 1677, 1678, 1679
 1677 itemp(8) = 3
      goto 1680
 1678 itemp(8) = 2
      goto 1680
 1679 itemp(8) = 1
 1680 if (istapt) 711, 712, 711
  711 do 710 i = 1, istapt
      k = istart + i
  710 itemp(k) = iptemp(i)
  712 nout1 = nout + 1
      nout2 = nout + 14
      nout3 = (nout / 2) + 8
      nout4 = (nout / 2) + 9
c     TO PREVENT USE OF JOINT LIKELIHOOD WHEN THERE ARE POINTERS        
c         
      itemp(1) = (10 * item1m) + item2m
      if (((iporel(1) .ne. 0) .or. (iporel(2) .ne. 0)) .or. (iporel(3)
     & .ne. 0)) itemp(13) = 0
      do 6676 i = 1, lend
      nout = nout + 1
 6676 ib(nout) = itemp(i)
      if (idebug .ne. 1) goto 678
      if (itypaf - 1) 41, 41, 42
   42 write(unit=61, fmt=43) (ib(i),i = nout1, nout2), b(nout3), nfam, 
     &nout1
   43 format(1x,8hOUTPUT: ,12i5,i10,i5,5x,f10.0,3x,7hNFAM = ,i5,4x,
     &10h1ST ADDR.=,i4)
      goto 44
   41 write(unit=61, fmt=677) (ib(i),i = nout1, nout2), cfam, nfam, 
     &nout1
  677 format(1x,8hOUTPUT: ,12i5,i10,i5,5x,a6,3x,7hNFAM = ,i5,4x,
     &10h1ST ADDR.=,i4)
   44 if (itypgr .eq. 1) write(unit=61, fmt=441) cgr
      if (itypgr .eq. 2) write(unit=61, fmt=442) b(nout4)
  441 format(100x,14hPEDIGREE NO. =,a6)
  442 format(100x,14hPEDIGREE NO. =,f10.0)
      j = nout4
      do 673 i = 1, nr
      k = (2 * j) + 3
      write(unit=61, fmt=672) b(j + 1), ib(k), ib(k + 1)
  672 format(6x,f12.6,2i4)
c     GET READY FOR A NEW FAMILY                                        
c         
  673 j = j + 3
  678 do 671 i = 1, 154
  671 itemp(i) = 0
      isell = iccc
      itemp(13) = iccc
      item1m = 2
      item2m = 2
      istar2 = 9
      istart = 18
      istap2 = 0
c     CHANGES FOR COMPATIBILITY                                         
c         
      istapt = 0
      if (itypaf .eq. 2) then
      ifam = icfam
      else
      cfam = ccfam
c     END OF CHANGES                                                    
c         
      end if
      if (ilast) 6677, 651, 6677
c                                                                       
c         
c     BUFFER OUT A BLOCK                                                
c         
c                                                                       
c         
 6677 ilast1 = 1
  
  680 ib(7501) = nfamb
      write(unit=48) (ib(i),i = 1, 7502)
  682 do 683 i = 1, 7502
  683 ib(i) = nulls
      nread = nread + 1
      nout = 0
      nfamb = 0
      if (ilast) 5002, 675, 5002
5002      if (ilast1 .eq. 1) goto 351
      goto 675
  777 goto 20
  350 ilast = 1
      goto 675
  351 endfile(unit=49) 
  
  435 do 436 i = 1, 114
  436 ib(i) = 0
      ib(1) = nread
      ib(2) = nia
      ib(3) = npie
      ib(4) = nfi
      do 490 i = 1, 27
c 490 IB(I+26)=ITITLE(I)                                                
c         
  490 ib(i + 26) = 0
      ib(54) = itypgr
      do 491 i = 1, 10
  491 b(i + 27) = pie(i)
      do 492 i = 1, 10
  492 b(i + 37) = fi(i)
      do 493 i = 1, 10
  493 b(i + 47) = zi(i)
      ib(25) = jli
      ib(26) = itypaf
      if (idebug .eq. 1) then
      do 494 i = 1, 6
  494 write(unit=61, fmt=495) (ib(j),j = 1 + (12 * (i - 1)), 12 * i)
  495 format(1h ,12i10)
      end if
      do 437 i = 1, 10
  437 b(i + 2) = yi(i)
      inulls = nulls
      write(unit=48) inulls, (ib(i),i = 1, 114), (0,i = 1, 7387)
  439 continue
      write(unit=62, fmt=465) icount, nrej, ntr, itab, nfam, nread
  465 format(5x,22hTOTAL INPUT RECORDS = ,i10//10x,11hREJECTED BY//15x,
     &8hSI/SD = ,i14//15x,5hTR = ,10x,i7/10x,20hRECORDS PROCESSED = 
     &,i7//12x,18hNO. OF FAMILIES = ,i7//7x,23hNO. OF OUTPUT BLOCKS = 
     &,i7//)
      write(unit=62, fmt=1599) 
 1599 format(1x,11hFAMILY SIZE)
      write(unit=62, fmt=1600) (i,i = 0, 18)
 1600 format(1x,21(2x,i2,1h:))
      write(unit=62, fmt=1601) ny
 1601 format(/1x,19i5//)
      write(unit=62, fmt=1602) 
 1602 format(21x,15hNO. OF FAMILIES,6x,15hNO. OF CHILDREN)
      write(unit=62, fmt=1603) itab1
 1603 format(1x,11hNO PARENTS ,i15,i20/1x,11hONE PARENT ,i15,i20/1x,
     &11hTWO PARENTS,i15,i20//)
      write(unit=62, fmt=1643) 
 1643 format(//1x,22hCOMPLETE ASCERTAINMENT)
      write(unit=62, fmt=1645) ((itab2(i,j,1),i = 1, 9),j = 1, 4)
 1645 format(/1x,10x,45h  NXN  NXA  NXU  AXN  AXA  AXU  UXN  UXA  UXU
     &//1x,7x,3hQXQ,9i5/1x,7x,3hQX?,9i5/1x,7x,3h?XQ,9i5/1x,7x,3h?X?,9i5)
      write(unit=62, fmt=1644) 
 1644 format(//1x,24hINCOMPLETE ASCERTAINMENT)
      write(unit=62, fmt=1645) ((itab2(i,j,2),i = 1, 9),j = 1, 4)
500      return 
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::::::: 
      subroutine sidcx()
c                                                                       
c         
c     SI (I) (J1, J2, J3, . . )                                         
c         
c     SD (I) (J1, J2, J3, . . )                                         
c         
c     THE FIRST SET OF PARENTHESES GIVES THE FIELD SEQUENCE FOR         
c         
c     SELECTION AND REJECTION.                                          
c         
c     THE SECOND SET GIVES THE ALPHANUMERIC CRITERIA FOR SELECTION OR   
c         
c     REJECTION.                                                        
c         
c                                                                       
c         
      implicit double precision (o-z, a-h)
      character a*1, icc*1, d*1
      character icdc*6, ibb*6
      character ifmt*108
      dimension d(45)
      common /convr1/ in(25), inf(60), isel(12), kid(3, 12), ntrn, itrn(
     &7, 60), itr(60), cc(500), jp, inv(60), ip, iq, itot
      common /convr2/ ifmt(2), ibb(100), icdc(500), icc(501)
      common /inter1/ index, icol, i1, i2, nchar
      common /inter2/ a(800)
      data kk /0/
c     KID(1,IP) = 1, IF SI CARD                                         
c         
c               = 2, IF SD CARD                                         
c         
c     KID(2,IP) = IDENTIFIES WHICH SELECTION/REJECTION FIELD            
c         
c     KID(3,IP) = TOTAL NUMBER OF SELECTION (OR DELETION) CRITERIA ON   
c         
c                 IP CARD.                                              
c         
c     IBB = STORES CRITERIA FOR SELECTION OR DELETION.                  
c         
c                                                                       
c         
c     IP COUNTS NO. OF SI/SD CARDS.                                     
c         
c     KID(1,IP) = 1 IF SI CARD AND = 2 IF SD CARD.                      
c         
c                                                                       
c         
      data d / '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', ' ', 
     &',', '.', '-', '(', ')', '=', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 
     &'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
     &'U', 'V', 'W', 'X', 'Y', 'Z', '*', ' ' /
      ip = ip + 1
      if (a(icol) .ne. d(26)) goto 30
   25 kid(1,ip) = 1
      goto 40
c                                                                       
c         
c     STORE INDEX IN FIRST SET OF PARENTHESES INTO KID(2,IP)            
c         
c                                                                       
c         
   30 kid(1,ip) = 2
   40 call parn
      jhel = xno2n(icol)
      kid(2,ip) = jhel
      if ((jhel .le. 0) .or. (jhel .gt. 12)) goto 160
   60 if (i2) 70, 65, 70
   65 icol = icol - 1
   70 call parn
      if (index) 75, 160, 75
c                                                                       
c         
c     IQ IS CUMULATIVE COUNT OF ALL SEL./REJ. CRITERIA.                 
c         
c                                                                       
c         
c     AN * INDICATES A BLANK                                            
c         
c                                                                       
c         
   75 na = iq
   80 j = 0
      na = na + 1
   85 if (icol .gt. i2) goto 140
   90 if (a(icol) .eq. d(12)) goto 110
   95 if (a(icol) .ne. d(44)) goto 105
  100 a(icol) = d(11)
  105 j = j + 1
      ibb(na)(j:j) = a(icol)
      kk = 1
      call nspacen
      goto 85
  110 kk = 2
      call nspacen
  115 if (j - 6) 120, 130, 160
  120 j = j + 1
  130 if (kk .lt. 3) goto 80
      goto 155
  140 if (kk - 1) 150, 145, 150
  145 kk = 3
      goto 115
c                                                                       
c         
c     STORE NO. OF CRITERIA ON CARD IN KID(3,IP)                        
c         
c                                                                       
c         
  150 na = na - 1
  155 kid(3,ip) = na - iq
      iq = na
      index = 1
      return 
  160 index = 0
      return 
      end
      subroutine tr()
      implicit double precision (o-z, a-h)
      character icc*1
      character icdc*6, ibb*6, ixx*6
      character ifmt*108
      dimension ixx(60)
      common /convr1/ in(25), inf(60), isel(12), kid(3, 12), ntrn, itrn(
     &7, 60), itr(60), cc(500), jp, inv(60), ip, iq, itot
      common /convr2/ ifmt(2), ibb(100), icdc(500), icc(501)
      common /add/ ib(7502), x(60), index
      common /xes/ ixx
      x9 = -9999.
      nn = 0
      do 410 izip = 1, ntrn
      i = izip
      nt = itrn(1,i)
      ik = itrn(2,i)
      ii = itrn(3,i)
      nc = itrn(5,i)
      na = itrn(6,i)
      nn1 = nn + 1
      i = itrn(4,i)
      if (x(ii) .ne. x9) goto 5001
 5000 x(ik) = x9
      goto 380
 5001 if (i) 5002, 5010, 5002
c     CHECK TRANSFORMATION CODE                                         
c         
 5002 if (x(i) .eq. x9) goto 5000
 5010 if (nt .lt. 8) goto 20
    5 if (nt .lt. 15) goto 21
   10 if (nt .lt. 20) goto 22
c     TR (1 - 7)                                                        
c         
      goto 200
c     TR (8 - 14)                                                       
c         
   20 goto (25, 30, 35, 35, 65, 90, 110), nt
   21 nt = nt - 7
c     TR (15 - 19)                                                      
c         
      goto (115, 120, 125, 125, 155, 160, 165), nt
   22 nt = nt - 14
c     GENERAL OR ARITHMETIC TRANSFORMATIONS                             
c         
c     MULTIPLICATION ( TR = 1 OR 15)                                    
c         
      goto (26, 31, 61, 116, 121, 420), nt
   25 x(ik) = x(ii) * cc(nn1)
      goto 400
   26 x(ik) = x(ii)
c     MULTIPLICATION (TR = 2 OR 16)                                     
c         
      goto 400
   30 x(ik) = (x(ii) * x(i)) * cc(nn1)
      goto 400
   31 x(ik) = x(ii) * x(i)
c     DIVISION (TR = 3, 4, OR 17)                                       
c         
      goto 400
   35 if (x(ii)) 50, 40, 50
   40 if (nc .ne. 2) goto 420
   45 x(ik) = cc(nn1 + 1)
      goto 400
   50 if (nt .ne. 3) goto 60
   55 x(ik) = cc(nn1) / x(ii)
      goto 400
   60 x(ik) = (cc(nn1) * x(i)) / x(ii)
      goto 400
   61 if (x(ii)) 62, 420, 62
   62 x(ik) = x(i) / x(ii)
c     POWER (TR = 5)                                                    
c         
      goto 400
   65 if (x(ii) .le. 0.) goto 40
   85 x(ik) = exp(cc(nn1) * log(x(ii)))
c     POWER (TR = 6)                                                    
c         
      goto 400
   90 if (x(ii)) 105, 95, 105
   95 if (nc .ne. 3) goto 420
  100 x(ik) = cc(nn1 + 2)
      goto 400
  105 x(ik) = exp((cc(nn1 + 1) * x(ii)) * log(cc(nn1)))
c     EXPONENT (TR = 7)                                                 
c         
      goto 400
  110 x(ik) = exp(x(ii))
c     ADDITION (TR = 8 OR 18)                                           
c         
      goto 400
  115 x(ik) = x(ii) + (x(i) * cc(nn1))
      goto 400
  116 x(ik) = x(ii) + x(i)
c     ADDITION (TR = 9 OR 19)                                           
c         
      goto 400
  120 x(ik) = cc(nn1) + (cc(nn1 + 1) * x(ii))
      goto 400
  121 x(ik) = x(ii) + cc(nn1)
c     LOGARITHMS (TR = 10, OR 11)                                       
c         
      goto 400
  125 if (x(ii) .gt. 0.) goto 140
  130 if (nc) 135, 420, 135
  135 x(ik) = cc(nn1)
      goto 400
  140 if (nt .ne. 3) goto 150
  145 x(ik) = log(x(ii))
      goto 400
  150 x(ik) = log10(x(ii))
c     TRIGONOMETRIC FUNCTIONS (TR = 12, 13, OR 14)                      
c         
      goto 400
  155 x(ik) = sin(x(ii) * cc(nn1))
      goto 400
  160 x(ik) = cos(x(ii) * cc(nn1))
      goto 400
  165 x(ik) = asin(x(ii) * cc(nn1))
c     SPECIAL TRANSFORMATION                                            
c         
      goto 400
  200 if (nt .ge. 23) goto 250
  201 i = (nn + nc) + 1
c     NUMERICAL DICHOTOMY (TR = 21)                                     
c         
      if (nt .ne. 21) goto 220
  205 if (x(ii) .le. cc(i)) goto 135
c     NUMERICAL TRICHOTOMY (TR = 22)                                    
c         
      goto 45
  220 if (x(ii) .lt. cc(i)) goto 135
  230 if (x(ii) .ge. cc(i + 1)) goto 100
      goto 45
  250 nt = nt - 22
      n = nn + nc
      if (nt - 6) 251, 340, 251
  251 if (nt .eq. 19) goto 395
  252 j = 0
      k = 0
      y = 0.
  255 j = j + 1
      k = k + 1
      i = izip
c     ALPHABETIC SPECIFICATIONS                                         
c         
      if (itrn(7,i)) 260, 265, 260
  260 i = n + j
      if (ixx(ii) .eq. icdc(i)) goto 310
c     NUMERICAL SPECIFICATIONS                                          
c         
      goto 285
  265 i = n + j
c     NUMERICAL INTERVAL                                                
c         
      if (icc(i) .ne. '-') goto 275
  270 j = j + 1
      if ((x(ii) .ge. cc(i)) .and. (x(ii) .lt. cc(i + 1))) goto 310
c     NON-INTERVAL                                                      
c         
      goto 285
c     NON-MATCH                                                         
c         
  275 if (x(ii) .eq. cc(i)) goto 310
  285 if (j - na) 255, 288, 255
c     REJECTION OR X(IK) = C(2)                                         
c         
  288 goto (290, 300, 315, 312, 330, 420, 420), nt
  290 i = izip
c     X(K) = C(2)                                                       
c         
      if (itr(i)) 420, 295, 420
  295 x(ik) = cc(nn1 + 1)
c     X(K) = C(1)                                                       
c         
      goto 380
  300 x(ik) = cc(nn1)
c     MATCH                                                             
c         
      goto 380
c     REJECTION OR X(IK) = C(1)                                         
c         
  310 goto (300, 290, 312, 315, 320, 420, 370), nt
  312 i = izip
c     X(K) = X(I)                                                       
c         
      if (itr(i)) 420, 300, 420
  315 x(ik) = x(ii)
c     X(K) = C(J)                                                       
c         
      goto 380
  320 i = nn + k
  325 x(ik) = cc(i)
c     REJECTION OR X(K) = C(N+1)                                        
c         
      goto 380
  330 i = izip
      if (itr(i)) 420, 335, 420
  335 i = n
c     EQUAL INTERVALS                                                   
c         
      goto 325
  340 i = n
      if (x(ii) .lt. cc(i + 1)) goto 300
  345 j = 0
  350 j = j + 1
      y = cc(i + 1) + (j * cc(i + 3))
      if (y .ge. cc(i + 2)) goto 365
  355 if (x(ii) .ge. y) goto 350
  360 x(ik) = cc(nn1) + (j * cc(nn1 + 1))
      goto 400
  365 if (x(ii) .lt. cc(i + 2)) goto 360
c     X(I) = A(N)                                                       
c         
      goto 100
  370 if (j .eq. na) goto 372
  371 j1 = nn + 2
      goto 373
  372 j1 = nn + 3
  373 n1 = n + 1
      n2 = (n + na) - 1
      ii = ik
      do 375 i = n1, n2
      if (icc(i) .eq. '-') goto 375
  374 x(ii) = cc(j1)
      ii = ii + 1
  375 continue
      if (j .eq. na) goto 380
  376 ik = (ik + k) - 1
      goto 300
  380 i = izip
      if (itrn(7,i)) 390, 400, 390
  390 nn = (nn + nc) + na
      goto 410
  395 beta = 0.0087266463 * (x(i) * cc(nn1 + 2))
      x(ik) = (x(ii) - cc(nn1 + 1)) * (((111.4175 * cos(beta)) - (.0940
     & * cos(3. * beta))) + (.0002 * cos(5. * beta)))
      ik = cc(nn1)
      x(ik) = (x(i) - cc(nn1 + 2)) * ((111.1363 - (.5623 * cos(2. * beta
     &))) + (.0011 * cos(4. * beta)))
  400 nn = (nn + nc) + na
  410 continue
      index = 1
      goto 430
  420 index = 0
  430 return 
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::::::: 
      subroutine trcx()
c                                                                       
c         
c     A TR CARD FOR TRANSFORMATION IMPLIES AN OPERATION ON VARIABLES    
c         
c     FROM INPUT DATA OR FROM A PRECEDING TRANSFORMATION TO FORM NEW    
c         
c     VALUES REPLACING THE ORIGINAL VALUES OR TO CREATE NEW VARIABLES IN
c         
c     ADDITION TO THE ORIGINAL VARIABLES.                               
c         
c     TR (N) (K) (I,J) (C-VALUES) (A-VALUES)                            
c         
c     N INDICATES WHAT TRANSFORMATION IS TO BE PERFORMED.               
c         
c     K INDICATES THE VARIABLE WHERE THE TRANSFORMED VALUE IS TO BE     
c         
c     STORED.                                                           
c         
c     I (AND J IF NEEDED) INDICATES ON WHICH VARIABLE THE TRANSFORMATION
c         
c     IS TO BE PERFORMED.                                               
c         
c                                                                       
c         
      implicit double precision (o-z, a-h)
      character a*1, icc*1, d*1
      character icdc*6, ibb*6
      character ifmt*108
      dimension d(45)
      common /convr1/ in(25), inf(60), isel(12), kid(3, 12), ntrn, itrn(
     &7, 60), itr(60), cc(500), jp, inv(60), ip, iq, itot
      common /convr2/ ifmt(2), ibb(100), icdc(500), icc(501)
      common /inter1/ index, icol, i1, i2, nchar
      common /inter2/ a(800)
c     ITRN (7,NTRN) = 0 IF A VALUES ARE NUMBERIC                        
c         
c                   = 1 IF A VALUES ARE ALPHABETIC                      
c         
c     ICC = STORES - FOR INTERVAL IN A-VALUE                            
c         
c     ITR = 1 IF TRANSFORMATION HAS REJECTION OPTION                    
c         
c         0 IF IT DOESNT                                                
c         
c                                                                       
c         
c     INTERPRET FIRST FOUR ( ) ON TR CARD AND STORE VALUES              
c         
c                                                                       
c         
      data d / '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', ' ', 
     &',', '.', '-', '(', ')', '=', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 
     &'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
     &'U', 'V', 'W', 'X', 'Y', 'Z', '*', ' ' /, nc/0/, n2/0/
# 37 "trcx.for"
      ntrn = ntrn + 1
      k = 1
      nn = itot
      na = itot
   25 call parn
      if (index) 30, 215, 30
   27 call nspacen
   30 itrn(k,ntrn) = xno2n(icol)
      if (k .ne. 2) goto 40
   35 if (itrn(1,ntrn) .ne. 41) goto 40
   50 if (a(icol) .ne. d(12)) goto 215
   52 na = na + 1
      cc(na) = xno1n(icol)
   40 k = k + 1
      if (k - 4) 25, 45, 65
   45 if (a(icol) .eq. d(12)) goto 27
   60 k = 5
   65 call parn
c                                                                       
c         
c     INTERPRET C-VALUES - 5TH ( ) AND NUMERIC A-VALUES - 6TH ( ) AND   
c         
c          STORE CRITERIA                                               
c         
c                                                                       
c         
# 55 "trcx.for"
      if (index) 70, 110, 70
# 60 "trcx.for"
   70 ii = 0
   75 if (icol .gt. i2) goto 105
   80 if (ii) 85, 100, 85
   85 if (a(icol) .eq. d(12)) goto 97
   90 if (a(icol) .ne. d(14)) goto 215
   95 icc(na) = d(14)
   97 call nspacen
      goto 70
  100 na = na + 1
      cc(na) = xno2n(icol)
      ii = 1
      if (index) 75, 215, 75
  105 if (k - 5) 120, 110, 120
  110 itrn(5,ntrn) = na - itot
      nc = itrn(5,ntrn)
      itot = na
      call parn
c     TRANSFORMATION 28 HAS ONLY A-VALUES, C-VALUES CALCULATED.         
c         
# 77 "trcx.for"
      if (itrn(1,ntrn) .eq. 28) then
# 79 "trcx.for"
      if (index .ne. 0) goto 215
      na = 3
      cc(itot + 1) = cc(itot - 2)
      cc(itot + 2) = cc(itot - 1)
      cc(itot + 3) = cc(itot)
      cc(itot - 2) = 1.
      cc(itot - 1) = 1.
      temp = ((cc(itot + 2) - cc(itot + 1)) / cc(itot + 3)) + 2.
      itemp = temp
      semp = itemp
      if (semp .eq. temp) cc(itot) = itemp
      if (semp .ne. temp) cc(itot) = itemp + 1
      itot = itot + 3
      goto 210
c                                                                       
c         
# 93 "trcx.for"
      end if
# 95 "trcx.for"
      if (index) 115, 114, 115
  114 if (itrn(1,ntrn) .gt. 20) goto 119
      goto 120
  115 k = 6
      if (a(icol) .eq. d(18)) goto 125
      goto 70
  119 if (itrn(1,ntrn) .ne. 41) goto 215
  120 ij = na - itot
      itot = na
      na = ij
c                                                                       
c         
c     INTERPRET ALPHANUMERIC A-VALUES - 6TH ( ) - AND STORE CRITERIA BY 
c         
c          CHARACTERS                                                   
c         
c     AN * INDICATES A BLANK                                            
c         
c                                                                       
c         
# 105 "trcx.for"
      goto 210
# 111 "trcx.for"
  125 na = itot
  130 k = 0
      na = na + 1
  135 call nspacen
      if (icol .gt. i2) goto 190
  140 if (a(icol) .eq. d(12)) goto 160
  145 if (a(icol) .ne. d(44)) goto 155
  150 a(icol) = d(11)
  155 k = k + 1
      icdc(na)(k:k) = a(icol)
      kk = 1
      goto 135
  160 if (k) 165, 135, 165
  165 kk = 2
  170 if (k - 6) 175, 185, 215
  175 k = k + 1
  185 if (kk - 2) 200, 130, 200
  190 if (k) 195, 205, 195
  195 kk = 3
      goto 170
  200 na = na
  205 na = na - itot
      itot = itot + na
      itrn(7,ntrn) = 1
c                                                                       
c         
c     CHECK FOR ALL NECESSARY INFORMATION FOR TRANSFORMATION            
c         
c                                                                       
c         
# 135 "trcx.for"
  210 itrn(6,ntrn) = na
# 139 "trcx.for"
      nt = itrn(1,ntrn)
      kk = itrn(2,ntrn)
      ii = itrn(3,ntrn)
      ij = itrn(4,ntrn)
      if ((((kk .le. 0) .or. (kk .gt. 60)) .or. (ii .le. 0)) .or. (ii
     & .gt. 60)) goto 215
  305 if (nt .ge. 20) goto 375
c                                                                       
c         
c     CHECK FOR CONSTANTS - ARITHMETIC TRANSFORMATIONS 1 TO 14          
c         
c                                                                       
c         
# 146 "trcx.for"
  308 if (nt .gt. 14) goto 600
# 150 "trcx.for"
  310 goto (320, 320, 315, 315, 315, 325, 345, 320, 335, 340, 340, 320
     &), nt
# 151 "trcx.for"
  315 n2 = 2
  320 n1 = 1
      goto 350
  325 if (cc(nn + 1)) 330, 215, 330
  330 n2 = 3
  335 n1 = 2
      goto 350
  340 n2 = 1
  345 n1 = 0
  350 if (nc - n1) 215, 355, 365
  355 goto (375, 375, 360, 360, 360, 360, 375, 375, 375, 360, 360, 375
     &), nt
# 162 "trcx.for"
  360 itr(ntrn) = 1
      goto 375
  365 goto (215, 215, 370, 370, 370, 370, 215, 215, 215, 370, 370, 215
     &), nt
c                                                                       
c         
c     CHECK FOR VARIABLE INDEX J                                        
c         
c                                                                       
c         
# 165 "trcx.for"
  370 if (nc - n2) 215, 375, 215
# 169 "trcx.for"
  375 if (ij) 215, 380, 390
  380 if (nt .ge. 20) goto 405
  385 if (((nt .eq. 2) .or. (nt .eq. 4)) .or. (nt .eq. 8)) goto 215
      goto 600
  390 if (nt .lt. 20) goto 395
  393 if (nt .ne. 41) goto 215
  394 if (nc - 3) 215, 600, 215
  395 if (((nt .ne. 2) .and. (nt .ne. 4)) .and. (nt .ne. 8)) goto 215
  400 if (ij .gt. 60) goto 215
c                                                                       
c         
c     SPECIAL TRANSFORMATIONS 21 - 29 - CHECK FOR NO. OF C-VALUES.      
c         
c                                                                       
c         
# 178 "trcx.for"
      goto 600
# 182 "trcx.for"
  405 nt = nt - 20
      goto (410, 415, 420, 420, 425, 425, 460, 415), nt
  410 n1 = 2
      goto 430
  415 n1 = 3
      goto 430
  420 n1 = 1
      n2 = 2
      goto 430
  425 n1 = 0
      n2 = 1
  430 if (nc - n1) 215, 435, 445
  435 if ((nt .le. 2) .or. (nt .ge. 7)) goto 500
  440 itr(ntrn) = 1
      goto 500
  445 if ((nt .le. 2) .or. (nt .ge. 7)) goto 215
  450 if (nc - n2) 215, 500, 215
  460 if (itrn(7,ntrn)) 465, 480, 465
  465 k = na
  470 if (nc - k) 215, 440, 475
  475 if ((nc - k) - 1) 215, 600, 215
  480 k = 0
      n1 = (nn + nc) + 1
      n2 = (n1 + na) - 1
      do 490 i = n1, n2
      if (icc(i) .eq. d(14)) goto 490
  485 k = k + 1
  490 continue
c                                                                       
c         
c     SPECIAL TRANSFORMATION 21 - 29 - CHECK FOR NO. OF A-VALUES.       
c         
c                                                                       
c         
# 210 "trcx.for"
      goto 470
# 214 "trcx.for"
  500 goto (505, 510, 520, 520, 520, 520, 600, 515, 520), nt
  505 n1 = 1
      goto 525
  510 n1 = 2
      goto 525
  515 n1 = 3
      goto 525
  520 if (na) 600, 215, 600
  525 if (na - n1) 215, 600, 215
  215 index = 0
      write(unit=61, fmt=292) ntrn, (itrn(i,ntrn),i = 1, 7)
  292 format(//1x,11hERROR IN TR,i8,5x,7i4)
      return 
  600 index = 1
      return 
      end
      function xnon(icol)
      implicit double precision (o-z, a-h)
      character a*1, d*1
      character finput*30
      character format*7
      dimension d(10)
      common /inter1/ index, jcol, i1, i2, nchar
      common /inter2/ a(800)
      data d / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
      call nsp1n
c                                                                       
c         
      goto 5
c                                                                       
c         
      entry xno1n(icol)
c                                                                       
c         
      call nspacen
c                                                                       
c         
      entry xno2n(icol)
    5 finput = ' '
      format = '(F  . )'
      index = 1
      n = 0
      l = 0
c                                                                       
c         
c     SEARCHES FOR  NUMBER - REAL OR INTEGER, POSITIVE OR NEGATIVE.     
c         
c          ICOL ON RETURN = 1ST NON-BLANK COLUMN AFTER NUMBER.          
c         
c          INDEX ON RETURN = 1 IF OKAY  ;  = 0 IF SOMETHING'S WRONG.    
c         
c                                                                       
c         
      m = 0
      if (a(icol) .eq. '-') then
      n = 1
      finput(1:1) = a(icol)
      icol = icol + 1
      end if
    6 if (((((((((((a(icol) .eq. '0') .or. (a(icol) .eq. '1')) .or. (a(
     &icol) .eq. '2')) .or. (a(icol) .eq. '3')) .or. (a(icol) .eq. '4'))
     & .or. (a(icol) .eq. '5')) .or. (a(icol) .eq. '6')) .or. (a(icol)
     & .eq. '7')) .or. (a(icol) .eq. '8')) .or. (a(icol) .eq. '9'))
     & .or. (a(icol) .eq. '.')) then
      n = n + 1
      finput(n:n) = a(icol)
      icol = icol + 1
      if (icol .gt. nchar) goto 38
      if (l .eq. 1) m = m + 1
      if (a(icol) .eq. '.') l = 1
      else
      goto 7
      end if
      goto 6
    7 if (n .ne. 0) goto 40
   38 xnon = 0.
      index = 0
      return 
   40 m1 = 0
      if (n .gt. 9) then
      m1 = n / 10
      n = n - (10 * m1)
      end if
      format(3:3) = d(m1 + 1)
      format(4:4) = d(n + 1)
      format(6:6) = d(m + 1)
c                                                                       
c         
c     BRINGS ICOL TO NEXT NON-BLANK CHARACTER                           
c         
      read(unit=finput(:(m1 * 10) + n), fmt=format) xnon
      jcol = icol
      call nsp1n
      end

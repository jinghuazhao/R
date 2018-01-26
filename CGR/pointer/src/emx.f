c         
c     PROGRAM EMX                                                       
c         
      subroutine emx(profile, terfile, isex, fm1, fm2, fm3, fm4)
      implicit double precision (o-z, a-h)
c     IZERO EQUALS '0' IN THE HARRIS COLLATING SEQUENCE                 
c         
c     NULLS IS THE HIGHEST INTEGER                                      
c         
c     DEBUG=1 FOR PRINTOUT OF DEBUGGING INFORMATION, ELSE=0             
c         
c                                                                       
c         
      parameter (izero = 48, nulls = -8388607)
      character sex*6, cfam*6, cgr*6, famnum*6, gr*6
      character jz*3, isex*1, is1*3, is2*3, is3*3, is4*3, ir1*3, ir2*3, 
     &ir3*3, ir4*3, isel4*3
      character title*12,profile*30,terfile*30
      character fm1*81,fm2*81,fm3*81,fm4*81
      logical select
      integer idebug
      dimension ival(30), isel1(9), isel2(4), isel3(5)
      dimension itemp(202), a(3751), title(7), cfam(2, 500), cgr(2, 500)
      dimension jy(9), jz(9), jw(4), sex(2, 2)
      dimension ijp(5),ijz(9),ijw(4)
      equivalence (a, ia)
      dimension is1(9), is2(4), is3(5), ir1(9), ir2(4), ir3(5)
      dimension ia(7502), ib1(7502), ib2(7502)
      data jy / 0, 1, 2, 10, 11, 12, 20, 21, 22 /
      data jz / '00', '01', '02', '10', '11', '12', '20', '21', '22' /
      data jw / 0, 4, 5, 6 /
      data sex / 'N AUTO', 'SOMAL ', ' SEX-L', 'INKED ' /
      data lend /0/
c     call getarg(1,profile)
c     call getarg(2,terfile)
c     call getarg(3,isex)
      open(unit=48,file='INTERMED.1',status='OLD',form='UNFORMATTED') 
      open(unit=34,file='INTERMED.2',status='UNKNOWN',form=
     &'UNFORMATTED') 
      open(unit=35,file='INTERMED.3',status='UNKNOWN',form=
     &'UNFORMATTED') 
      open(unit=52,file='INTERMED.ALP',status='OLD') 
      open(unit=61,file=profile,status='OLD',access='append')
cc JH Zhao 6/6/2001 from fileopt='eof'
      open(unit=62,file=terfile,status='OLD',access='append')
cc JH Zhao 6/6/2001 from fileopt='eof'
 
c     ijp, ijz and ijw have been added to cope with a printing problem
c     that occurs when the program is informing the user of what selection
c     rejection parameters have been selected.

      read(unit=52, fmt=377) i, idebug
  377 format(2i6)
      do 371 i = 1, 9
      ijz(i) = 0
      is1(i) = '   '
  371 ir1(i) = '   '
      do 372 i = 1, 4
      is2(i) = '   '
      ijw(i) = 0
  372 ir2(i) = '   '
      do 373 i = 1, 5
      ijp(i) =0
      is3(i) = '   '
  373 ir3(i) = '   '
      is4 = 'ALL'
      ir4 = 'NIL'
c                                                                       
c         
c     ASK SELECTION/REJECTION QUESTIONS                                 
c         
c                                                                       
c         
c---- INTERMEDIATE FILE AUTOSOMAL(A) OR SEX-LINKED(S)?                  
c         
      ir4 = 'NIL'
c FOR VAX                                                               
c         
      jsex = 1
c      This question is asked by the shell script and passed to this 
c      routine as a parameter.
c      write(unit=*, fmt=3000) 
c 3000 format(31h AUTOSOMAL(A) OR SEX-LINKED(S)?)
c      read(unit=*, fmt=1) isex
      if ((isex .ne. 'A') .and. (isex .ne. 'S').and.(isex.ne.'a')
     &     .and.(isex.ne.'s')) goto 15
      if ((isex .eq. 'S').or.(isex.eq.'s')) jsex = 2
      inin = 114
      if (jsex .eq. 2) inin = 136
      rewind(unit=48) 
  369 read(unit=48) inulls, (ia(i),i = 1, inin)
      if (inulls .eq. nulls) goto 370
      goto 369
  370 rewind(unit=48) 
      if (idebug .eq. 1) then
      do 374 i = 1, 6
  374 write(unit=61, fmt=376) (ia(j),j = 1 + (12 * (i - 1)), 12 * i)
  376 format(1h ,12i10)
      end if
      nblock = ia(1)
      itypaf = ia(26)
      if (jsex .eq. 2) itypaf = ia(46)
      itypgr = ia(54)
      if (jsex .eq. 2) itypgr = ia(74)
      write(unit=61, fmt=375) itypaf
  375 format(8h ITYPAF=,i9)
      write(unit=*, fmt=3002) 
c                                                                       
c         
 3002 format(51h MATING TYPE SELECTION(00,01,02,10,11,12,20,21,22)?)
      call freem(fm1, nsel1, ival, ierr)
      if (ierr .eq. 1) goto 500
      if (nsel1 .gt. 9) goto 17
      do 40 i = 1, nsel1
      jj = ival(i) / 10
      kk = ival(i) - (jj * 10)
      if ((jj .gt. 2) .or. (kk .gt. 2)) goto 15
  
   40 isel1(i) = ival(i)
      write(unit=*, fmt=3003)      
 3003 format(28h POINTER SELECTION(0,4,5,6)?)
      call freem(fm2, nsel2, ival, ierr)
      if (ierr .eq. 1) goto 500
      if (nsel2 .gt. 4) goto 17
      do 41 i = 1, nsel2
      if ((((ival(i) .ne. 0) .and. (ival(i) .ne. 4)) .and. (ival(i)
     & .ne. 5)) .and. (ival(i) .ne. 6)) goto 15     
   41 isel2(i) = ival(i)
      if((isex .eq. 'A').or.(isex.eq.'a')) write(unit=*, fmt=3004) 
      if((isex .eq. 'S').or.(isex.eq.'s')) write(unit=*, fmt=3005) 
 3004 format(37h POINTER DEGREE SELECTION(1,2,3,4,5)?)         
 3005 format(35h POINTER DEGREE SELECTION(1,2,3,4)?)
      call freem(fm3, nsel3, ival, ierr)
      if (ierr .eq. 1) goto 500
      if (((isex .eq. 'A') .and. (nsel3 .gt. 5)) .or. ((isex .eq. 'S')
     & .and. (nsel3 .gt. 4))) goto 17
      do 42 i = 1, nsel3
      if (((isex .eq. 'A') .and. (ival(i) .gt. 5)) .or. ((isex .eq. 'S')
     & .and. (ival(i) .gt. 4))) goto 15      
   42 isel3(i) = ival(i)
      write(unit=*, fmt=3001) 
 3001 format(53h ASCERTAINMENT COMPLETE(C) OR MULTIPLE INCOMPLETE(M)?)
c     read(unit=*, fmt=1) isel4
      read(fm4,1) isel4
    1 format(a1)
       if(isel4.eq.'c')isel4='C'
       if(isel4.eq.'m')isel4='M'
       if (((isel4 .ne. '   ') .and. (isel4 .ne. 'C')) .and. (isel4 .ne. 
     &'M')) goto 15

c     Create the output information on which Mating types have been 
c     selected.
      if ((nsel1 .ne. 0) .and. (nsel1 .ne. 9))then
        i1 = 0
        i2 = 0

c	  Selection list for mating type.
          do 430 i = 1, 9
            do 431 j = 1, nsel1
              if(jy(i) .eq. isel1(j))then
                 i1 = i1 + 1
                 is1(i1) = jz(i)
		 ijz(i) = 1
		 goto 430
              end if
  431       continue
  430     continue

c         Rejection list for mating type.
          do 432 i=1,9
            if(ijz(i).eq.0)then
               i2 = i2 + 1
               ir1(i2) = jz(i)
            end if
  432      continue
      else
        is1(1) = 'ALL'
        ir1(1) = 'NIL'
      end if

c     Create the output information on which Pointer Selection values
c     have been requested.
      if ((nsel2 .ne. 0) .and. (nsel2 .ne. 4))then
        i1 = 0
        i2 = 0

c         Selection for pointers
          do 435 i = 1, 4
            do 406 j = 1, nsel2
              if (jw(i) .eq. isel2(j))then
                i1 = i1 + 1
                is2(i1) = char(jw(i) + 48)
		ijw(i) = 1
                goto 435
              end if
  406      continue
  435    continue

c        Rejection for pointers.
         do 407 i=1,4
              if(ijw(i).eq.0)then
                i2 = i2 + 1
                ir2(i2) = char(jw(i) + 48)
              end if
  407       continue
      else
        is2(1) = 'ALL'
        ir2(1) = 'NIL'
      end if

c     Create the output information on which Pointer Degree Selection
c     values are requested.
      nnn = 5
      if (isex .eq. 'S  ') nnn = 4

      if ((nsel3 .ne. 0) .and. (nsel3 .ne. nnn))then
        i1 = 0
        i2 = 0

c        Selection for pointer degree.
 
         do 410 i = 1, nnn
            do 414 j = 1, nsel3
              if (i .eq. isel3(j))then
                 i1 = i1 + 1
                 is3(i1) = char(i + 48)
		 ijp(i) = 1
                 goto 410
              end if
  414       continue
  410    continue

c        Rejection for pointer degree.

	 do 419 i=1,nnn
	     if(ijp(i).eq.0)then
                 i2 = i2 + 1
                 ir3(i2) = char(i + 48)
             end if
419      continue
      else
         is3(1) = 'ALL'
         ir3(1) = 'NIL'
      end if

c     Split C or M as needs be.
  394 if (isel4 .eq. '   ') goto 395
      is4 = isel4
      if (is4 .eq. 'C  ') ir4 = 'M  '
      if (is4 .eq. 'M  ') ir4 = 'C  '

c     Report the values that will split the records.
  395 write(unit=62, fmt=420) sex(1,jsex), sex(2,jsex)
  420 format(///1x,130(1h=)/1x,31hEMX SELECTION/REJECTION PROGRAM///1x,
     &13hTHIS IS FOR A,2a6,4hFILE)
      read(unit=52, fmt=3) title
    3 format(a12)
      write(unit=62, fmt=2) title
    2 format(/6h FILE ,7a12)
      write(unit=62, fmt=421) 
  421 format(//5x,28hSELECTION/REJECTION CRITERIA//5x,15hSELECTION FILE:
     &/)
      write(unit=62, fmt=422)(is1(ij),ij=1,9),(is2(ij),ij=1,4),
     &     (is3(ij),ij=1,5),is4
  422 format(10x,15hMATING TYPE:   ,9(2x,a3)//10x,15hPOINTER:       
     &,4(2x,a3)//10x,15hPOINTER DEGREE:,5(2x,a3)//10x,15hASCERTAINMENT: 
     &,2x,a3/)
      write(unit=62, fmt=423) 
  423 format(5x,15hREJECTION FILE:/)
      write(unit=62, fmt=422) (ir1(ij),ij=1,9),(ir2(ij),ij=1,4),
     &     (ir3(ij),ij=1,5),ir4
      write(unit=61, fmt=420) sex(1,jsex), sex(2,jsex)
      write(unit=61, fmt=2) title
      write(unit=61, fmt=421) 
      write(unit=61, fmt=422)(is1(ij),ij=1,9),(is2(ij),ij=1,9)
     &     ,(is3(ij),ij=1,5),is4
      write(unit=61, fmt=423) 
c                                                                       
c         
c     START PROCESSING                                                  
c         
c                                                                       
c         
      write(unit=61, fmt=422)(ir1(ij),ij=1,9),(ir2(ij),ij=1,4)
     &     ,(ir3(ij),ij=1,5),ir4
      do 5 i = 1, 7502
      ib1(i) = nulls
    5 ib2(i) = nulls
      nfam1 = 0
      nfam2 = 0
      nfb1 = 0
      nfb2 = 0
      nout1 = 0
      nout2 = 0
      icount = 0
      ilast = 0
    9 if (nblock .eq. 0) goto 350
      nblock = nblock - 1
      read(unit=48, end=350) (ia(i),i = 1, 7502)
      if (ia(1) .eq. nulls) goto 350
   12 nfamb = ia(7501)
      nout = 0
c mask the following loop
c     do 240 ii = 1, nfamb
c to avoid jump into a for block
c JH Zhao 16/4/2005
      ii = 1
999   continue
      icount = icount + 1
      nr = ia(nout + 2)
      lend = 18 + (6 * nr)
      if (jsex .eq. 2) lend = 18 + (8 * nr)
      do 201 i = 1, lend
c     DETERMINE IF FAMILY IS TO BE SELECTED OR REJECTED                 
c         
c     ----------------------------------------------------              
c         
  201 itemp(i) = ia(nout + i)
      if ((nsel1 .eq. 0) .or. (nsel1 .eq. 9)) goto 251
      do 250 i = 1, nsel1
      if (itemp(1) .eq. isel1(i)) goto 251
  250 continue
      goto 209
  251 if ((nsel2 .eq. 0) .or. (nsel2 .eq. 4)) goto 257
      do 252 i = 1, nsel2
      if (isel2(i) - 5) 253, 255, 256
  253 if (isel2(i) .eq. 4) goto 254
      if (((itemp(5) + itemp(6)) + itemp(7)) .eq. 0) goto 257
      goto 252
  254 if (itemp(5) .ne. 0) goto 257
      goto 252
  255 if (itemp(6) .ne. 0) goto 257
      goto 252
  256 if (itemp(7) .ne. 0) goto 257
  252 continue
      goto 209
  257 if ((nsel3 .eq. 0) .or. (nsel3 .eq. nnn)) goto 267
      do 260 i = 1, 3
      iporel = itemp(i + 4)
      if ((isex .eq. 'A').or.(isex.eq.'a')) then
      if (i .eq. 3) iporel = iporel - 35
      iporel = (iporel / 7) + 1
      else
      iporel = iporel - ((i - 1) * 28)
      iporel = (iporel / 7) + 1
      end if
      do 260 j = 1, nsel3
      if (iporel .eq. isel3(j)) goto 267
  260 continue
      goto 209
  267 if (isel4 .eq. '   ') goto 203
c     FAMILY IS SELECTED                                                
c         
      if (((isel4 .eq. 'C') .and. (itemp(14) .ne. 0)) .or. ((isel4 .eq. 
     &'M') .and. (itemp(14) .le. 0))) goto 209
  203 if ((nout1 + lend) .gt. 7500) goto 300
  204 nfam1 = nfam1 + 1
      nfb1 = nfb1 + 1
      do 205 i = 1, lend
      nout1 = nout1 + 1
  205 ib1(nout1) = itemp(i)
      select = .true.
c     FAMILY IS REJECTED                                                
c         
      goto 200
  209 if ((nout2 + lend) .gt. 7500) goto 400
      nfam2 = nfam2 + 1
      nfb2 = nfb2 + 1
      do 210 i = 1, lend
      nout2 = nout2 + 1
  210 ib2(nout2) = itemp(i)
      select = .false.
  200 nout = nout + lend
      if ((itypaf .eq. 1) .or. (itypgr .eq. 1)) then
      if (select) then
      read(unit=52, fmt=4) cfam(1,nfam1), cgr(1,nfam1)
      else
      read(unit=52, fmt=4) cfam(2,nfam2), cgr(2,nfam2)
      end if
      end if
    4 format(2a6)
c 240 continue
      ii = ii + 1
240   if (ii <= nfamb) goto 999 
      if (ia(nout) .eq. nulls) goto 190
      if (ia(nout + 1) .eq. nulls) goto 9
      if (nout .eq. 7500) goto 9
  190 write(unit=61, fmt=191) nout, nfamb, ii
  191 format(1x,28hSOMETHING IS WRONG WITH NOUT,3i5)
      goto 500
  300 ib1(7501) = nfb1
      write(unit=34) (ib1(i),i = 1, 7502)
  302 do 305 i = 1, 7502
  305 ib1(i) = nulls
      nfb1 = 0
      nout1 = 0
      if (ilast) 400, 204, 400
  400 ib2(7501) = nfb2
      write(unit=35) (ib2(i),i = 1, 7502)
  402 do 405 i = 1, 7502
  405 ib2(i) = nulls
      nfb2 = 0
      nout2 = 0
      if (ilast) 351, 209, 351
  350 ilast = 1
      goto 300
  351 rewind(unit=48) 
  352 read(unit=48) inulls, (ia(i),i = 1, inin)
      if (inulls .ne. nulls) goto 352
      write(unit=34) inulls, (ia(i),i = 1, inin), (0,i = 1, 7501 - inin)
      write(unit=35) inulls, (ia(i),i = 1, inin), (0,i = 1, 7501 - inin)
      if ((itypaf .eq. 1) .or. (itypgr .eq. 1)) then
      rewind(unit=52) 
      write(unit=52, fmt=377) 2, nfam1
      write(unit=52, fmt=3) title
      do 360 i = 1, nfam1
  360 write(unit=52, fmt=4) cfam(1,i), cgr(1,i)
      do 365 i = 1, nfam2
  365 write(unit=52, fmt=4) cfam(2,i), cgr(2,i)
      end if
      write(unit=62, fmt=411) icount, nfam1, nfam2
  411 format(//5x,36hNO. OF FAMILIES IN     INPUT FILE = ,i5//5x,
     &36hNO. OF FAMILIES IN SELECTION FILE = ,i5//5x,
     &36hNO. OF FAMILIES IN REJECTION FILE = ,i5)
      write(unit=61, fmt=411) icount, nfam1, nfam2
      goto 666
   15 write(unit=61, fmt=16) 
   16 format(1x,25hINPUT VALUES OUT OF RANGE)
      goto 500
   17 write(unit=61, fmt=18) 
   18 format(1x,27hTOO MANY SELECTION CRITERIA)
      goto 500
  666 if (idebug .ne. 1) goto 500
      rewind(unit=34) 
      rewind(unit=35) 
      if (itypaf .eq. 1) then
      rewind(unit=52) 
      read(unit=52, fmt=377) i, i
      read(unit=52, fmt=3) title
      end if
      do 700 ii = 1, 2
      inp = 33 + ii
      nfam = 0
  709 read(unit=inp, end=700) (ia(i),i = 1, 7502)
      if (ia(1) .eq. nulls) goto 700
  702 nfamb = ia(7501)
      nout = 0
      do 720 jj = 1, nfamb
      nfam = nfam + 1
      nr = ia(nout + 2)
      lend = 18 + (6 * nr)
      if (jsex .eq. 2) lend = 18 + (8 * nr)
      nout1 = nout + 1
      nout2 = nout + 14
      nout3 = (nout / 2) + 8
      nout4 = (nout / 2) + 9
      if (itypaf - 1) 741, 741, 742
  742 write(unit=61, fmt=743) (ia(i),i = nout1, nout2), a(nout3), nfam, 
     &nout1
  743 format(1x,8hOUTPUT: ,12i5,i10,i5,5x,f6.0,3x,7hNFAM = ,i5,4x,
     &10h1ST ADDR.=,i4)
      goto 745
  741 read(unit=52, fmt=4) famnum, gr
      write(unit=61, fmt=744) (ia(i),i = nout1, nout2), famnum, nfam, 
     &nout1
  744 format(1x,8hOUTPUT: ,12i5,i10,i5,5x,a6,3x,7hNFAM = ,i5,4x,
     &10h1ST ADDR.=,i4)
  745 j = nout4
      do 746 i = 1, nr
      if (isex .eq. 'A') then
      k = (2 * j) + 3
      write(unit=61, fmt=747) a(j + 1), ia(k), ia(k + 1)
  747 format(6x,f12.6,3i7)
      j = j + 3
      else
      k = (2 * j) + 4
      write(unit=61, fmt=747) a(j + 1), ia(k), ia(k + 1), ia(k + 2)
      j = j + 4
      end if
  746 continue
  720 nout = nout + lend
      goto 709
  700 continue
  500 endfile(unit=34) 
      endfile(unit=35) 
      endfile(unit=52) 
      write(unit=62, fmt=502) 
      write(unit=61, fmt=502) 
  502 format(1x,130(1h=))
      close(unit=34) 
      close(unit=35) 
      close(unit=52) 
      close(unit=61) 
      close(unit=62) 
      stop 
      end
      subroutine freem(fm, n, ival, ierr)
      implicit double precision (o-z, a-h)
c     IZERO IS '0' IN THE HARRIS COLLATING SEQUENCE                     
c         
c                                                                       
c         
      parameter (izero = 48)
      dimension a(81), ival(30), itemp(8)
c FOR VAX                                                               
c         
      character fm*81, a*1
      ncol = 0
c                                                                       
c         
c FOR HARRIS                                                            
c         
c     READ(33,1)A                                                       
c         
c                                                                       
c         
# 9  "freem.for"
      read(unit=fm, fmt=1) a
# 14 "freem.for"
    1 format(81a1)
      do 10 i = 1, 30
   10 ival(i) = 0
      ierr = 0
      n = 0
      jack = 1
      do 100 icol = 1, 81
      if ((a(icol) .eq. ' ') .or. (a(icol) .eq. ',')) goto 90
      ix = ichar(a(icol)) - izero
      if ((ix .ge. 0) .and. (ix .le. 9)) goto 110
      write(unit=61, fmt=19) 
   19 format(1x,13hERROR IN DATA)
      ierr = 1
      goto 900
  110 if (jack .eq. 0) goto 113
  111 n = n + 1
      jack = 0
      do 112 i = 1, 8
  112 itemp(i) = 0
  113 ncol = ncol + 1
      itemp(ncol) = ix
      goto 100
   90 if (jack .eq. 1) goto 100
      jack = 1
      isum = 0
      kval = 1
      do 80 i = 1, ncol
      jval = kval * itemp((ncol + 1) - i)
      isum = isum + jval
   80 kval = kval * 10
      ival(n) = isum
  100 continue
  900 return 
      end

!MANUAL TO BE UPDATED
    
!Generates initial (unrelaxed) FCC atomic coordinates for half of torus (HT) + wetting layer (WL) (HTWL) in a matrix box.
!Sorts atoms and prepares the initial atom configuration.
!Subroutine(SR) FCCbox generates FCC structure (files 11->At1 & 12->At2)
    !in a box with 8*na=8*(nxmin+nxmax+1)*(nymin+nymax+1)*(nzmin+nzmax+1) atoms (type 1 & 2) and latt_ct=1;
    !there are two kinds of atoms, Atom1 and Atom2. At11 & At21 rearrange arrays At1 & At2 in C style.
!With At11 and At21 two boxes are built, one with latt_ct=a1 (for OUT material, GaAs) 
!   and the other one with latt_ct=a2 (for IN material, GaSb):
!xIN1, yIN1, zIN1 are coords of Atom1 (Ga, inside) in box2; xOUT1, yOUT1, zOUT1 are coords of Atom1 (Ga, outside) in box1.
!xIN2, yIN2, zIN2 are coords of Atom2 (Sb, inside) in box2; xOUT3, yOUT3, zOUT3 are coords of Atom3 (As, outside) in box1.
!Any other orientation of the BULK can EASILY be obtained from this code with translation+rotations. 
    !*********
!SR QDHTF generates the geometrical function for HTWL and interogates if atoms of box 1 are OUT HTWL & atoms of box 2 are IN HTWL.
!SR HTWLB locates (coords) the atoms of: Ga in HTWL in box (HTWLB) in file 140;
!   Sb IN HTWL in file 160; As OUT HTWL in file 180.
!*********
!SR Cross_section  locates atoms in a 'cross-sections' of thickness=layer for HTWLB
    !layer = cross-section layer thickness
!   Below files collect coordinates for HTWLB atoms in a layer centered in origin, parallel with yz-plan, and of thicknes=layer:
!   13 (GayzIN.dat) collects IN Atom1 (Ga) coords; 14 (Sbyz.dat) collects IN Atom2 (Sb) coords
!   15 (GayzOUT.dat) collects OUT Atom1 (Ga) coords; 16 (Asyz.dat) collects Atom3 (OUT, As) coords
!**************************************************
!The loop 'Do n=1, nt' ends by counting #atoms of type 1,2,3, i.e. Nat1,Nat2,Nat3 & Nat23 (for Nat2+Nat3) in HTWLB.
!*********
!IniConfig uses file 140,160,180 and generates initial atom coords of HTWLB in file 200.
!*********
!SR Sortx123(x) finds distances and indexes (ordering #) of all atoms and their first 4 neighbors by use of file 200.
    !It needs SR Sorting4X which finds distances to closes 4 neighbors
    !and indexes (order #) of all atoms of given type and their first 4 neighbors
    !Thus for
  !Ga:in files 240(Tip23c-distances & neighbor indexes permutation(NIP)),241(Nb1-NIP and 4 initial indexes of closest 4 neighbors)
  !Sb:in files 260(Tip21c----------------II----------------------------),261(Nb1-----------------------II------------------------)
  !As:in files 280(Tip31c----------------II----------------------------),281(Nb1-----------------------II------------------------)
!*********
!SR EnConfig(x) rearrange data in files 240,...281 and generates 300 to be used in U_DU.   
!**************************************************
!TO DO:
!   to be implement in L-BFGS
!**********To Be UPDATED!**************
   
    Program Main
    use DataType !necessary for a1,a2,Rt,Rq,ratio,Nat1,..., Nat23 type of data and storing their values
!---------------------------------------------------------------------------------------------
    Real*8:: xIN1, yIN1, zIN1, xOUT1, yOUT1, zOUT1 !coords: xIN1,...Ga IN; xOUT1,...Ga IN
    Real*8:: xIN2, yIN2, zIN2, xOUT3, yOUT3, zOUT3 !coords: xIN2,...Sb IN;  xOUT3,...As OUT
    Real*8:: layer,WL
    Integer:: Nw, nxmin, nxmax, nymin, nymax, nzmin, nzmax, na, nt, n
    Logical:: logIN1, logOUT1, logIN2, logOUT3
    Real*8, allocatable:: At1(:,:),At2(:,:) , At11(:,:),At21(:,:)
    Integer::i1,i2,i3
    Integer::QD
!---------------------------------------------------------------------------------------------
    Open (15,file='Parameters.inp',status='old') !reading input parameters 
    Read (15,*) a1, a2
    Read (15,*) ratio    
    Close (15)
 Print*,'Type 1 for cone QD or any integer for hemi-torus QD'
    Read(*,*) QD
    if (QD==1)then
     Rc=30; h=30
     Nw=2;nxmin=9;nxmax=9;nymin=9;nymax=9;nzmin=3;nzmax=7;     
    Print*,'One considers a cone QD. Parameters may be changed, see Readme'  
    else
     Rt=30; Rq=13
     Nw=2;nxmin=12;nxmax=12;nymin=12;nymax=12;nzmin=2;nzmax=3; 
     Print*,'One considers hemi-torus QD. Parameters may be changed, see Readme'  
    end  if
continue
!*********************************************
!following na, nt are # atoms in a FCC box

layer=a2*ratio; ! cross-section layer thickness in plane yz. 
na=(nxmin+nxmax+1)*(nymin+nymax+1)*(nzmin+nzmax+1) !half of the total # atoms of one kind;
            !generally, there are 4 kinds of atoms: 2 kinds for IN and 2 kinds for OUT HT;
            !one atom may be common for IN and OUT (Ga, here).    
nt=4*na !total # Atom1 in box: some are IN the rest OUT HT.
            !8*na total # atoms (Atom1 & Atom2) in box: some are IN the rest OUT HT.
WL=Nw*a2/2  !wetting layer thickness; Nw # MLs.    
!*********************************************
    
!FCC atom1 & atom2 coords in a box with latt_ct=1.    
    !SR FCCbox gives FCC atom coords of box ((-nxmin,nxmax),(-nymin,nymax),(-nzmin,nzmax)) with nt atoms and latt_ct=1.
Call FCCbox(nxmin,nxmax,nymin,nymax,nzmin,nzmax) !generates files 11 & 12 
    Allocate (At1(na,20),At2(na,20),At11(nt,5),At21(nt,5)) 
    Rewind(11)            !file 11 collects sets: type of atom(1), index(ordering), FCC bulk atom coordinates and latt_ct=1.
    Read (11,*) At1
        !If stack overflow: Project>Torus_Coord Properties...>Fortran>Optimization>Heap Arrays>0
    At11=RESHAPE(At1, (/nt,5/), ORDER = (/2,1/)) !Rearrange array At1 in C style:
        !At11(n,1)=n(index), At11(n,2)=1(Atom type), At11(n,3)=x1(n), At11(n,4)=y1(n), At11(n,5)=z1(n) for latt_ct=1
    
    Rewind(12)            !file 12 collects sets: type of atom(2), index(ordering), FCC bulk atom coordinates.
    Read (12,*) At2
    At21=RESHAPE(At2, (/nt,5/), ORDER = (/2,1/)) !Rearrange array At2 in C style:
        !At21(n,1)=n(index), At21(n,2)=2(Atom type), At21(n,3)=x2(n), At21(n,4)=y2(n), At21(n,5)=z2(n) for latt_ct=1
    continue
!*********************************************
!I.2 FCC atom coordinates of HTWL in box
!At11 & At21 (which are built to get FCC with latt_ct=1) are used to build FCC bulk atom coordinates in 2 boxes, one is obtained from 
!   box1 with latt_ct=a1 (GaSb), the other is obtained box2 with latt_ct=a2 (GaAs); e.g., xIN1,.... 
!   In the below files 38,39,140.., 180 & 13,..., 16 the atoms are arranged IN and OUT FTWL-matrix and may be used for GRAPHs.
!   110 ,...310 would be used for SORTING.
    !--------------------------------------------- 

i1=0; i2=0; i3=0;
!--------------------------------------------------------------------------------
!This section sets atom coords FCC in HTWL in box ((-nxmin,nxmax),(-nymin,nymax),(-nzmin,nzmax))*a1 centered at (0,0,0).
Do n=1, nt  !n counts atoms (type 1 (e.g., Ga) or 2 (e.g.,As or Sb))
    !First, one sets atom coordinates FCC symmetry in 2 boxes with latt_ct = a1,a2, respectively:
    
    !box2 FCC filled with common Atom1 (Ga); latt_ct=a2 (IN):  
     xIN1=a2*At11(n,3); yIN1=a2*At11(n,4); zIN1=a2*At11(n,5)  !x,y,z positions
    
    !box1 FCC filled with common Atom1 (Ga), latt_ct=a1 (OUT): 
     !§xOUT1=xIN1/a2*a1; yOUT1=yIN1/a2*a1; zOUT1=zIN1/a2*a1 !§or
     xOUT1=a1*At11(n,3); yOUT1=a1*At11(n,4); zOUT1=a1*At11(n,5)
     
    !box2 FCC filled with Atom2 (Sb); latt_ct=a2 (IN):    
     xIN2=a2*At21(n,3); yIN2=a2*At21(n,4); zIN2=a2*At21(n,5)
       
    !box1 FCC filled with Atom3 (As); latt_ct=a1 (OUT):
     !§xOUT3=xIN2/a2*a1; yOUT3=yIN2/a2*a1; zOUT3=zIN2/a2*a1 !§or
      xOUT3=a1*At21(n,3); yOUT3=a1*At21(n,4); zOUT3=a1*At21(n,5)
      
    !Second, one distruíbutes atoms IN and OUT HT (a1<a2):
   if (QD==1) then     
    Call ConeQD(xIN1,yIN1,zIN1,WL,logIN1)      !interogates if xIN1,... are IN HTWL- for Atom1(Ga)  
    Call ConeQD(xOUT1,yOUT1,zOUT1,WL,logOUT1)      !interogates if xIN1,... are IN HTWL- for Atom1(Ga)
    Call ConeQD(xIN2,yIN2,zIN2,WL,logIN2)      !interogates if xIN2,... are IN HTWL - for Atom2(Sb)
    Call ConeQD(xOUT3,yOUT3,zOUT3,WL,logOUT3)  !interogates if xOUT3,... are OUT HTWL - for Atom3(As)
   else   
    Call QDHTF(xIN1,yIN1,zIN1,WL,logIN1)      !interogates if xIN1,... are IN HTWL- for Atom1(Ga)
    Call QDHTF(xOUT1,yOUT1,zOUT1,WL,logOUT1)  !interogates if xOUT1,... are OUT HTWL - for Atom1(Ga)
    Call QDHTF(xIN2,yIN2,zIN2,WL,logIN2)      !interogates if xIN2,... are IN HTWL - for Atom2(Sb)
    Call QDHTF(xOUT3,yOUT3,zOUT3,WL,logOUT3)  !interogates if xOUT3,... are OUT HTWL - for Atom3(As)
    end if
    
    !Third, one collects FCC atom coords in HTWLB:
   ! Call HTWLB(xIN1,yIN1,zIN1,xOUT1,yOUT1,zOUT1,xIN2,yIN2,zIN2,xOUT3,yOUT3,zOUT3, &
    !           logIN1,logOUT1,logIN2,logOUT3,nxmin,nxmax,nymin,nymax,i1,i2,i3) !atom coordinates for IN and OUT HTWL in box;
                                                             !i1,i2,i3 are loop counters for kind of atoms (Ga,Sb,As resppectively);
                                                             !writes files 140,160,180 which store unrelaxed atom coords.      
    
    Call HTWLB(xIN1,yIN1,zIN1,xOUT1,yOUT1,zOUT1,xIN2,yIN2,zIN2,xOUT3,yOUT3,zOUT3, &
               logIN1,logOUT1,logIN2,logOUT3,nxmin,nxmax,nymin,nymax,i1,i2,i3) !atom coordinates for IN and OUT HTWL in box;
                                                             !i1,i2,i3 are loop counters for kind of atoms (Ga,Sb,As resppectively);
                                                             !writes files 140,160,180 which store unrelaxed atom coords.      
       
    Call Cross_section(xIN1, yIN1, zIN1, xOUT1, yOUT1, zOUT1,xIN2, yIN2, zIN2, xOUT3, &
                      yOUT3,zOUT3,logIN1,logOUT1,logIN2,logOUT3,nxmin,nxmax,nymin,nymax,layer) !atom coords in 'cross-section' for HTWL in box
End Do
!------------------------------
    Nat1=i1 !# atom1(Ga)
    Nat2=i2 !# atom2(Sb)
    Nat3=i3 !# atom3(As)
    Nat23=Nat2+Nat3 !Nat1,..., Nat23 type is stored in module DataType
!Obs.: (Nat1+Nat2+Nat3) .ne. 8*na; some of 8*na atoms are out of the box resulting by superposition of the 2 boxes with latt_ct = a1, a2. 
!------------------------------
Print*,''
Print*,'Unrelaxed configuration for GaSb/GaAs QD done.'
Print*,'The atoms positions are stored in files:'
Print*,'    140(Ga),160(Sb),180(As) for QD volume.'
Print*,'    13(Ga-in),14(Sb),15(Ga-out),16(As) for QD "cross-section".'
!READ(*,*)
    End Program Main    
!=======================================================
 
    Subroutine FCCbox(nxmin,nxmax,nymin,nymax,nzmin,nzmax)
    ! Generates FCC atom positions in a box with latt_ct=a=1; for each set (i,j,k) there are 8 atoms posed in the box.
    ! One adds atoms in layer types (LTs) as follows. With LT1 as the plane xy having z=0 one has (for k>0):
    !       LT3 is upwardly shifted by k*a/2 with respect to LT1; LT1 and LT3 are filled with Atom1.
    !       LT2 is upwardly shifted by k*3*a/4 with respect to LT1.
    !       LT4 is upwardly shifted by k*a/4 with respect to LT1; LT2 and LT4 are filled with Atom2.
    !       In addition, LT2, LT3, LT4 are horizontally shifthed with respect to LT1, see ZincBlende-structure-planes.nb.
    ! Any other orientation of the BULK atoms can EASILY be obtained with this subroutine by implementing translation+rotations.   
    ! Atom1-Atom2; ex: Atom1=Ga, Atom2=Sb or As for BULK
       
    Integer::nxmin,nxmax,nymin,nymax,nzmin,nzmax,i,j,k, n, m ! n is index
   
    n=0; m=0;
    Do 1 i=-nxmin,nxmax
        Do 1 j=-nymin,nymax
                Do 1 k=-nzmin,nzmax
            n=n+1
             m=n*4-3
        
    Write(11,*) 1, m, i, j, k                      !Atom1 in LT1
    Write(11,*) 1, m+1, (2*i+1)/2d0, (2*j-1)/2d0, k  !Atom1 in LT1
    Write(11,*) 1, m+2, i, (2*j+1)/2d0, (2*k+1)/2d0  !Atom1 in LT3
    Write(11,*) 1, m+3, (2*i+1)/2d0, j, (2*k+1)/2d0  !Atom1 in LT3
    
    Write(12,*) 2, m, (4*i-1)/4d0, (4*j+1)/4d0, (4*k+3)/4d0 !Atom2 in LT2
    Write(12,*) 2, m+1, (4*i+1)/4d0, (4*j-1)/4d0, (4*k+3)/4d0 !Atom2 in LT2
    Write(12,*) 2, m+2, (4*i+1)/4d0, (4*j+1)/4d0, (4*k+1)/4d0 !Atom2 in LT4
    Write(12,*) 2, m+3, (4*i+3)/4d0, (4*j-1)/4d0, (4*k+1)/4d0 !Atom2 in LT4
    
1   continue 
    End Subroutine FCCbox
!*******************************************************************
    
    Subroutine QDHTF(x,y,z,WL,log)
    !QDHTF-quantum dot HT function
    !Generates HT function in Cartesian
    !Rq-torus radius; Rt1-internal radius of the HT; Rt2-outer radius of the HT
   use DataType
    Real*8:: Rt1, Rt2, x, y, z, rho, WL
    Logical:: log, log1, log2
    Rt1=Rt-Rq
    Rt2=Rt+Rq  
	rho=dsqrt(x**2+y**2)
        
    log1=((rho.ge.Rt1).and.(rho.le.Rt2).and.((rho-Rt)**2+z**2.le.Rq**2).and.(z.ge.0d0))! log1=true for (x,y,z) inside HT    
    log2=((z.le.0d0).and.(z.ge.-WL)) ! log2=true for (x,y,z) inside WL    
    log=log1.or.log2 ! log=true for (x,y,z) inside HTWL
!    continue
    End Subroutine QDHTF
!*********************************************  
    
    Subroutine ConeQD(x,y,z,WL,log)
    !QDHTF-quantum dot HT function
    !Generates HT function in Cartesian
    !Rc-con radius; h-cone height
   use DataType
      Real*8::  x, y, z, rho, WL
    Logical:: log, log1, log2
	rho=dsqrt(x**2+y**2)
   
    log1=(rho.le.(h-z)*Rc/h) .and.(z.le.h).and.(z.ge.0d0)! log1=true for (x,y,z) inside Cone    
    log2=((z.le.0d0).and.(z.ge.-WL)) ! log2=true for (x,y,z) inside WL    
    log=log1.or.log2 ! log=true for (x,y,z) inside Cone+WL
  !  continue
    End Subroutine ConeQD
!*********************************************   
  
Subroutine  HTWLB(xIN1,yIN1,zIN1,xOUT1,yOUT1,zOUT1,xIN2,yIN2,zIN2,xOUT3,yOUT3,zOUT3, &
                     logIN1,logOUT1,logIN2,logOUT3,nxmin,nxmax,nymin,nymax,i1,i2,i3) !generates FCC coordinates for IN and OUT HTWL in box
    use DataType
    Logical:: logIN1, logOUT1, logIN2, logOUT3
    Integer:: nxmin,nxmax,nymin,nymax,i1,i2,i3 !, nzmin, nzmax, na, nt, n
    Real*8:: xIN1, yIN1, zIN1, xOUT1, yOUT1, zOUT1 !xIN1,...Ga IN; xOUT1,...Ga IN
    Real*8:: xIN2, yIN2, zIN2, xOUT3, yOUT3, zOUT3 !xIN2,...Sb IN;  xOUT3,...As OUT
    
    If(logIN1) then !IN (Ga) HT;if xIN1,.. are IN HTWL inside xy square of sides ((nxmin+nxmax+1),(nymin+nymax+1)) they are stored in 38&140
        if ((dAbs(xIN1).le.nxmax*a1).and.(dAbs(yIN1).le.nymax*a1)) then !makes the boxes IN and OUT of same size in xy plane (a1<a2, usually).
        i1=i1+1;
        write(38,*) xIN1,yIN1,zIN1 !only IN Atom1 (Ga)
        write(140,*) xIN1,yIN1,zIN1
       ! write(110,*) 1,i1, xIN1,yIN1,zIN1
        end if
    End If
     
    If(logOUT1) then !OUT (Ga) HT
        continue       
    else !if xOUT1,.. are OUT HT they are stored in 140
        if ((dAbs(xOUT1).le.nxmax*a1).and.(dAbs(yOUT1).le.nymax*a1)) then !makes the boxes IN and OUT of same size in xy plane (a1<a2, usually).
        i1=i1+1
        write(39,*)  xOUT1,yOUT1,zOUT1 !only OUT Atom1 (Ga)
        write(140,*) xOUT1,yOUT1,zOUT1   
      !  write(110,*) 1, i1, xOUT1,yOUT1,zOUT1 
        end if
    End If
        
        If(logIN2) then !IN (Sb) HT
        if ((dAbs(xIN2).le.nxmax*a1).and.(dAbs(yIN2).le.nymax*a1)) then !makes the boxes IN and OUT of same size in xy plane (a1<a2, usually).
        i2=i2+1
        write(160,*) xIN2,yIN2,zIN2  
       ! write(210,*) 2, i2, xIN2,yIN2,zIN2  
        end if
        End If
               
    If(logOUT3) then !OUT (As) HT
        continue
    else 
        if ((dAbs(xOUT3).le.nxmax*a1).and.(dAbs(yOUT3).le.nymax*a1)) then !makes the boxes IN and OUT of same size in xy plane (a1<a2, usually).
        i3=i3+1
        write(180,*) xOUT3,yOUT3,zOUT3  
       ! write(310,*) 3,i3, xOUT3,yOUT3,zOUT3  
        end if
    End If
    End Subroutine HTWLB
!******************************************************
    
Subroutine Cross_section(xIN1, yIN1, zIN1, xOUT1, yOUT1, zOUT1,xIN2, yIN2, zIN2, xOUT3, &
    yOUT3, zOUT3,logIN1, logOUT1, logIN2, logOUT3,nxmin,nxmax,nymin,nymax,layer)
!Sets 'cross-sections' atom coordinates for HTWL in box:layer centered on origin,
    !parallel with yz-plan, and of thicknes=layer
use DataType
Logical:: logIN1, logOUT1, logIN2, logOUT3
Integer:: nxmin,nxmax,nymin,nymax !, nzmin, nzmax, na, nt, n
 Real*8:: xIN1, yIN1, zIN1, xOUT1, yOUT1, zOUT1 , layer     !xIN1,...Ga IN; xOUT1,...Ga IN
    Real*8:: xIN2, yIN2, zIN2, xOUT3, yOUT3, zOUT3 !, WL !xIN2,...Sb IN;  xOUT3,...As OUT
    continue
If (dAbs(xIN1).lt.layer) then !IN Atom1 (Ga) HT
    if ((dAbs(xIN1).le.nxmax*a1).and.(dAbs(yIN1).le.nymax*a1)) then !boxes of IN and OUT type made to have the same size in xy plane; a1<a2.
        if(logIN1) then 
        write(13,*) yIN1, zIN1 
        end if
    end if
    End If
    
    If (dAbs(xIN2).lt.layer) then !IN Atom2 (Sb) HT    
        if(logIN2) then 
        if ((dAbs(xIN2).le.nxmax*a1).and.(dAbs(yIN2).le.nymax*a1)) then !boxes of IN and OUT type made to have the same size in xy plane; a1<a2.
        write(14,*) yIN2, zIN2 
        end if
    end if
    End If
         
    If (dAbs(xOUT1).lt.layer) then !OUT  Atom1 (Ga) HT
    if(logOUT1) then
        continue
    else
        if ((dAbs(xOUT1).le.nxmax*a1).and.(dAbs(yOUT1).le.nymax*a1)) then !makes the boxes IN and OUT of same size in xy plane (a1<a2, usually).
        write(15,*) yOUT1,zOUT1 
        end if
    end if
    End If
    
    If (dAbs(xOUT3).lt.layer) then !OUT  Atom3 (As) HT
    if(logOUT3) then
        continue
    else
        if ((dAbs(xOUT3).le.nxmax*a1).and.(dAbs(yOUT3).le.nymax*a1)) then !makes the boxes IN and OUT of same size in xy plane (a1<a2, usually).
        write(16,*) yOUT3,zOUT3 
        end if
    end if
    End If
    End Subroutine Cross_section

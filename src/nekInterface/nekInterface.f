c-----------------------------------------------------------------------
c
c NEK5000 Interface
c
c-----------------------------------------------------------------------
      subroutine nekf_ptr(ptr,id,len)

      implicit none

      integer len 
      character*(len) id
      integer i8
      integer*8 ptr
      pointer(ptr,i8)

      include 'SIZE'
      include 'TOTAL'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nrs_vd/ glo_num((2**ldim)*lelt), ngv
      integer*8 glo_num, ngv

      COMMON /SCNRS/ SC_NRS(LX1*LY1*LZ1*LELT*7)
      real SC_NRS

      if (id .eq. 'nelv') then 
         ptr = loc(nelv)
      elseif (id .eq. 'nekcomm') then 
         ptr = loc(nekcomm)
      elseif (id .eq. 'istep') then 
         ptr = loc(istep)
      elseif (id .eq. 'param') then 
         ptr = loc(param(1))
      elseif (id .eq. 'nelt') then 
         ptr = loc(nelt)
      elseif (id .eq. 'nelv') then 
         ptr = loc(nelv)
      elseif (id .eq. 'ndim') then 
         ptr = loc(ndim)
      elseif (id .eq. 'nx1') then 
         ptr = loc(nx1)
      elseif (id .eq. 'cb_scnrs') then
         ptr = loc(sc_nrs(1)) 
      elseif (id .eq. 'glo_num') then
         ptr = loc(glo_num(1)) 
      elseif (id .eq. 'ngv') then
         ptr = loc(ngv) 
      elseif (id .eq. 'xc') then
         ptr = loc(XC(1,1)) 
      elseif (id .eq. 'yc') then
         ptr = loc(YC(1,1)) 
      elseif (id .eq. 'zc') then
         ptr = loc(ZC(1,1)) 
      elseif (id .eq. 'xm1') then
         ptr = loc(XM1) 
      elseif (id .eq. 'ym1') then
         ptr = loc(YM1) 
      elseif (id .eq. 'zm1') then
         ptr = loc(ZM1) 
      elseif (id .eq. 'unx') then
         ptr = loc(UNX(1,1,1,1)) 
      elseif (id .eq. 'uny') then
         ptr = loc(UNY(1,1,1,1)) 
      elseif (id .eq. 'unz') then
         ptr = loc(UNZ(1,1,1,1)) 
      elseif (id .eq. 'cbc') then
         ptr = loc(cbc(1,1,1)) 
      elseif (id .eq. 'eface1') then
         ptr = loc(eface1) 
      elseif (id .eq. 'eface') then
         ptr = loc(eface) 
      elseif (id .eq. 'icface') then
         ptr = loc(icface(1,1)) 
      elseif (id .eq. 'ifield') then
         ptr = loc(ifield) 
      elseif (id .eq. 'zgm1') then
         ptr = loc(zgm1(1,1)) 
      elseif (id .eq. 'wxm1') then
         ptr = loc(wxm1(1)) 
      elseif (id .eq. 'zgm2') then
         ptr = loc(zgm2(1,1)) 
      elseif (id .eq. 'wxm2') then
         ptr = loc(wxm2(1)) 
      elseif (id .eq. 'boundaryID') then
         ptr = loc(boundaryID(1,1)) 
      elseif (id .eq. 'vx') then
         ptr = loc(vx(1,1,1,1)) 
      elseif (id .eq. 'vy') then
         ptr = loc(vy(1,1,1,1)) 
      elseif (id .eq. 'vz') then
         ptr = loc(vz(1,1,1,1)) 
      elseif (id .eq. 'pr') then
         ptr = loc(pr(1,1,1,1)) 
      elseif (id .eq. 't') then
         ptr = loc(t(1,1,1,1,1)) 
      elseif (id .eq. 'time') then
         ptr = loc(time)
      else
         write(6,*) 'ERROR: nek_ptr cannot find ', id
         call exitt 
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_setup(comm_in,path_in,session_in)

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
c
      include 'OPCTR'
      include 'CTIMER'

      integer comm_in
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /rdump/ ntdump

      common /nrs_vd/ glo_num((2**ldim)*lelt), ngv
      integer*8 glo_num, ngv
      common /ivrtx/ vertex ((2**ldim),lelt)
      integer vertex

      character session_in*(*),path_in*(*)
      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest 
      logical ifbswap

      ! set word size for REAL
      wdsize = sizeof(rtest)
      ! set word size for INTEGER
      isize = sizeof(itest)
      ! set word size for INTEGER*8
      isize8 = sizeof(itest8) 
      ! set word size for LOGICAL
      lsize = sizeof(ltest) 
      ! set word size for CHARACTER
      csize = sizeof(ctest)

      call setupcomm(comm_in,newcomm,newcommg,path_in,session_in)
      call iniproc()

      etimes = dnekclock()
      istep  = 0

      call initdim         ! Initialize / set default values.
      call initdat
      call files

      call setDefaultParam
      param(1)  = 1.0
      param(2)  = 1.0
      param(27) = 1  ! 1st-order in time
      param(32) = 1  ! read only vel BC from re2
      param(99) = -1 ! no dealiasing

      ifflow = .true.
      iftran = .true.
      ifheat = .false.
      ifvo   = .true.
      ifpo   = .true.

      call bcastParam

      call usrdat0

      call read_re2_hdr(ifbswap)
      call chkParam
      call mapelpr 
      call read_re2_data(ifbswap)
      do iel = 1,nelt
      do ifc = 1,2*ndim
         boundaryID(ifc,iel) = bc(5,ifc,iel,1)
      enddo
      enddo

      call setvar          ! Initialize most variables

      igeom = 2
      call setup_topo      ! Setup domain topology  
      call genwz           ! Compute GLL points, weights, etc.

      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat' 

      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if(nio.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2' 

      call fix_geom
      call geom_reset(1)    ! recompute Jacobians, etc.

      call vrdsmsh          ! verify mesh topology
      call mesh_metrics     ! print some metrics

      call setlog(.false.)  ! Initalize logical flags

      call bcmask  ! Set BC masks for Dirichlet boundaries.

      call set_vert(glo_num,ngv,2,nelv,vertex,.false.)

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

      ntdump=0

      call setics

      call dofcnt

      tinit = dnekclock_sync() - etimes

      call flush(6)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_outfld()

      include 'SIZE'
      include 'TOTAL'

      common /scrcg/ pm1(lx1,ly1,lz1,lelv)

      call copy(pm1,pr,nx1*ny1*nz1*nelv)
      call outfld('   ')

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_restart(rfile,l, getu, getp)

      character*(l) rfile
      integer getu, getp

      include 'SIZE'
      include 'RESTART'
      include 'INPUT'

      call blank(initc(1),132)
      call chcopy(initc(1),rfile,l)
      call setics()

      getu = 0
      getp = 0
      if(ifgetu) getu = 1
      if(ifgetp) getp = 1

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_end()

      call nek_end()

      return
      end
c-----------------------------------------------------------------------
      real function nekf_cfl(u,v,w,dt)

      real u(*), v(*), w(*), dt
      real cfl

      call compute_cfl(cfl,u,v,w,dt)
      nekf_cfl = cfl

      return
      end
c-----------------------------------------------------------------------
      real function nekf_uf(u,v,w)

      real u(*), v(*), w(*)

      call nekuf(u,v,w)

      return
      end
c-----------------------------------------------------------------------

      integer function nekf_lglel(e)

      integer e

      include 'SIZE'
      include 'PARALLEL'

      nekf_lglel = lglel(e)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_uic(ifld)

      include 'SIZE'
      include 'TSTEP'

      ifield_ = ifield
      ifield = ifld
      call nekuic
      ifield = ifield_

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_ifoutfld(iswitch)

      include 'SIZE'
      include 'TSTEP'

      ifoutfld = .true.
      if (iswitch .eq. 0) ifoutfld = .false. 

      return
      end
c-----------------------------------------------------------------------

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
      include 'NEKINTF'

      if (id .eq. 'nelv') then 
         ptr = loc(nelv)
      elseif (id .eq. 'lelt') then 
         ptr = loc(llelt)
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
       elseif (id .eq. 'boundaryIDt') then
         ptr = loc(boundaryIDt(1,1)) 
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
      elseif (id .eq. 'qtl') then
         ptr = loc(qtl(1,1,1,1)) 
      elseif (id .eq. 'time') then
         ptr = loc(time)
      elseif (id .eq. 'ifgetu') then
         ptr = loc(getu)
      elseif (id .eq. 'ifgetp') then
         ptr = loc(getp)
      elseif (id .eq. 'ifgett') then
         ptr = loc(gett)
      elseif (id .eq. 'ifgetps') then
         ptr = loc(getps)
      else
         write(6,*) 'ERROR: nek_ptr cannot find ', id
         call exitt 
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_setup(comm_in,path_in, session_in, npscal_in) 

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'NEKINTF'

      integer comm_in, iftmsh_in
      character session_in*(*),path_in*(*)

      common /rdump/ ntdump

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

      llelt = lelt

      call setupcomm(comm_in,newcomm,newcommg,path_in,session_in)
      call iniproc()

      etimes = dnekclock()
      istep  = 0

      call initdim ! Initialize / set default values.
      call initdat
      call files

      call read_re2_hdr(ifbswap)

      call setDefaultParam
      param(1)  = 1.0
      param(2)  = 1.0
      param(7)  = 1.0
      param(8)  = 1.0
      param(27) = 1  ! torder 1 to save mem
      param(32) = 1  ! read only vel BC from re2
      param(99) = -1 ! no dealiasing to save mem

      ifflow = .true.
      iftran = .true.
      ifheat = .false.
      ifvo   = .true.
      ifpo   = .true.

      if (npscal_in .gt. 0) then
        ifheat = .true.
        if(nelgt.ne.nelgv) then
          param(32) = 2
          iftmsh(2) = .true.
        endif
        npscal = npscal_in - 1
        param(23) = npscal
        ifto   = .true.       
        do i = 1,npscal
          ifpsco(i) = .true.
        enddo 
      endif

      call bcastParam

      call usrdat0

      call chkParam
      call mapelpr 
      call read_re2_data(ifbswap)

      call izero(boundaryID, size(boundaryID))
      do iel = 1,nelv
      do ifc = 1,2*ndim
         boundaryID(ifc,iel) = bc(5,ifc,iel,1)
      enddo
      enddo
      call izero(boundaryIDt, size(boundaryIDt))
      if(nelgt.ne.nelgv) then 
        do iel = 1,nelt
        do ifc = 1,2*ndim
           boundaryIDt(ifc,iel) = bc(5,ifc,iel,2)
        enddo
        enddo
      endif

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

      call findSYMOrient

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

      call dofcnt

      p0thn = p0th
      ntdump=0

      call flush(6)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_resetio()

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      real ts
      integer npscals, p63s
      logical ifxyos, ifvos, ifpos, iftos, ifpscos(ldimt1)
      common /ros/  ts
      common /ios/  npscals, p63s
      common /ifos/ ifxyos, ifvos, ifpos, iftos, ifpscos

      time = ts

      param(63) = p63s

      npscal = npscals
      ifxyo  = ifxyos 
      ifvo   = ifvos 
      ifpo   = ifpos 
      ifto   = iftos 
      do i = 1,ldimt1
        ifpsco(i) = ifpscos(i) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_setio(ttime, xo, vo, po, so, ns, fp64)

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      real ttime
      integer xo, vo, po, so, fp64

      real ts
      integer npscals, p63s
      logical ifxyos, ifvos, ifpos, iftos, ifpscos(ldimt1)

      common /ros/  ts
      common /ios/  npscals, p63s
      common /ifos/ ifxyos, ifvos, ifpos, iftos, ifpscos

      if(ns.gt.ldimt) call exitti('nekf_setifo: ns > ldimt$',ns) 

      ts = time
      time = ttime

      p63s = param(63)
      param(63) = fp64 
 
      npscals = npscal
      ifxyos  = ifxyo
      ifvos   = ifvo
      ifpos   = ifpo
      iftos   = ifto
      do i = 1,ldimt1
        ifpscos(i)  = ifpsco(i) 
      enddo

      npscal = ns-1 
      ifxyo  = .false.
      ifvo   = .false.
      ifpo   = .false.
      ifto   = .false.
      do i = 1,ldimt1
        ifpsco(i) = .false.
      enddo
 
      if(xo.ne.0) ifxyo = .true.
      if(vo.ne.0) ifvo  = .true.
      if(po.ne.0) ifpo  = .true.
      if(so.ne.0) then
        ifto = .true.
        do i = 1,npscal
          ifpsco(i) = .true.
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_outfld(suffix)

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      character*3 suffix
      common /scrcg/ pm1(lx1,ly1,lz1,lelv)

      call copy(pm1,pr,nx1*ny1*nz1*nelv)
      call outfld(suffix)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_restart(rfile,l)

      character*(l) rfile

      include 'SIZE'
      include 'RESTART'
      include 'INPUT'

      call blank(initc(1),132)
      call chcopy(initc(1),rfile,l)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_end()

      call nek_end()

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
      subroutine nekf_setics()

      include 'SIZE'
      include 'RESTART'
      include 'NEKINTF'

      call setics()
      getu = 1
      getp = 1
      gett = 1
 
      if (.not. ifgetu) getu = 0 
      if (.not. ifgetp) getp = 0
      if (.not. ifgett) gett = 0

      return
      end
c-----------------------------------------------------------------------
      integer function nekf_bcmap(bID, ifld)

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      integer bID, ifld
      character*3 c

      if (bID < 1) then ! not a boundary
        nekf_bcmap = 0
        return 
      endif 

      ibc = 0 
      c = cbc_bmap(bID, ifld)

      if (ifld.eq.1) then
        if (c.eq.'W  ') then 
          ibc = 1
        else if (c.eq.'v  ') then 
          ibc = 2
        else if (c.eq.'o  ' .or. c.eq.'O  ') then 
          ibc = 3
        else if (c.eq.'SYX') then 
          ibc = 4
        else if (c.eq.'SYY') then 
          ibc = 5
         else if (c.eq.'SYZ') then 
          ibc = 6
        endif
      else if(ifld.gt.1) then
        if (c.eq.'t  ') then 
          ibc = 1
        else if (c.eq.'f  ') then 
          ibc = 2
        else if (c.eq.'o  ' .or. c.eq.'O  ' .or. c.eq.'I  ') then 
          ibc = 3
        endif
      endif

c      write(6,*) ifld, 'bcmap: ', bID, c, ibc

      if (ibc.eq.0) then
        write(6,*) 'Found unsupport BC type:', c
        call exitt 
      endif

      nekf_bcmap = ibc

      return
      end
c-----------------------------------------------------------------------
      subroutine findSYMOrient

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      integer bID
      logical ifalgn,ifnorx,ifnory,ifnorz

      do iel=1,nelt
      do ifc=1,2*ndim
         if (cbc(ifc,iel,1).eq.'SYM') then
           bID = boundaryID(ifc,iel) 
           call chknord(ifalgn,ifnorx,ifnory,ifnorz,ifc,iel) 
           if (ifnorx) cbc_bmap(bID, 1) = 'SYX'
           if (ifnory) cbc_bmap(bID, 1) = 'SYY'
           if (ifnorz) cbc_bmap(bID, 1) = 'SYZ'
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function nekf_nbid(isTmsh)

      include 'SIZE'
      include 'TOTAL'

      n = 0
      if(isTmsh.eq.1) then
        do iel = 1,nelt
        do ifc = 1,2*ndim
           n = max(n,boundaryIDt(ifc,iel))
        enddo
        enddo
      else
        do iel = 1,nelv
        do ifc = 1,2*ndim
           n = max(n,boundaryID(ifc,iel))
        enddo
        enddo
      endif

      nekf_nbid = iglmax(n,1)

      return
      end
c-----------------------------------------------------------------------
      integer*8 function nekf_set_vert(npts, isTmsh)

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      integer npts, isTmsh

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer vertex

      integer*8 ngv

      nx  = npts**(1./ndim)
      nel = nelt
      if (isTmsh.eq.0) nel = nelv
      call set_vert(glo_num,ngv,nx,nel,vertex,.false.)

      nekf_set_vert = ngv

      return
      end
c-----------------------------------------------------------------------

#include "bcType.h"

c-----------------------------------------------------------------------
c
c NEK5000 Interface
c
c-----------------------------------------------------------------------
      subroutine nekf_ptr(ptr,id,len)

      implicit none

      integer len 
      character*(len) id
      integer*8 i8
      pointer(ptr,i8)

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'
      include 'HSMG'

      if (id .eq. 'nelv') then 
         ptr = loc(nelv)
      elseif (id .eq. 'lelt') then 
         llelt = lelt
         ptr = loc(llelt)
      elseif (id .eq. 'ldimt') then 
         lldimt = ldimt            
         ptr = loc(lldimt)
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
      elseif (id .eq. 'vmult') then
         ptr = loc(vmult)
      elseif (id .eq. 'tmult') then
         ptr = loc(tmult)
      elseif (id .eq. 'cb_scnrs') then
         ptr = loc(sc_nrs(1))
      elseif (id .eq. 'p0th') then
         ptr = loc(p0th)
      elseif (id .eq. 'dp0thdt') then
         ptr = loc(dp0thdt)
      elseif (id .eq. 'wx') then
         ptr = loc(wx(1,1,1,1))
      elseif (id .eq. 'wy') then
         ptr = loc(wy(1,1,1,1))
      elseif (id .eq. 'wz') then
         ptr = loc(wz(1,1,1,1))
      elseif (id .eq. 'bfx') then
         ptr = loc(bfx(1,1,1,1))
      elseif (id .eq. 'bfy') then
         ptr = loc(bfy(1,1,1,1))
      elseif (id .eq. 'bfz') then
         ptr = loc(bfz(1,1,1,1))
      elseif (id .eq. 'bq') then
         ptr = loc(bq(1,1,1,1,1))
      else
         write(6,*) 'ERROR: nek_ptr cannot find ', id
         call exitt 
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_bootstrap(comm_in,path_in,session_in,mesh_in)

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'NEKINTF'

      integer comm_in
      character session_in*(*),path_in*(*)
      character mesh_in*(*)

      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest 

      character*1  re2fle1(132)
      equivalence  (RE2FLE,re2fle1)

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

      call usrdat0

      istep  = 0
      call initdim ! Initialize / set default values.
      call initdat
      call files

      lp = 0 !ltrunc(PATH,132)
      call chcopy(re2fle1(lp+1),mesh_in,len(mesh_in))
      ls = lp + len(mesh_in)
      call blank(re2fle1(ls+1),len(re2fle)-ls)


      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_setup(ifflow_in,
     $                      npscal_in, idpss_in, p32, mpart, contol,
     $                      rho, mue, rhoCp, lambda, stsform) 

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'NEKINTF'

      integer iftmsh_in, ifflow_in, mpart, p32
      integer idpss_in(*)
      real rho, mue, rhoCp, lambda, contol
      integer stsform

      common /rdump/ ntdump

      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer*8 vertex

      etimes = dnekclock_sync()

      call read_re2_hdr(ifbswap, .true.)

      if(ndim.eq.2) call exitti('Mesh has to be 3D!$', ndim) 

      call setDefaultParam
      loglevel   = 1
      cpfld(1,2) = rho
      cpfld(1,1) = mue
      cpfld(2,2) = rhoCp
      cpfld(2,1) = lambda

      param(27) = 1  ! torder 1 to save mem
      param(32) = p32 ! number of BC fields read from re2
      param(99) = -1 ! no dealiasing to save mem

      meshPartitioner = mpart 
      connectivityTol = contol

      ifflow = .true.
      if(ifflow_in.eq.0) ifflow = .false.
      iftran = .true.
      ifheat = .false.
      ifvo   = .true.
      ifpo   = .true.
      if(stsform.eq.1) ifstrs = .true.

      if (npscal_in .gt. 0) then
        ifheat = .true.
        if(nelgt.ne.nelgv) iftmsh(2) = .true.
        if(nelgt.ne.nelgv .and. param(32).eq.1) param(32) = 2 
        npscal = npscal_in - 1
        param(23) = npscal
        ifto = .true.    
        call icopy(idpss, idpss_in, npscal+1)
        do i = 1,npscal
          ifpsco(i) = .true.
        enddo 
      endif

      call bcastParam

      call chkParam
      call mapelpr 
      call read_re2_data(ifbswap, .true., .true., .true.)

      ifld_bId = 2
      if(ifflow) ifld_bId = 1
      do iel = 1,nelv
      do ifc = 1,2*ndim
         boundaryID(ifc,iel) = -1
         if(bc(5,ifc,iel,ifld_bId).gt.0)
     $     boundaryID(ifc,iel) = bc(5,ifc,iel,ifld_bId)
      enddo
      enddo
      if(nelgt.ne.nelgv) then 
        do iel = 1,nelt
        do ifc = 1,2*ndim
         boundaryIDt(ifc,iel) = -1
         if(bc(5,ifc,iel,2).gt.0)
     $     boundaryIDt(ifc,iel) = bc(5,ifc,iel,2)
        enddo
        enddo
      endif


      call setvar          ! Initialize most variables

      igeom = 2
      call setup_topo      ! Setup domain topology

      if(.not. ifflow) then
        call rone(vmult,lx1*ly1*lz1*nelv)
        ifield = 1
        call dssum(vmult,lx1,ly1,lz1)
        call invcol1(vmult,lx1*ly1*lz1*nelv)
      endif

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

      call setlog(.false.)  ! Initalize logical flags

      call bcmask  ! Set BC masks for Dirichlet boundaries.

      ifield = 1

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

      call dofcnt

      p0thn = p0th
      ntdump=0

      etimeSetup = dnekclock_sync() - etimes
      if(nio.eq.0) write(6,999) etimeSetup 
 999  format(' nek setup done in ', 1p1e13.4, ' s')
      if(nio.eq.0) write(6,*) 
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

      character suffix*(*)

      common /scrcg/ pm1(lx1,ly1,lz1,lelv)

      call copy(pm1,pr,nx1*ny1*nz1*nelv)
      call outfld(suffix)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_storesol()

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      parameter(ltot=lx1*ly1*lz1*lelt)
      common /outtmp/  w1(ltot),w2(ltot),w3(ltot),wp(ltot)
     &                ,wt(ltot,ldimt)

      ntot1  = lx1*ly1*lz1*nelt

      call copy(w1,vx,ntot1)
      call copy(w2,vy,ntot1)
      call copy(w3,vz,ntot1)
      call copy(wp,pr,ntot1)
      do i = 1,ldimt
         call copy(wt(1,i),t(1,1,1,1,i),ntot1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_restoresol()

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      parameter(ltot=lx1*ly1*lz1*lelt)
      common /outtmp/  w1(ltot),w2(ltot),w3(ltot),wp(ltot)
     &                ,wt(ltot,ldimt)

      ntot1  = lx1*ly1*lz1*nelt

      call copy(vx,w1,ntot1)
      call copy(vy,w2,ntot1)
      call copy(vz,w3,ntot1)
      call copy(pr,wp,ntot1)
      do i = 1,ldimt
         call copy(t(1,1,1,1,i),wt(1,i),ntot1)
      enddo

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

      include 'SIZE'
      include 'DPROCMAP'

#ifdef DPROCMAP
#ifdef MPI
      call MPI_Win_free(dProcmapH, ierr)
#endif
#endif 
      !call nek_end()

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
      integer function nekf_bcmap(bID, ifld, ismesh)

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      integer bID, ifld, ismesh
      character*3 c

      if (bID < 1) then ! not a boundary
        nekf_bcmap = 0
        return 
      endif 

      ibc = 0 
      c = cbc_bmap(bID, ifld)

      if (ifld.eq.1) then
        if (c.eq.'W  ') then 
          ibc = p_bcTypeW 
        else if (c.eq.'v  ') then 
          ibc = p_bcTypeV 
          if(ismesh.eq.1) then
            ibc = p_bcTypeW 
          endif
        else if (c.eq.'int') then 
          ibc = p_bcTypeINT
          if(ismesh.eq.1) then
            ibc = p_bcTypeW 
          endif
        else if (c.eq.'o  ' .or. c.eq.'O  ') then 
          ibc = p_bcTypeO 
          if(ismesh.eq.1) then
            ! outflow remaps to SYM bounds for mesh solver
            ibc = p_bcTypeSYM 
          endif
        else if (c.eq.'on ' .or. c.eq.'ON ') then 
          ibc = p_bcTypeON 
          if(ismesh.eq.1) then
            ! outflow remaps to SYM bounds for mesh solver
            ibc = p_bcTypeSYM 
          endif
        else if (c.eq.'onx') then 
          ibc = p_bcTypeONX 
          if(ismesh.eq.1) then
            ! outflow remaps to SYM bounds for mesh solver
            ibc = p_bcTypeSYM 
          endif
        else if (c.eq.'ony') then 
          ibc = p_bcTypeONY 
          if(ismesh.eq.1) then
            ! outflow remaps to SYM bounds for mesh solver
            ibc = p_bcTypeSYM 
          endif
        else if (c.eq.'onz') then 
          ibc = p_bcTypeONZ 
          if(ismesh.eq.1) then
            ! outflow remaps to SYM bounds for mesh solver
            ibc = p_bcTypeSYM 
          endif
        else if (c.eq.'SYX') then 
          ibc = p_bcTypeSYMX 
        else if (c.eq.'SYY') then 
          ibc = p_bcTypeSYMY 
         else if (c.eq.'SYZ') then 
          ibc = p_bcTypeSYMZ 
         else if (c.eq.'SYM') then 
          ibc = p_bcTypeSYM
         else if (c.eq.'shx') then 
          ibc = p_bcTypeSHLX 
         else if (c.eq.'shy') then 
          ibc = p_bcTypeSHLY 
         else if (c.eq.'shz') then 
          ibc = p_bcTypeSHLZ 
         else if (c.eq.'shl') then 
          ibc = p_bcTypeSHL 
          if(ismesh.eq.1) then
            ! outflow remaps to SYM bounds for mesh solver
            ibc = p_bcTypeSYM 
          endif
         else if (c.eq.'mv ') then 
          ibc = p_bcTypeV 
        endif
      else if(ifld.gt.1) then
        if (c.eq.'t  ') then 
          ibc = p_bcTypeS 
        else if (c.eq.'int') then 
          ibc = p_bcTypeINTS 
        else if (c.eq.'o  ' .or. c.eq.'O  ' .or. c.eq.'I  ') then 
          ibc = p_bcTypeF0 
        else if (c.eq.'f  ') then 
          ibc = p_bcTypeF 
        endif
      endif

c      write(6,*) ifld, 'bcmap: ', bID, 'cbc: ', c, 'ibc: ', ibc

      if (ibc.eq.0) then
        write(6,*) 'Found unsupport BC type:', c
        call exitt 
      endif

      nekf_bcmap = ibc

      return
      end
c-----------------------------------------------------------------------
      integer function nekf_nbid(isTmsh)

      include 'SIZE'
      include 'TOTAL'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer sum
      integer*8 bid8(2*ldim*lelt)
      integer maxbid

#if 0
      if(isTmsh.eq.1) then
        n = 2*ndim*nelt
        do i = 1,n
           bid8(i) = boundaryIDt(i,1)
        enddo
      else
        n = 2*ndim*nelv
        do i = 1,n
           bid8(i) = boundaryID(i,1)
        enddo
      endif

      call fgslib_gs_unique(bid8, n, nekcomm, np)

      sum = 0 
      do i = 1,n
        if(bid8(i).gt.0) sum = sum + 1
      enddo

      nekf_nbid = iglsum(sum,1)
#else
      maxbid = 0
      if(isTmsh.eq.1) then
        n = 2*ndim*nelt
        do i = 1,n
           if(boundaryIDt(i,1) .gt. maxbid) maxbid = boundaryIDt(i,1) 
        enddo
      else
        n = 2*ndim*nelv
        do i = 1,n
           if(boundaryID(i,1) .gt. maxbid) maxbid = boundaryID(i,1) 
        enddo
      endif

      nekf_nbid = iglmax(maxbid,1)
#endif

      return
      end
c-----------------------------------------------------------------------
      integer*8 function nekf_set_vert(nx, isTmsh)

      include 'SIZE'
      include 'TOTAL'
      include 'NEKINTF'

      integer npts, isTmsh

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer*8 vertex

      integer*8 ngv

      nel = nelt
      if (isTmsh.eq.0) nel = nelv
      call set_vert(glo_num,ngv,nx,nel,vertex,.false.)

      nekf_set_vert = ngv

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_dssum(u)
      include 'SIZE'
      include 'TOTAL'

      ifld = ifield
      ifield = 1
      call dssum(u,lx1,ly1,lz1)
      ifield = ifld 

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_gen_bcmap()
c
c     generate cbc_bmap mapping a boundaryID to boundary condition
c     note: 
c       * boundaryID is index-1 and contiguous
c       * each boundary condition defines a boundaryID 
c
      include 'SIZE'
      include 'TOTAL'

      integer bID, bcID
      integer map(p_velNBcType)
      integer ibc_bmap(lbid, ldimt1) 

      logical ifalg,ifnorx,ifnory,ifnorz
      character*3 cb 

      call ifill(boundaryID, -1, size(boundaryID))
      call ifill(boundaryIDt,-1, size(boundaryIDt))

      call izero(map, size(map))

      if(.not.ifflow .and. .not.ifheat) return 

      if(ifflow) then
        do iel = 1,nelv
        do ifc = 1,2*ndim
           cb = cbc(ifc,iel,1) 
           call chknord(ifalg,ifnorx,ifnory,ifnorz,ifc,iel)

           if(cb.eq.'W  ') map(p_bcTypeW) = 1
           if(cb.eq.'int') map(p_bcTypeINT) = 1
           if(cb.eq.'v  ') map(p_bcTypeV) = 1
           if(cb.eq.'mv ') map(p_bcTypeMV) = 1

           if(cb.eq.'o  ' .or. cb.eq.'O  ') map(p_bcTypeO) = 1
           if(cb.eq.'on ' .or. cb.eq.'ON ') then
             if (ifnorx) map(p_bcTypeONX) = 1 
             if (ifnory) map(p_bcTypeONY) = 1 
             if (ifnorz) map(p_bcTypeONZ) = 1 
             if (.not.ifalg) map(p_bcTypeON) = 1
           endif

           if(cb.eq.'SYM') then
             if (ifnorx) map(p_bcTypeSYMX) = 1 
             if (ifnory) map(p_bcTypeSYMY) = 1 
             if (ifnorz) map(p_bcTypeSYMZ) = 1 
             if (.not.ifalg) map(p_bcTypeSYM) = 1
           endif

           if(cb.eq.'shl') then
             if (ifnorx) map(p_bcTypeSHLX) = 1 
             if (ifnory) map(p_bcTypeSHLY) = 1 
             if (ifnorz) map(p_bcTypeSHLZ) = 1 
             if (.not.ifalg) map(p_bcTypeSHL) = 1
           endif
        enddo
        enddo
      elseif(ifheat) then
        do iel = 1,nelv
        do ifc = 1,2*ndim
           cb = cbc(ifc,iel,1)
           if(cb.eq.'int') map(p_bcTypeINTS) = 1
           if(cb.eq.'t  ') map(p_bcTypeS) = 1
           if(cb.eq.'I  ' .or. cb.eq.'O  ') map(p_bcTypeF0) = 1
           if(cb.eq.'f  ') map(p_bcTypeF) = 1
        enddo
        enddo
      endif

      ! assign each bcType a bID
      bID = 1
      do i = 1,p_velNBcType
        map(i) = iglmax(map(i),1)
        if(map(i).gt.0) then
          map(i) = bID
          bID = bID + 1
        endif 
      enddo
 
      ierr = 0
      ifld = 1
      do iel = 1,nelv
      do ifc = 1,2*ndim
         cb = cbc(ifc,iel,ifld)
         call chknord(ifalg,ifnorx,ifnory,ifnorz,ifc,iel)
 
         if(cb.eq.'W  ') then
           boundaryID(ifc,iel) = map(p_bcTypeW) 
         else if(cb.eq.'int') then
           boundaryID(ifc,iel) = map(p_bcTypeINT) 
         else if(cb.eq.'v  ') then
           boundaryID(ifc,iel) = map(p_bcTypeV) 
         else if(cb.eq.'mv ') then
           boundaryID(ifc,iel) = map(p_bcTypeMV) 
         else if(cb.eq.'o  ' .or. cb.eq.'O  ') then
           boundaryID(ifc,iel) = map(p_bcTypeO) 
         else if(cb.eq.'on ' .or. cb.eq.'ON ') then
           if (ifnorx) boundaryID(ifc,iel) = map(p_bcTypeONX) 
           if (ifnory) boundaryID(ifc,iel) = map(p_bcTypeONY) 
           if (ifnorz) boundaryID(ifc,iel) = map(p_bcTypeONZ) 
           if (.not.ifalg) boundaryID(ifc,iel) = map(p_bcTypeON) 
         else if(cb.eq.'SYM') then
           if (ifnorx) boundaryID(ifc,iel) = map(p_bcTypeSYMX) 
           if (ifnory) boundaryID(ifc,iel) = map(p_bcTypeSYMY) 
           if (ifnorz) boundaryID(ifc,iel) = map(p_bcTypeSYMZ) 
           if (.not.ifalg) boundaryID(ifc,iel) = map(p_bcTypeSYM) 
         else if(cb.eq.'shl') then
           if (ifnorx) boundaryID(ifc,iel) = map(p_bcTypeSHLX) 
           if (ifnory) boundaryID(ifc,iel) = map(p_bcTypeSHLY) 
           if (ifnorz) boundaryID(ifc,iel) = map(p_bcTypeSHLZ) 
           if (.not.ifalg) boundaryID(ifc,iel) = map(p_bcTypeSHL) 
         else
           if(cb.ne.'E  ' .and. cb.ne.'P  ') then
             ierr = 1
             write(6,*) 'b/c(1) of type ', cb
             goto 99
           endif
         endif
      enddo
      enddo
 99   call err_chk(ierr, 'Invalid velocity boundary condition type!$')

      if(map(p_bcTypeW).gt.0)
     $  cbc_bmap(map(p_bcTypeW), ifld) = 'W  '
      if(map(p_bcTypeINT).gt.0)
     $  cbc_bmap(map(p_bcTypeINT), ifld) = 'int'
      if(map(p_bcTypeV).gt.0)
     $  cbc_bmap(map(p_bcTypeV), ifld) = 'v  '
      if(map(p_bcTypeMV).gt.0)
     $   cbc_bmap(map(p_bcTypeMV), ifld) = 'mv '
      if(map(p_bcTypeO).gt.0)
     $   cbc_bmap(map(p_bcTypeO), ifld) = 'o  '
      if(map(p_bcTypeON).gt.0)
     $   cbc_bmap(map(p_bcTypeON), ifld) = 'on '
      if(map(p_bcTypeONX).gt.0)
     $   cbc_bmap(map(p_bcTypeONX), ifld) = 'onx'
      if(map(p_bcTypeONY).gt.0)
     $   cbc_bmap(map(p_bcTypeONY), ifld) = 'ony'
      if(map(p_bcTypeONZ).gt.0)
     $   cbc_bmap(map(p_bcTypeONZ), ifld) = 'onz'
      if(map(p_bcTypeSYMX).gt.0)
     $   cbc_bmap(map(p_bcTypeSYMX), ifld) = 'SYX'
      if(map(p_bcTypeSYMY).gt.0)
     $   cbc_bmap(map(p_bcTypeSYMY), ifld) = 'SYY'
      if(map(p_bcTypeSYMZ).gt.0)
     $   cbc_bmap(map(p_bcTypeSYMZ), ifld) = 'SYZ'
      if(map(p_bcTypeSYM).gt.0)
     $   cbc_bmap(map(p_bcTypeSYM), ifld) = 'SYM'
      if(map(p_bcTypeSHLX).gt.0)
     $   cbc_bmap(map(p_bcTypeSHLX), ifld) = 'shx'
      if(map(p_bcTypeSHLY).gt.0)
     $   cbc_bmap(map(p_bcTypeSHLY), ifld) = 'shy'
      if(map(p_bcTypeSHLZ).gt.0)
     $   cbc_bmap(map(p_bcTypeSHLZ), ifld) = 'shz'
      if(map(p_bcTypeSHL).gt.0)
     $   cbc_bmap(map(p_bcTypeSHL), ifld) = 'shl'

c      write(6,*) 'vel cbc_bmap: ', (cbc_bmap(i,1), i=1,p_velNBcType)
 
      do ifld = 2,nfield
        ierr = 0
        if(idpss(ifld-1).lt.0 .or. iftmsh(ifld)) goto 199
        call izero(ibc_bmap, size(ibc_bmap))
        do iel  = 1,nelv
        do ifc  = 1,2*ndim
          bID = boundaryID(ifc,iel)
          if(bID.gt.0) then
            cb = cbc(ifc,iel,ifld) 
            if(cb.eq.'t  ') then
              bcID = p_bcTypeS 
            else if(cb.eq.'int') then
              bcID = p_bcTypeINTS
            else if(cb.eq.'I  ' .or. cb.eq.'O  ') then
              bcID = p_bcTypeF0 
            else if(cb.eq.'f  ') then
              bcID = p_bcTypeF 
            else
              if(cb.ne.'E  ' .and. cb.ne.'P  ') then
                ierr = 1
                write(6,*) 'b/c(2) of type ', cb
              endif
            endif 
            ibc_bmap(bID, ifld) = bcID 
          endif          
        enddo
        enddo
        call err_chk(ierr, 'Invalid scalar boundary condition type!$')

        do bID = 1,p_scalNBcType
           bcID = iglmax(ibc_bmap(bID, ifld),1)
           if(bcID.eq.p_bcTypeINTS) cbc_bmap(bID, ifld) = 'int' 
           if(bcID.eq.p_bcTypeS) cbc_bmap(bID, ifld) = 't  ' 
           if(bcID.eq.p_bcTypeF0) cbc_bmap(bID, ifld) = 'I  ' 
           if(bcID.eq.p_bcTypeF) cbc_bmap(bID, ifld) = 'f  ' 
        enddo

c        write(6,*) ifld, 't cbc_bmap: ', (cbc_bmap(i,ifld), i=1,6)
 199    continue
      enddo

      ! cht
      ifld = 2
      if(idpss(ifld-1).gt.-1 .and. iftmsh(ifld)) then
        ierr = 0
        call izero(map, size(map))
 
        do iel = 1,nelt
        do ifc = 1,2*ndim
           cb = cbc(ifc,iel,ifld) 
           if(cb.eq.'t  ') then
             map(p_bcTypeS) = 1
           else if(cb.eq.'int') then 
             map(p_bcTypeINTS) = 1
           else if(cb.eq.'I  ' .or. cb.eq.'O  ') then 
             map(p_bcTypeF0) = 1
           else if(cb.eq.'f  ') then
             map(p_bcTypeF) = 1
           else 
             if(cb.ne.'E  ' .and. cb.ne.'P  ') then
               ierr = 1
               write(6,*) 'b/c(3) of type ', cb
             endif
           endif
        enddo
        enddo
        call err_chk(ierr, 'Invalid temp boundary condition type!$')
 
        bid = 1
        do i = 1,p_scalNBcType
          map(i) = iglmax(map(i),1)
          if(map(i).gt.0) then
            map(i) = bid
            bid = bid + 1
          endif 
        enddo
 
        do iel = 1,nelt
        do ifc = 1,2*ndim
           cb = cbc(ifc,iel,ifld) 
           if(cb.eq.'int')  boundaryIDt(ifc,iel) = map(p_bcTypeINTS) 
           if(cb.eq.'t  ')  boundaryIDt(ifc,iel) = map(p_bcTypeS) 
           if(cb.eq.'I  ' .or. cb.eq.'O  ') 
     $       boundaryIDt(ifc,iel) = map(p_bcTypeF0) 
           if(cb.eq.'f  ') boundaryIDt(ifc,iel) = map(p_bcTypeF)  
        enddo
        enddo

        if(map(p_bcTypeINTS).gt.0) cbc_bmap(map(1), ifld) = 'int'
        if(map(p_bcTypeS).gt.0) cbc_bmap(map(1), ifld) = 't  '
        if(map(p_bcTypeF0).gt.0) cbc_bmap(map(2), ifld) = 'I  '
        if(map(3).gt.p_bcTypeF) cbc_bmap(map(3), ifld) = 'f  '

c        write(6,*) 'cht t cbc_bmap:', (cbc_bmap(i,ifld), i=1,6)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nekf_scptr(id,ptr)

      implicit none

      integer id
      integer*8 i8
      pointer(ptr,i8)

      include 'SIZE'
      include 'NEKINTF'
     
      ptr = nrs_scptr(id)  

      return
      end
C----------------------------------------------------------------------
C
C     Generate geometric factors without updating coords
C
C----------------------------------------------------------------------
      subroutine nekf_updggeom()
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'GEOM'
      include 'WZ'

      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)

      ifld_save = ifield
      ifield = 1

      CALL LAGMASS
      CALL GEOM1 (XM3,YM3,ZM3)
      CALL GEOM2
      CALL UPDMSYS (1)
      CALL VOLUME
      CALL SETINVM
      CALL SETDEF
      CALL SFASTAX

      ifield = ifld_save

      return
      end

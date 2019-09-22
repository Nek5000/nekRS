      subroutine vertical_mean(ff,i_name)
      include 'SIZE'
      include 'TOTAL'

      parameter(nelz=16)

      real ff(lx1,ly1,lz1,lelt)
      real fbar(lz1,nelz)
      real zbar(lz1,nelz)
      real wght(lz1,nelz)
      real work(lz1,nelz)

      integer e,eg,ex,ey,ez,f

      nelxy = nelgv/nelz 

      call rzero(fbar,lz1*nelz)
      call rzero(zbar,lz1*nelz)
      call rzero(wght,lz1*nelz)

      f = 5 ! to evaluate surface Jacobian 
      do e=1,nelv
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelxy,1,nelz)
         do k=1,nz1
         do i=1,nx1*ny1
             fbar(k,ez) = fbar(k,ez)+area(i,1,f,e)*ff(i,1,k,e)
             zbar(k,ez) = zbar(k,ez)+area(i,1,f,e)*zm1(i,1,k,e)
             wght(k,ez) = wght(k,ez)+area(i,1,f,e)
         enddo
         enddo
      enddo

      call gop(fbar,work,'+  ',lz1*nelz)
      call gop(zbar,work,'+  ',lz1*nelz)
      call gop(wght,work,'+  ',lz1*nelz)

      do i=1,lz1*nelz
        fbar(i,1)=fbar(i,1)/wght(i,1) 
        zbar(i,1)=zbar(i,1)/wght(i,1)
      enddo

      if(nid.eq.0)then
        if(i_name.eq.1)OPEN(10,file="ver_uzte.dat",position="append")
        if(i_name.eq.2)OPEN(10,file="ver_epst.dat",position="append")
        if(i_name.eq.3)OPEN(10,file="ver_epsv.dat",position="append")
        if(i_name.eq.4)OPEN(10,file="ver_temp.dat",position="append")
        if(i_name.eq.5)OPEN(10,file="ver_dtdz.dat",position="append")

        do k=1,nelz
        do i=1,lz1-1
          write(10,*)zbar(i,k),fbar(i,k)
        enddo   
        enddo
        CLOSE(10)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pdf_calc(ff,step,i_offset,i_name)
      include 'SIZE'
      include 'TOTAL'

      parameter(npdf=2401)

      real ff(lx1,ly1,lz1,lelt)
      real pdf(npdf)
      real work(npdf)
      real val, vol, offset, wght, norm
      integer e 

      call rzero(pdf,npdf)
      call rzero(work,npdf)

      if(i_offset==0) offset=0.0 
      if(i_offset==1) offset=int(npdf/2)*step 

!-----Volume for aspect ratio=1
      vol=atan(1.0)

      do e=1,nelv
        do k=1,nz1
          do j=1,ny1
            do i=1,nx1

               wght= bm1(i,j,k,e) 
               val= ff(i,j,k,e)

               do ipdf=1,npdf
                dm1=(ipdf -1 )*step-offset
                dm2= ipdf     *step-offset

                if((val.ge.dm1).and.(val.lt.dm2))then
                  pdf(ipdf)=pdf(ipdf)+wght
                endif 

               enddo    
            enddo
          enddo
        enddo
      enddo

      do i=1,npdf 
        pdf(i)=pdf(i)/vol
      enddo

      call gop(pdf,work,'+  ',npdf)

!------------------------------------------------------------ 
      if(nid.eq.0)then

        if(i_name.eq.1)OPEN(10,file="pdf_dudx.dat",position="append")
        if(i_name.eq.2)OPEN(10,file="pdf_epst.dat",position="append")
        if(i_name.eq.3)OPEN(10,file="pdf_epsv.dat",position="append")
        if(i_name.eq.4)OPEN(10,file="pdf_temp.dat",position="append")
        if(i_name.eq.5)OPEN(10,file="pdf_dtdx.dat",position="append")

        do ipdf=1,npdf
          write(10,*) (ipdf-1)*step-offset, pdf(ipdf)
        enddo
        CLOSE(10)

        norm=0.0
        do ipdf=1,npdf
          norm=norm+pdf(ipdf)
        enddo
        if(i_name.eq.1)write(77,*)'dudx',norm
        if(i_name.eq.2)write(77,*)'epsT',norm
        if(i_name.eq.3)write(77,*)'epsv',norm
        if(i_name.eq.4)write(77,*)'Temp',norm
        if(i_name.eq.5)write(77,*)'dTdx',norm

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pdf_calc_bulk(ff,step,i_offset,i_name)
      include 'SIZE'
      include 'TOTAL'

      parameter(npdf=2401)

      real ff(lx1,ly1,lz1,lelt)
      real pdf(npdf)
      real work(npdf)
      real val, vol, offset, wght, norm, height1, height2, radius, radiusc
      integer e 

!-----Set arrays to zero
      call rzero(pdf,npdf)
      call rzero(work,npdf)

!-----Offset
      if(i_offset==0) offset=0.0 
      if(i_offset==1) offset=int(npdf/2)*step 

!-----Volume for aspect ratio=1
      radiusc=0.3
      height1=0.2
      height2=0.8       
      vol=4*atan(1.0)*radiusc**2*(height2-height1)

      do e=1,nelv
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
               wght= bm1(i,j,k,e) 
               val= ff(i,j,k,e)

	       radius=sqrt(xm1(i,j,k,e)**2+ym1(i,j,k,e)**2)
               IF (radius.le.radiusc) THEN 

!-----------------------------------------------------------------------
      if ((zm1(i,j,k,e).ge.height1).and.(zm1(i,j,k,e).le.height2)) then
                
               do ipdf=1,npdf
                dm1=(ipdf -1 )*step-offset
                dm2= ipdf     *step-offset

                if((val.ge.dm1).and.(val.lt.dm2))then
                  pdf(ipdf)=pdf(ipdf)+wght
                endif 

               enddo

               endif    
!-----------------------------------------------------------------------

               ENDIF
      enddo
      enddo
      enddo
      enddo

      do i=1,npdf 
        pdf(i)=pdf(i)/vol
      enddo

      call gop(pdf,work,'+  ',npdf)

      if(nid.eq.0) then
        if(i_name.eq.1)OPEN(10,file="pdf_dudx_b.dat",position="append")
        if(i_name.eq.2)OPEN(10,file="pdf_epst_b.dat",position="append")
        if(i_name.eq.3)OPEN(10,file="pdf_epsv_b.dat",position="append")
        if(i_name.eq.4)OPEN(10,file="pdf_temp_b.dat",position="append")
        if(i_name.eq.5)OPEN(10,file="pdf_dtdx_b.dat",position="append")

        do ipdf=1,npdf
          write(10,*) (ipdf-1)*step-offset, pdf(ipdf)
        enddo
        CLOSE(10)

        norm=0.0
        do ipdf=1,npdf
          norm=norm+pdf(ipdf)
        enddo
        if(i_name.eq.1)write(78,*)'dudx',norm
        if(i_name.eq.2)write(78,*)'epsT',norm
        if(i_name.eq.3)write(78,*)'epsv',norm
        if(i_name.eq.4)write(78,*)'Temp',norm
        if(i_name.eq.5)write(78,*)'dTdx',norm
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine volume_mean(ff,i_name)
      include 'SIZE'
      include 'TOTAL'

      real ff(lx1,ly1,lz1,lelt)
      real work1
      real vol, vol_check, mean_val1, mean_val2, mean_val3, mean_val4
      integer e 

      call rzero(mean_val1,1)
      call rzero(mean_val2,1)
      call rzero(mean_val3,1)
      call rzero(mean_val4,1)
      call rzero(vol_check,1)
      call rzero(work1,1)

!-----Volume for aspect ratio=1
      vol=atan(1.0)
      do e=1,nelv
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         vol_check=vol_check + bm1(i,j,k,e)
         mean_val1=mean_val1 + ff(i,j,k,e)*bm1(i,j,k,e)
         mean_val2=mean_val2 + (ff(i,j,k,e)**2)*bm1(i,j,k,e)
         mean_val3=mean_val3 + (ff(i,j,k,e)**3)*bm1(i,j,k,e)	       
         mean_val4=mean_val4 + (ff(i,j,k,e)**4)*bm1(i,j,k,e)
      enddo
      enddo
      enddo
      enddo

!-----Normalization      
      mean_val1=mean_val1/vol
      mean_val2=mean_val2/vol
      mean_val3=mean_val3/vol
      mean_val4=mean_val4/vol

!-----Gather over all processes (-> mpi_allreduce)
      call gop(mean_val1,work1,'+  ',1)
      call gop(mean_val2,work1,'+  ',1)
      call gop(mean_val3,work1,'+  ',1)
      call gop(mean_val4,work1,'+  ',1)
      call gop(vol_check,work1,'+  ',1)

!------------------------------------------------------------ 
      if(nid.eq.0)then

      if(i_name.eq.1)OPEN(10,file="mean_epst.dat",position="append")
      if(i_name.eq.2)OPEN(10,file="mean_epsv.dat",position="append")
      if(i_name.eq.3)OPEN(10,file="mean_uzte.dat",position="append")
      write(10,*)mean_val1,mean_val2,mean_val3,mean_val4,vol_check,vol
      CLOSE(10)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bulk_mean(ff,i_name)
      include 'SIZE'
      include 'TOTAL'

      real ff(lx1,ly1,lz1,lelt)
      real work1
      real vol, vol_check, mean_val1, mean_val2, mean_val3, mean_val4
      integer e 
      real height1, height2, radius, radiusc

      call rzero(mean_val1,1)
      call rzero(mean_val2,1)
      call rzero(mean_val3,1)
      call rzero(mean_val4,1)
      call rzero(vol_check,1)
      call rzero(work1,1)

      radiusc=0.3
      height1=0.2
      height2=0.8
      vol=4*atan(1.0)*radiusc**2*(height2-height1)

      do e=1,nelv
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         radius=sqrt(xm1(i,j,k,e)**2+ym1(i,j,k,e)**2)
         if (radius.le.radiusc) then
            if (zm1(i,j,k,e).ge.height1 .and.
     $          zm1(i,j,k,e).le.height2) then
                   vol_check=vol_check + bm1(i,j,k,e)
	           mean_val1=mean_val1 + ff(i,j,k,e)*bm1(i,j,k,e)
	           mean_val2=mean_val2 + (ff(i,j,k,e)**2)*bm1(i,j,k,e)
	           mean_val3=mean_val3 + (ff(i,j,k,e)**3)*bm1(i,j,k,e)	       
	           mean_val4=mean_val4 + (ff(i,j,k,e)**4)*bm1(i,j,k,e)
	    endif
         endif
      enddo
      enddo
      enddo
      enddo

      mean_val1=mean_val1/vol
      mean_val2=mean_val2/vol
      mean_val3=mean_val3/vol
      mean_val4=mean_val4/vol

      call gop(mean_val1,work1,'+  ',1)
      call gop(mean_val2,work1,'+  ',1)
      call gop(mean_val3,work1,'+  ',1)
      call gop(mean_val4,work1,'+  ',1)
      call gop(vol_check,work1,'+  ',1)

      if(nid.eq.0)then
        if(i_name.eq.1)OPEN(10,file="mean_epst_b.dat",position="append")
        if(i_name.eq.2)OPEN(10,file="mean_epsv_b.dat",position="append")
        if(i_name.eq.3)OPEN(10,file="mean_uzte_b.dat",position="append")
        write(10,*)mean_val1,mean_val2,mean_val3,mean_val4,vol_check,vol
        CLOSE(10)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nusselt
      include 'SIZE'  
      include 'TOTAL' 

      parameter (lt=lx1*ly1*lz1*lelt)
      real dtdx(lt),dtdy(lt),dtdz(lt),w(lt)

      common /chkcmnr/ atime,timel,flux(2)
      common /chkcmni/ iesurf(0:2*lelt,2),ifsurf(0:2*lelt,2)

      integer icalld
      save    icalld
      data    icalld  /0/

      if (icalld.eq.0) then
         icalld = icalld + 1
         atime  = 0.
         timel  = time
         call find_srf(iesurf(0,1),ifsurf(0,1),5)    ! bottom face
         call find_srf(iesurf(0,2),ifsurf(0,2),6)    ! top    face
      endif

!-----Compute temperature gradients for surface fluxes
      call gradm1(dtdx,dtdy,dtdz,t)
      do k = 1,2
         flux(k)  = 0.0
         do isurf = 1,iesurf(0,k)
            ie    = iesurf(isurf,k)
            iface = ifsurf(isurf,k) 
            if (CBC(IFACE,IE,1) .NE .'E  ') then
               call surface_flux(dq,dtdx,dtdy,dtdz,ie,iface,w)
               flux(k)  = flux(k)  + dq
            end if
         enddo
      enddo
      call gop(flux,w,'+  ',2)

!-----Output of Nusselt number at top/bottom plate
      dtime = time  - timel
      atime = atime + dtime
      pi4=atan(1.0)
      if (nid.eq.0) then
        write(6,1) istep,time,atime,flux(1)/pi4,-flux(2)/pi4
      endif
    1 format(i6,' Nusselt',1p4e14.6)

      timel = time

      return
      end
c-----------------------------------------------------------------------
      subroutine find_srf(iesurf,ifsurf,inface)
c
c     Find list of surfaces over which the flux should be computed
c     Number of such surfaces is returned in iesurf(0).
c
      include 'SIZE'
      include 'TOTAL'

      integer iesurf(0:lelt),ifsurf(0:lelt)

      nsurf = 0
      nfaces = 2*ndim
      do ie=1,nelv
         nsurf         = nsurf+1
         iesurf(nsurf) = ie
         ifsurf(nsurf) = inface
      enddo
      iesurf(0) = nsurf
      ifsurf(0) = nsurf
      return
      end
c-----------------------------------------------------------------------


      PROGRAM rpaPlasmons
      implicit none

      REAL timeStart, timeFinish,timeStart_i,timeFinish_i

      COMPLEX*16 ima,g,chi_c
      REAL*8 Pi,kx,mu
      REAL*8 t,tz,tt,tpp,dop,temp,wmax,gamma,w,qx,qy,qz
      REAL*8 Vc,Aq,Vq,alpha
      REAL*8 muCalculo


      INTEGER ikx,iw,Nk,i,Nw,Nkz
      INTEGER narg,ppos,jj

      PARAMETER (Nk=400)
      PARAMETER (Nkz=1)
      PARAMETER (Nw=1500)

      CHARACTER(LEN=50) :: var, comm !var contiene las variables pasados por el bash
      CHARACTER(LEN=100) :: name_file='',comentarios='',name_file2=''
      CHARACTER(LEN=4) :: char_aux='', char_aux2=''
      CHARACTER(LEN=8) :: date
      CHARACTER(LEN=10) :: time

      REAL*8 cosk(-Nk:Nk), sink(-Nk:Nk)
      COMPLEX*16 p(0:Nw)

      call cpu_time(timeStart)
      call cpu_time(timeStart_i)

      ima=dcmplx(0.0d0,1.0d0)
      Pi  = 3.14159265359d0

c--------CREATING THE FILE NAME OVER THE BASH ARGUMENTS - OPTIONAL: YOU CAN ADD COMMENTS AT THE END-------------------
      narg=command_argument_count()
      IF(narg.GT.0) THEN
        DO i=1,narg
            CALL get_command_argument(i,var)
            SELECT CASE(i)
            CASE(1)
                READ(var,'(I5)') jj
                IF(jj.LT.10) THEN
                  name_file=TRIM(name_file)//'00000'//var
                ELSE IF (jj.LT.100) THEN
                    name_file=TRIM(name_file)//'0000'//var
                ELSE if(jj.LT.1000) then
                    name_file=TRIM(name_file)//'000'//var
				else if(jj.LT.10000) then
					name_file=TRIM(name_file)//'00'//var
				else
					name_file=TRIM(name_file)//'0'//var
                ENDIF
			CASE(2)
               READ(var,'(D11.6)') qx
            CASE(3)
                READ(var,'(D11.6)') qy
			CASE(4)
                READ(var,'(D11.6)') qz
            CASE(5)
               comm=var
            END SELECT
        ENDDO
      ENDIF

      IF(comm.EQ.'com') THEN
        WRITE(*,*)'Comentarios en el archivo='
        READ(*,'(a)') comentarios
      ENDIF

c--------CLEANING THE FILE NAME-------------------

      DO WHILE(.TRUE.)
        ppos=SCAN(name_file,'.',BACK=.true.)
        IF (ppos.GT.0) THEN
         name_file=TRIM(name_file(1:ppos-1))//TRIM(name_file(ppos+1:))
        ELSE
           EXIT
        ENDIF
      ENDDO

      WRITE(*,*)'Files: ', TRIM(name_file)

      OPEN(2,FILE=TRIM(name_file),STATUS='unknown')


c--------SETTING THE PARAMETERS-------------------

      t=0.4d0 
      tt=0.09d0*t 
      tpp=0.07d0*t

      tz=0.01d0*t

      dop=0.16d0 

	  temp=1.0d-4
      gamma=0.04d0

	  Vc=6.91d0

	  alpha=3.1

	  wmax=10.0d0*t


c--------PRECALCULATING COSINE AND SIN---------------------

      DO ikx=-Nk,Nk-1
        kx=DBLE(ikx)*Pi/DBLE(Nk)
        cosk(ikx)=dcos(kx)
        sink(ikx)=dsin(kx)
      ENDDO

c--------SHOWING PARAMETERS---------------------

	  write(*,*)dop,t,tt,tpp,tz,temp,Vc,alpha,Vq,qx,qy,qz

c--------CALCULATING MU---------------------

	  mu=muCalculo(dop,t,tt,tz,tpp,temp,Nk,Nkz)
	  write(*,*) 'MU=',mu
	  
c--------CALCULATING POLARIZATION---------------------

	  Aq=alpha*(2.0d0-dcos(qx)-dcos(qy))+1.0d0
      Vq=Vc/(Aq-dcos(qz))

	  CALL calcP (t,tt,tz,tpp,temp,gamma,
     &            qx,qy,qz,wmax,mu,p,cosk,sink,Nk,Nkz,Nw)

c--------CALCULATING SUSCEPTIBILITY---------------------

      do iw=0,Nw
          w=dble(iw)*wmax/dble(Nw)
          chi_c=2.0d0*p(iw)/(1.0d0-2.0d0*Vq*p(iw))

	  WRITE(2,22) jj,qx,qy,qz,w,dimag(p(iw)),dreal(p(iw)),
     &            dimag(chi_c),dreal(chi_c),Vq

      enddo

      call cpu_time(timeFinish_i)
      WRITE(*,'(a,F9.2,a)') 'Tiempo=', timeFinish_i-timeStart_i, 's'
      call cpu_time(timeStart_i)





 100  call cpu_time(timeFinish)
      WRITE(*,'(a,F9.2,a)') 'Tiempo total de ejecucion=',
     &  timeFinish-timeStart, 'seg'
      WRITE(*,'(a,I5,a,F6.2,a)') 'Tiempo total de ejecucion=',
     &          INT((timeFinish-timeStart)/60.0),'m',
     &(timeFinish-timeStart)-(INT((timeFinish-timeStart)/60.0)*60),'s'
	 
   10 format(6e14.5)
   22 format(I5,9e14.5)

      END



      SUBROUTINE calcP (t,tt,tz,tpp,temp,gamma,
     &                  qx,qy,qz,wmax,mu,p,cosk,sink,Nk,Nkz,Nw)

      REAL*8 t,mu,w,qx,qy,qz,temp,gamma,tt,tz,wmax,tpp
      REAL*8 ek,ekmq,nf,nfk,nfkpq,ekperp,ekpar,ekmqpar,ekmqperp,ekpp
      REAL*8 cosqy,cosq2y,sinqy,sinq2y,cosqx,sinqx,cosqz,sinqz,sz,cz,kz
      REAL*8 ekmqpp,kx,ky
      COMPLEX*16 ima, nume, denom,g

      INTEGER  Nw,Nk,Nkz
      INTEGER iw,ikx,iky,ikz

      REAL*8 cosk(-Nk:Nk-1), sink(-Nk:Nk-1)
      COMPLEX*16 p(0:Nw)

      normak=((2.0d0*DBLE(Nk))**2)*2.0d0*dble(Nkz)
		
      ima=dcmplx(0.0d0,1.0d0)
      Pi  = 3.14159265359d0


	  cosqx=dcos(qx)
	  cosqy=dcos(qy)
	  cosqz=dcos(qz)

	  sinqx=dsin(qx)
	  sinqy=dsin(qy)
	  sinqz=dsin(qz)

	  do iw=0,Nw
		p(iw)=dcmplx(0.0d0,0.0d0)
	  enddo

	  DO ikz=-Nkz, Nkz-1
		kz=dble(ikz)*Pi/dble(Nkz)
		cz = dcos(kz)
		sz = dsin(kz)
		DO ikx=-Nk, Nk-1
			kx=dble(ikx)*Pi/dble(Nk)
			DO iky=-Nk, Nk-1
				ky=dble(iky)*Pi/dble(Nk)
				ekpar=-2.0d0*t*(cosk(ikx)+cosk(iky))
     &                +4.0d0*tt*cosk(ikx)*cosk(iky)-mu
				ekperp=0.25d0*tz*((cosk(ikx)-cosk(iky))**2)!*cz
				ekpp=-2.0d0*tpp*(dcos(2.0d0*kx)+dcos(2.0d0*ky))!(c2x+c2y)

				ek=ekpar+ekperp+ekpp
				ek=ek*0.5d0 !effective mass

				ekmqpar=-2.0d0*t*((cosk(ikx)*cosqx+
     &                  sink(ikx)*sinqx)+(cosk(iky)*cosqy+
     &                  sink(iky)*sinqy))
     &                          +4.0d0*tt*(cosk(ikx)*cosqx+
     &                  sink(ikx)*sinqx)*(cosk(iky)*cosqy+
     &                  sink(iky)*sinqy)-mu
                        ekmqperp=0.25d0*tz*(((cosk(ikx)*cosqx+
     &                  sink(ikx)*sinqx)-(cosk(iky)*cosqy+
     &                  sink(iky)*sinqy))**2)!*&                  (cz*cosqz+sz*sinqz)
                ekmqpp=-2.0d0*tpp*(dcos(2.0d0*(kx-qx))+
     &                 dcos(2.0d0*(ky-qy)))
                        

				ekmq=ekmqpar+ekmqperp+ekmqpp
				ekmq=ekmq*0.5d0

				nfk=nf(ek,temp)
				nfkpq=nf(ekmq,temp)

				nume=nfkpq-nfk

				do iw=0,Nw
					w=dble(iw)*wmax/dble(Nw)
					denom=-w-ima*gamma+ekmq-ek

					g=nume/denom

					p(iw)=p(iw)+g
				enddo
            enddo
		ENDDO
	  ENDDO 
	  do iw=0,Nw
	    	p(iw)=p(iw)/normak
	  enddo

      RETURN
      END SUBROUTINE


      FUNCTION nf (ek,temp)

            REAL*8 temp,ek,nf,aa

            IF (temp.GT.0.0d0) THEN
                    aa=ek/temp
                    IF(aa.GT.60.0d0) THEN
                        nf=0.0d0
                    ELSE
                        IF(aa.LT.-60.0d0) THEN
                            nf=1.0d0
                        ELSE
                            nf = 1.0d0/(dexp(aa)+1.0d0)
                        ENDIF
                    ENDIF
                ELSE
                    aa=ek
                    IF(aa.GT.0.0d0) THEN
                        nf=0.0d0
                    ELSE
                        IF(aa.LT.0.0d0) THEN
                            nf=1.0d0
                        ELSE
                            nf = 0.5d0
                        ENDIF
                    ENDIF
                ENDIF


      END FUNCTION



      FUNCTION muCalculo (dop,t,tt,tz,tpp,temp,Nmu,Nkz)
        implicit none
        REAL*8 muA,muB,mu,muCalculo,Pi,normakDico,t,tpp
        REAL*8 dop,temp,kx,ky,cx,cy,eka,ekb,nf,delta,tt
        REAL*8 tz,cz,kz,ekpar,ekperp,nfk
        REAL*8 sumaA,sumaC,sumaB,c2x,c2y,ekpp

        INTEGER Nmu, Nd, ikx,iky,i,Nkz,ikz
        PARAMETER (Nd=30)

        REAL*8 cosk(-Nmu:Nmu),sink(-Nmu:Nmu)

        normakDico=(2.0d0*dble(Nkz))*(2.0d0*DBLE(Nmu))**2

        Pi  = 3.14159265359d0

        DO ikx=-Nmu,Nmu
          kx=DBLE(ikx)*Pi/DBLE(Nmu)
          cosk(ikx)=dcos(kx)
          sink(ikx)=dsin(kx)
        ENDDO

        muA=-10.0d0
        muB=10.0d0
        sumaA=0.0d0
        sumaB=0.0d0
        sumaC=0.0d0

        do ikz=-Nkz,Nkz-1
            kz=dble(ikz)*Pi/dble(Nkz)
            cz = dcos(kz)
            DO ikx=-Nmu, Nmu-1
                cx = cosk(ikx)
                kx=dble(ikx)*Pi/dble(Nmu)
				DO iky=-Nmu, Nmu-1
					ky=dble(iky)*Pi/dble(Nmu)
					cy = cosk(iky)

					ekpar=-2.0d0*t*(cx+cy)+4.0d0*tt*cx*cy
					ekperp=0.25d0*tz*((cx-cy)**2)!*cz
					ekpp=-2.0d0*tpp*(dcos(2.0d0*kx)+dcos(2.0d0*ky))

					eka=ekpar-muA+ekperp+ekpp
					
					nfk=nf(eka,temp)
					sumaA=nfk+sumaA

					ekb=ekpar-muB+ekperp+ekpp
					
					nfk=nf(ekb,temp)
					sumaB=nfk+sumaB
				enddo
            ENDDO
        ENDDO

        sumaA=1.0d0-dop-2.0d0*(sumaA/normakDico)
        sumaB=1.0d0-dop-2.0d0*(sumaB/normakDico)

        DO i=1,Nd
            mu=(muA+muB)*0.5d0
            sumaC=0.0d0

            do ikz=-Nkz,Nkz-1
                kz=dble(ikz)*Pi/dble(Nkz)
                cz = dcos(kz)
                DO ikx=-Nmu, Nmu-1
                    cx = cosk(ikx)
                    kx=dble(ikx)*Pi/dble(Nmu)
                    DO iky=-Nmu, Nmu-1
                        ky=dble(iky)*Pi/dble(Nmu)
                        cy = cosk(iky)

                        ekpar=-2.0d0*t*(cx+cy)+4.0d0*tt*cx*cy
                        ekperp=0.25d0*tz*((cx-cy)**2)!*cz
                        ekpp=-2.0d0*tpp*(dcos(2.0d0*kx)+dcos(2.0d0*ky))

                        eka=ekpar-mu+ekperp+ekpp
						
                        nfk=nf(eka,temp)

                        sumaC=nfk+sumaC
					enddo
				ENDDO
			ENDDO

            sumaC=1.0d0-dop-2.0d0*(sumaC/normakDico)

            IF(sumaA*sumaC.GT.0.0d0)THEN
				sumaA=sumaC
                muA=mu
            ELSE
                sumaB=sumaC
                muB=mu
            ENDIF
        ENDDO

        muCalculo=mu

      END FUNCTION







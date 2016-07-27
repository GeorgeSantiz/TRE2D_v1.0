MODULE resistfm
!  ******************************************************************************
!  SUBRUTINA PARA LEER PARAMETROS DEL MODELO  
!
!  SUBRUTINE PARAMETROS
!    SIN PARAMETROS DE ENTRADA/SALIDA
!
!  ******************************************************************************
!  SUBRUTINA PARA GENERAR UN ARCHIVO CON COORDENADAS QUE LIMITAN A LOS BLOQUES 
!  DE CONDUCTIVIDAD  
!
!  SUBROUTINE plot_model
!    SIN PARAMETROS DE ENTRADA/SALIDA
!
!  ******************************************************************************
!  SUBRUTINA PARA CREAR GENERAR LAS MEDICONES POSIBLES CON EL DISPOSITIVO WENNER
!  SUBROUTINE WN
!    SIN PARAMETROS DE ENTRADA/SALIDA
!
!  ******************************************************************************
!  SUBRUTINA PARA CREAR GENERAR LAS MEDICONES POSIBLES CON EL DISPOSITIVO WENNER-SHCLUMBERGER
!  SUBROUTINE WS
!    SIN PARAMETROS DE ENTRADA/SALIDA
!
!  ******************************************************************************
!  SUBRUTINA PARA CREAR GENERAR LAS MEDICONES POSIBLES CON EL DISPOSITIVO DIPOLO-DIPOLO
!  SUBROUTINE DD
!    SIN PARAMETROS DE ENTRADA/SALIDA
!
!  ************************************************************************************
!  SUBRUTINA PARA GENERAR LA MATRIZ DEL SISTEMA
!  SUBROUTINE MATRIX_A
!    SIN PARAMETROS DE ENTRADA/SALIDA
!
!  ************************************************************************************
!  SUBRUTINA PARA GENERAR LA MATRIZ PARA EL CALCULO DEL POTENCIAL PRIMARIO
!  SUBROUTINE MATRIX_M(COND)
!    PARAMETROS DE ENTRADA
!      COND: MATRIZ QUE CONTIENE A LA CONDUCTIVIDAD PRIMARIA
!
!  ************************************************************************************
!  SUBRUTINA PARA CONVERTIR LA MATRIZ "A" EN FORMATO BAND STORAGE
!  SUBRUTINE SPARSE
!    SIN PARAMETROS DE ENTRADA/SALIDA
!
!  ************************************************************************************
!  FUNCION DE ANDERSON(1979) PARA CALCULO DE LA TRANFORMADA DE HANKEL MEDIANTE FILTRADO
!  LINEAL
!  COMPLEX FUNCTION ZHANKS(ORD,MULT,FUN,TOL,NF,NEW)    
!    DETALLE DE LOS PARAMETROS EN EL ENCABEZADO DE LA SUBRUTINA  
!
!    PROGRAMER             DATE         
!   GEORGE SANTIZ       27/07/2016      CODIGO (EXCEPTO FUNCTION ZHANKS) 
!  MODULO CON SUBRUTINAS PARA EL MODELADO DE MEDIDAS DE RESISTIVIDAD APARENTE 
!  UTILIZANDO UN STENCIL OBTENIDO MEDIANTE UNA APROXIMACION POR DIFERENCIAS FINITAS
!  LAS SUBRUTINAS COMPARTEN VARIABLES GLOBALES PARA EVITAR EL PASO DE PARAMETROS 
!
!********************************************************************
!                    VARIABLES GLOBALES
!*********************************************************************
  IMPLICIT NONE
  CHARACTER(LEN=2) :: array
  CHARACTER(LEN=1) :: RMM
  CHARACTER(LEN=20) :: FOR
  CHARACTER(LEN=150) :: arch

  INTEGER :: N_elec, N_rho, Nod_x, Nod_z, Nod_elec,j, N, l, i, k, subnivel, iz, ic1s, ic2s, ip1s, ip2s, ra2, rb2
  INTEGER :: N_tot, Nx, Nz, ic1, ic2, ip1, ip2, nivel, Niv, mover, KL, KU, INFO, LDAB, Nod_layer, L_B
  REAL(kind=4) :: D_elec, hx, hz, INF, sigm_m, xaux(12), hs, hn, he, hw, R1, R2, PI=ACOS(-1.0), AMP
  REAL(KIND=4) :: xmin, xmax, zmin, zmax, xr, zr, ar, zaux(12), finish, start, EI, xc, inc, zc, sdf
  REAL(kind=4), ALLOCATABLE :: rho(:), x1(:), z1(:), sigm(:,:), x(:), z(:), PHI_P(:), B(:)
  INTEGER, ALLOCATABLE :: AUX(:,:), id(:), IPIV(:)
  REAL(KIND=4), ALLOCATABLE :: A(:,:), M(:,:)

  REAL :: RHO1, RHO2, E, dx1, dx2, dz

CONTAINS
!********************************************************************
!           RUTINA DE LECTURA DE PARAMETROS DEL MEDIO
!********************************************************************
  SUBROUTINE PARAMETROS
    !leemos de un achivo externo leemos los parametros del medio
    CALL SYSTEM('clear')
    WRITE(*,*) 
    WRITE(*,*)'****************************************************'
    WRITE(*,'(A,$)')'  NOMBRE DEL ARCHIVO: '
    READ(*,11) arch
    11 FORMAT(A)
    WRITE(*,*)'****************************************************'
    WRITE(*,*)
    OPEN(1, FILE=TRIM(arch), STATUS='OLD', ACTION='READ')
      READ(1,*) N_elec                    !numero de electrodos
      READ(1,*) D_elec                    !distancia entre electrodos
      READ(1,*) Nod_elec                  !nodos entre electrodos
      READ(1,*) Nod_z                     !nodos en direccion z
      ALLOCATE(z1(Nod_z))                 !alojamos memoria para nodos en z
      READ(1,*) z1(:)                     !profundidades de los nodos en z
      READ(1,*) array                     !tipo de dispositivo electrodico
      READ(1,*) AMP                       !corriente[amp]
      READ(1,*) RMM, Nod_layer            !metodo de regularizacion
      READ(1,*) N_rho                     !numero de resistividades

      ALLOCATE(rho(N_rho), id(N_rho))    
      DO i=1, N_rho
        READ(1,*) id(i), rho(i)           !leemos una resistvidad y asignamos un identificador de 0-9
      END DO

      N=(N_elec-1)*(Nod_elec+1)+1         !calculamos # de bloques de resistividad en un renglon
      WRITE(FOR,*)N-1                     !escribimos en FOR para formato dinamico
      ALLOCATE(AUX(Nod_z-1,N-1), x1(N))
      DO k=1, Nod_z-1
        READ(1,"("// ADJUSTL(FOR)// "I1.1)" ) AUX(k,:)  !leemos los identificadores de las resistividades 
      END DO     
      CLOSE(1) 

    RETURN
  END SUBROUTINE


!********************************************************************
!         RUTINA PARA EL TRAZO DEL MODELO DE CONDUCTIVIDADES
!********************************************************************
  SUBROUTINE plot_model
   
    hx=D_elec/DBLE(Nod_elec+1)          !distancia entre nodos en la zona de interes
    x1=(/(hx*(i-1), i=1, N)/)           !generamos un vector con las distancias sobre el eje x
    
    !se genera un archivo con parametros para el trazo de la grafica del modelo
    OPEN(1, FILE='grafic.par', ACTION='WRITE', STATUS='REPLACE')
      WRITE(1,*)N_rho                                            !numero de resistividades
      WRITE(1,*)rho(:)                                           !resistividades
      WRITE(1,*)MINVAL(x1), MAXVAL(x1), MINVAL(z1), MAXVAL(z1)   !limites del grafico
      WRITE(1,*) Nod_z
      WRITE(1,*) z1(:)
      WRITE(1,*)SIZE(AUX(1,:))*SIZE(AUX(:,1))                    !cantidad de bloques a trazar
      !DEBEMOS GENERAR CADA OBJETO  
      DO k=1, Nod_z-1
         DO i=1, N-1
                 !cada bloque (rectangulo) encuentra limitado por 4 rectas
                 xmin=x1(i)                                      
                 xmax=x1(i+1)
                 zmin=z1(k)
                 zmax=z1(k+1)
             WRITE(1,*)AUX(k,i), xmin, xmax, zmin, zmax           !el grafico se traza con el indentificado
         END DO                                                   !el valor de resistividad se coloca en la
      END DO                                                      !etiqueta de la paleta de colores
      CLOSE(1)
    RETURN
  END SUBROUTINE


!*****************************************************************************
!           RUTINA CONFIGURACION WENNER
!*****************************************************************************
  SUBROUTINE WN
  CALL MATRIX_A
  ALLOCATE(M(N,N), B(N), IPIV(N))

  IF(RMM=='H')THEN
    sigm_m=SUM(sigm)/(SIZE(sigm(1,:))*SIZE(sigm(:,1)))
    sigm=sigm_m
    CALL MATRIX_M(sigm)
  ELSEIF(RMM=='L')THEN
    sigm(1:Nod_layer,:)= SUM(sigm(1:Nod_layer,:))/(Nod_layer*SIZE(sigm(1,:)))
    sigm(Nod_layer+1:Nz-1,:)= SUM(sigm(Nod_layer+1:Nz-1,:))/((Nz-1-Nod_layer)*SIZE(sigm(1,:)))
    CALL MATRIX_M(sigm) !calculamos la matriz M con la distribucion de resistividad de dos estratos
  END IF

  M=M-A
  
  PRINT*, 'SE CREO LA MATRIZ M'

  KL=NX-2
  KU=NX-2
  LDAB=N
  CALL SPARSE
  
  CALL CPU_TIME(start)
  CALL SGBTRF( N, N, KL, KU, A, LDAB, IPIV, INFO )
  PRINT*, 'TERMINE LA DESCOMPOSICION'


  !calculamos el numero maximo de niveles
  IF(MOD(N_elec,3)==0)THEN 
    Niv=(N_elec/3)-1                !si la division es entera 
  ELSE
    Niv=FLOOR(DBLE(N_elec)/3.0)    !sino bajamos a piso
  END IF
  
  !calculamos el numero total de lecturas posibles a realizas con este dispositivo electrodico
  N_tot = Niv*N_elec-(3*(Niv)*(Niv+1))/2 

  !alojamos memoria a un vector auxiliar, que almacene el potencial primario en superficie
  ALLOCATE(PHI_P(Nx-2))

  
  OPEN(1, FILE='aparente.res', ACTION='WRITE', STATUS='REPLACE')
  WRITE(1,*) N_tot
  WRITE(1,*) MINVAL(x1), MAXVAL(x1), 0.0, Niv*3.0*D_elec*0.173
  !iteramos por cada nivel
  DO nivel=1, Niv
    !en cada nivel se va mermando el numero de datos
    DO mover=1, N_elec-3*nivel

      !fijamos indices de posicion electrodica, luego se desplazan barriendo todos los electrodos
      ic1=13                      + (mover-1)*(Nod_elec+1)
      ip1=13+1*(Nod_elec+1)*nivel + (mover-1)*(Nod_elec+1)
      ip2=13+2*(Nod_elec+1)*nivel + (mover-1)*(Nod_elec+1) 
      ic2=13+3*(Nod_elec+1)*nivel + (mover-1)*(Nod_elec+1)

      !generamos el potencial primario
      IF(RMM=='H')THEN
        j=0 
        DO k=1, Nz-1
          DO i=2, Nx-1 
            j=j+1
            R1=SQRT( (x(ic1)-x(i))**2 + (0.0010-z(k))**2 )
            R2=SQRT( (x(ic2)-x(i))**2 + (0.0010-z(k))**2 )
            R1=R1/1000.0
            R2=R2/1000.0
            B(j)=(AMP/(2.0*PI*sigm_m))*(1.0/R1 - 1.0/R2)
            B(j)=B(j)/1000.0
          END DO
        END DO
      ELSEIF(RMM=='L')THEN
        j=0
        DO k=1, Nz-1
          DO i=2, Nx-1 
            j=j+1
            RHO1= 1.0/sigm(1,13)   !resistividad promedio de la capa
            RHO2= 1.0/sigm(Nod_layer+1,13)   !resistividad promedio de la capa
            E=z(Nod_layer+1)
            dx1=ABS(x(ic1)-x(i))
            dx2=ABS(x(ic2)-x(i))
            dz=ABS(z(k)- 0.00010)
            dx1=dx1/1000.0; dx2=dx2/1000.0; dz=dz/1000.0;E=E/1000.0
            IF(i==ic1)dx1=0.000001
            IF(i==ic2)dx2=0.000001
            !print*, dx1, dx2, dz, rho1, rho2, E
            B(j) = (U_LAYER(dx1,dz,E,AMP,RHO1,RHO2) - U_LAYER(dx2,dz,E,AMP,RHO1,RHO2))/1000.0
            !print*, b(j), i; read*, ;
            !if(isnan(b(j)))stop
          END DO
        END DO
      ELSE
        STOP
      END IF

      !guadamos en PHI_P el potencial primario en superficie
      PHI_P = B(1:Nx-2)
      !aplicamos el stencil sobre el potencial primario (ver ec.)
      B=MATMUL(M,B)


      !subrutina de lapack para realizar sustitucion hacia atras y hacia adelante
      CALL SGBTRS( 'No transpose', N, KL, KU, 1, A, LDAB, IPIV, B, N, INFO )

      !al potencial primario en superficie, sumamos el potencial secundario en superficie
      PHI_P = PHI_P + B(1:NX-2)

      !escribimos resultados a un archivo
      xr=(x(ic1)+x(ic2))*0.50                                        !posicion sobre el eje x
      zr=DBLE(nivel)                                                    !nivel de lectura
      ar=2.0*PI*D_elec*DBLE(nivel)*ABS(PHI_P(ip1-1)-PHI_P(ip2-1))/AMP   !resistividad aparente en [Ohm-m]
      ic1=(ic1-12-1)/(Nod_elec+1) + 1
      ic2=(ic2-12-1)/(Nod_elec+1) + 1
      ip1=(ip1-12-1)/(Nod_elec+1) + 1
      ip2=(ip2-12-1)/(Nod_elec+1) + 1
      WRITE(1,*) ic1, ic2, ip1, ip2, xr, zr*3.0*D_elec*0.173, ar   
    END DO
  END DO
  CLOSE(1)
  CALL CPU_TIME(finish)
  WRITE(*,*) 'RESUELTO EN', finish-start,'[segundos]'
  END SUBROUTINE

!***************************************************************************
! SUBRUTINA PARA EL CASO DE WENNER-SCHLUMBERGER
!***************************************************************************
  SUBROUTINE WS
    INTEGER :: N_tot
    REAL :: xr, zr, ar

  CALL MATRIX_A
  ALLOCATE(M(N,N), B(N), IPIV(N))


  IF(RMM=='H')THEN
    sigm_m=SUM(sigm)/(SIZE(sigm(1,:))*SIZE(sigm(:,1)))
    sigm=sigm_m
    CALL MATRIX_M(sigm)
  ELSEIF(RMM=='L')THEN
    sigm(1:Nod_layer,:)= SUM(sigm(1:Nod_layer,:))/(Nod_layer*SIZE(sigm(1,:)))
    sigm(Nod_layer+1:Nz-1,:)= SUM(sigm(Nod_layer+1:Nz-1,:))/((Nz-1-Nod_layer)*SIZE(sigm(1,:)))
    CALL MATRIX_M(sigm) !calculamos la matriz M con la distribucion de resistividad de dos estratos
  END IF
  M=M-A
  PRINT*, 'SE CREO LA MATRIZ M'

  KL=NX-2
  KU=NX-2
  LDAB=N
  CALL SPARSE
  
  CALL CPU_TIME(start)
  CALL SGBTRF( N, N, KL, KU, A, LDAB, IPIV, INFO )
  PRINT*, 'TERMINE LA DESCOMPOSICION'

    N_tot=0
    OPEN(1, FILE='aux.res', STATUS='REPLACE', ACTION='WRITE')
     WRITE(1,*) N_tot
     WRITE(1,*) MINVAL(x1), MAXVAL(x1), 1.D0, DBLE(N_tot)
     !valores iniciales de electrodos   
     nivel=0
     iz=0
     
     DO
       nivel=nivel+1
       DO subnivel=1, 4
       
         !ic1=13                          
         !ip1=13+(2*nivel-1)*(Nod_elec+1)                      +   (Nod_elec+1)*(subnivel-1)                        
         !ip2=13+(2*nivel-1)*(Nod_elec+1) + nivel*(Nod_elec+1) +   (Nod_elec+1)*(subnivel-1)  
         !ic2= ip2 + ip1 - ic1  
         ic1=13
         ip1=13 + (Nod_elec+1)*subnivel                         + (3*(nivel-1))*(Nod_elec+1)
         ip2=13 + (Nod_elec+1)*subnivel   + (Nod_elec+1)*nivel  + (3*(nivel-1))*(Nod_elec+1)
         ic2=13 + (Nod_elec+1)*subnivel*2 + (Nod_elec+1)*nivel  + (3*(nivel-1))*(Nod_elec+1)*2


         ic1s=ic1
         ic2s=ic2
         ip1s=ip1
         ip2s=ip2         
               
         IF(ic2>Nx-12)THEN
           EXIT
         END IF
         mover=0
         iz=iz+1
         !operacion auxiliar
         ic1=ic1-(Nod_elec+1)
         ic2=ic2-(Nod_elec+1)
         ip1=ip1-(Nod_elec+1)
         ip2=ip2-(Nod_elec+1)
         DO
           mover=mover+1
           ic1=ic1+(Nod_elec+1)
           ic2=ic2+(Nod_elec+1)
           ip1=ip1+(Nod_elec+1)
           ip2=ip2+(Nod_elec+1)
	     IF(ic2>Nx-12)EXIT
           !PRINT*, ic1-12, ic2-12, ip1-12, ip2-12, mover; read*,
        !GENERAR EL VECTOR CON POTENCIALES Y RECORRERLO HASTA LLEGAR AL ULTIMO ELECTRODO
      !generamos el potencial primario
      IF(RMM=='H')THEN
        j=0 
        DO k=1, Nz-1
          DO i=2, Nx-1 
            j=j+1
            R1=SQRT( (x(ic1)-x(i))**2 + (0.0010-z(k))**2 )
            R2=SQRT( (x(ic2)-x(i))**2 + (0.0010-z(k))**2 )
            R1=R1/1000.0
            R2=R2/1000.0
            B(j)=(AMP/(2.0*PI*sigm_m))*(1.0/R1 - 1.0/R2)
            B(j)=B(j)/1000.0
          END DO
        END DO
      ELSEIF(RMM=='L')THEN
        j=0
        DO k=1, Nz-1
          DO i=2, Nx-1 
            j=j+1
            RHO1= 1.0/sigm(1,13)   !resistividad promedio de la capa
            RHO2= 1.0/sigm(Nod_layer+1,13)   !resistividad promedio de la capa
            E=z(Nod_layer+1)
            dx1=ABS(x(ic1)-x(i))
            dx2=ABS(x(ic2)-x(i))
            dz=ABS(z(k)- 0.00010)
            dx1=dx1/1000.0; dx2=dx2/1000.0; dz=dz/1000.0;E=E/1000.0
            IF(i==ic1)dx1=0.000001
            IF(i==ic2)dx2=0.000001
            !print*, dx1, dx2, dz, rho1, rho2, E
            B(j) = (U_LAYER(dx1,dz,E,AMP,RHO1,RHO2) - U_LAYER(dx2,dz,E,AMP,RHO1,RHO2))/1000.0
            !print*, b(j), i; read*, ;
            !if(isnan(b(j)))stop
          END DO
        END DO
      ELSE
        STOP
      END IF
           !USAMOS UN VECTOR AUXILIAR DE DIMENSION SEGUN EL NUMERO DE NODOS EN SUPERFICIE
           PHI_P=B(1:NX-2)

           !MULTIPLICAMOS [M-A]*B !!M  YA CONTIENE LA RESTA M-A
           B=MATMUL(M,B)

           CALL SGBTRS( 'No transpose', N, KL, KU, 1, A, LDAB, IPIV, B, N, INFO )
           PHI_P=PHI_P + B(1:NX-2)
           xr=(x(ic1)+x(ic2))*0.5D0
           zr=iz
           ar= 2.D0*PI*ABS(1.D0/(1.D0/(x(ip1)-x(ic1))-1.D0/(x(ip2)-x(ic1)) &
              &             - 1.D0/(x(ic2)-x(ip1))+1.D0/(x(ic2)-x(ip2))))*ABS(PHI_P(ip1-1)-PHI_P(ip2-1))/AMP
           WRITE(1,*) (ic1-12-1)/(Nod_elec+1) + 1, (ic2-12-1)/(Nod_elec+1) + 1, &
                   &  (ip1-12-1)/(Nod_elec+1) + 1, (ip2-12-1)/(Nod_elec+1) + 1, xr, zr, ar
           N_tot=N_tot+1
         END DO   
       END DO
       IF(ic2s>Nx-12)EXIT
     END DO
     CLOSE(1)

     CALL CPU_TIME(finish)
     WRITE(*,*) 'RESUELTO EN', finish-start,'[segundos]'

     OPEN(1, FILE='aparente.res', STATUS='REPLACE')
     OPEN(2, FILE='aux.res', STATUS='OLD')
       WRITE(1,*) N_tot
       WRITE(1,*) MINVAL(x1), MAXVAL(x1), 1.D0, DBLE(iz)
       READ(2,*)
       READ(2,*)
     DO mover=1, N_tot
        READ(2,*) ic1, ic2, ip1, ip2, xr, zr, ar
        WRITE(1,*)ic1, ic2, ip1, ip2, xr, zr, ar
     END DO
     CLOSE(1)
     CLOSE(2, STATUS='DELETE')

    RETURN
  END SUBROUTINE

!***************************************************************************
! SUBRUTINA PARA EL CASO DIPOLO-DIPOLO
!***************************************************************************
  SUBROUTINE DD
    INTEGER :: N_tot
    REAL :: xr, zr, ar

  CALL MATRIX_A
  ALLOCATE(M(N,N), B(N), IPIV(N))


  IF(RMM=='H')THEN
    sigm_m=SUM(sigm)/(SIZE(sigm(1,:))*SIZE(sigm(:,1)))
    sigm=sigm_m
    CALL MATRIX_M(sigm)
  ELSEIF(RMM=='L')THEN
    sigm(1:Nod_layer,:)= SUM(sigm(1:Nod_layer,:))/(Nod_layer*SIZE(sigm(1,:)))
    sigm(Nod_layer+1:Nz-1,:)= SUM(sigm(Nod_layer+1:Nz-1,:))/((Nz-1-Nod_layer)*SIZE(sigm(1,:)))
    CALL MATRIX_M(sigm) !calculamos la matriz M con la distribucion de resistividad de dos estratos
  END IF
  M=M-A
  PRINT*, 'SE CREO LA MATRIZ M'

  KL=NX-2
  KU=NX-2
  LDAB=N
  CALL SPARSE
  
  CALL CPU_TIME(start)
  CALL SGBTRF( N, N, KL, KU, A, LDAB, IPIV, INFO )
  PRINT*, 'TERMINE LA DESCOMPOSICION'

    N_tot=0
    OPEN(1, FILE='aux.res', STATUS='REPLACE', ACTION='WRITE')
     WRITE(1,*) N_tot
     WRITE(1,*) MINVAL(x1), MAXVAL(x1), 1.D0, DBLE(N_tot)
     !valores iniciales de electrodos   
     nivel=0
     iz=0
     
     DO
       nivel=nivel+1
       DO subnivel=1, 4
         !ic1=13                          
         !ip1=13+(2*nivel-1)*(Nod_elec+1)                      +   (Nod_elec+1)*(subnivel-1)                        
         !ip2=13+(2*nivel-1)*(Nod_elec+1) + nivel*(Nod_elec+1) +   (Nod_elec+1)*(subnivel-1)  
         !ic2= ip2 + ip1 - ic1  
         ic1=13
         ic2=13 + (Nod_elec+1)*nivel
         ip1=13 + (Nod_elec+1)*nivel   +  (2*nivel-1)*(Nod_elec+1) + (subnivel-1)*(Nod_elec+1)
         ip2=13 + (Nod_elec+1)*nivel*2 +  (2*nivel-1)*(Nod_elec+1) + (subnivel-1)*(Nod_elec+1)


         ic1s=ic1
         ic2s=ic2
         ip1s=ip1
         ip2s=ip2         
               
         IF(ip2>Nx-12)THEN
           EXIT
         END IF
         mover=0
         iz=iz+1
         !operacion auxiliar
         ic1=ic1-(Nod_elec+1)
         ic2=ic2-(Nod_elec+1)
         ip1=ip1-(Nod_elec+1)
         ip2=ip2-(Nod_elec+1)
         DO
           mover=mover+1
           ic1=ic1+(Nod_elec+1)
           ic2=ic2+(Nod_elec+1)
           ip1=ip1+(Nod_elec+1)
           ip2=ip2+(Nod_elec+1)
	     IF(ip2>Nx-12)EXIT
           !PRINT*, ic1-12, ic2-12, ip1-12, ip2-12, mover; read*,
        !GENERAR EL VECTOR CON POTENCIALES Y RECORRERLO HASTA LLEGAR AL ULTIMO ELECTRODO
      !generamos el potencial primario
      IF(RMM=='H')THEN
        j=0 
        DO k=1, Nz-1
          DO i=2, Nx-1 
            j=j+1
            R1=SQRT( (x(ic1)-x(i))**2 + (0.0010-z(k))**2 )
            R2=SQRT( (x(ic2)-x(i))**2 + (0.0010-z(k))**2 )
            R1=R1/1000.0
            R2=R2/1000.0
            B(j)=(AMP/(2.0*PI*sigm_m))*(1.0/R1 - 1.0/R2)
            B(j)=B(j)/1000.0
          END DO
        END DO
      ELSEIF(RMM=='L')THEN
        j=0
        DO k=1, Nz-1
          DO i=2, Nx-1 
            j=j+1
            RHO1= 1.0/sigm(1,13)   !resistividad promedio de la capa
            RHO2= 1.0/sigm(Nod_layer+1,13)   !resistividad promedio de la capa
            E=z(Nod_layer+1)
            dx1=ABS(x(ic1)-x(i))
            dx2=ABS(x(ic2)-x(i))
            dz=ABS(z(k)- 0.00010)
            dx1=dx1/1000.0; dx2=dx2/1000.0; dz=dz/1000.0;E=E/1000.0
            IF(i==ic1)dx1=0.000001
            IF(i==ic2)dx2=0.000001
            !print*, dx1, dx2, dz, rho1, rho2, E
            B(j) = (U_LAYER(dx1,dz,E,AMP,RHO1,RHO2) - U_LAYER(dx2,dz,E,AMP,RHO1,RHO2))/1000.0
            !print*, b(j), i; read*, ;
            !if(isnan(b(j)))stop
          END DO
        END DO
      ELSE
        STOP
      END IF
           !USAMOS UN VECTOR AUXILIAR DE DIMENSION SEGUN EL NUMERO DE NODOS EN SUPERFICIE
           PHI_P=B(1:NX-2)

           !MULTIPLICAMOS [M-A]*B !!M  YA CONTIENE LA RESTA M-A
           B=MATMUL(M,B)

           CALL SGBTRS( 'No transpose', N, KL, KU, 1, A, LDAB, IPIV, B, N, INFO )
           PHI_P=PHI_P + B(1:NX-2)
           xr=(x(ic1)+x(ip2))*0.5D0
           zr=iz
           ar= 2.D0*PI*ABS(1.D0/(1.D0/(x(ip1)-x(ic1))-1.D0/(x(ip2)-x(ic1)) &
              &             - 1.D0/(x(ip1)-x(ic2))+1.D0/(x(ip2)-x(ic2))))*ABS(PHI_P(ip1-1)-PHI_P(ip2-1))/AMP
           WRITE(1,*) (ic1-12-1)/(Nod_elec+1) + 1, (ic2-12-1)/(Nod_elec+1) + 1, &
                   &  (ip1-12-1)/(Nod_elec+1) + 1, (ip2-12-1)/(Nod_elec+1) + 1, xr, zr, ar
           N_tot=N_tot+1
         END DO   
       END DO
       IF(ip2s>Nx-12)EXIT
     END DO
     CLOSE(1)

     CALL CPU_TIME(finish)
     WRITE(*,*) 'RESUELTO EN', finish-start,'[segundos]'

     OPEN(1, FILE='aparente.res', STATUS='REPLACE')
     OPEN(2, FILE='aux.res', STATUS='OLD')
       WRITE(1,*) N_tot
       WRITE(1,*) MINVAL(x1), MAXVAL(x1), 1.D0, DBLE(iz)
       READ(2,*)
       READ(2,*)
     DO mover=1, N_tot
        READ(2,*) ic1, ic2, ip1, ip2, xr, zr, ar
        WRITE(1,*)ic1, ic2, ip1, ip2, xr, zr, ar
     END DO
     CLOSE(1)
     CLOSE(2, STATUS='DELETE')

    RETURN
  END SUBROUTINE

!********************************************************************
!   RUTINA PARA CREAR LA MATRIZ DEL SISTEMA, POTENCIAL SECUNDARIO
!********************************************************************
  SUBROUTINE MATRIX_A
  !1.- CONSTRUIR LA MATRIZ DE CONDUCTIVIDADES
  ALLOCATE(sigm(Nod_z+11,N+23), x(N+24), z(Nod_z+12))  !alojamos memoria para la matriz de conductividades

  DO j=1, N_rho
    DO k=1, Nod_z-1
      DO i=1, N-1
        IF(AUX(k,i)==id(j)) sigm(k,12+i)=1.0/rho(j)
      END DO  
    END DO
  END DO  

  !extender el modelo en direccion -x y +x

  DO i=1, 12
    sigm(1:Nod_z-1,i)=sigm(1:Nod_z-1,13)
    sigm(1:Nod_z-1,N+11+i)=sigm(1:Nod_z-1,N+11)
  END DO

  !extender el modelo en direccion +z
  DO k=1, 12
     sigm(Nod_z-1+k,:)=sigm(Nod_z-1,:)
  END DO

  !generamos las distancias en direccion x, 12 nodos para extender hasta +1000 y 12 para -1000
  EI=1.0*(N-1)*hx
  inc=((2000.0-EI)/2 - hx)*(1.0/(66.0*hx))
  xaux(1)=hx
  xc=hx
  DO i=2, 12
    xc = (i-1)*inc*hx + xc
    xaux(i)=xc
  END DO
  x(1:12)=-1.D0*xaux(12:1:-1)

  DO i=1, N
    x(12+i)=DBLE(i-1)*hx
  END DO
  x(N+13:N+24)=xaux
  x(N+13:N+24)=x(N+13:N+24)+x(N+12)  

  !generamos las distancias en direccion z, 12 nodos para extender hasta +1000
  EI=z1(Nod_z)
  hz=z1(Nod_z)-z1(Nod_z-1)
  inc=(1000.0- EI - hz)*(1.0/(66.0*hz))

  zaux(1)=hz
  zc=hz
  DO i=2, 12
    zc = zc + (i-1)*inc*hz
    zaux(i)=zc
  END DO
  
  z(1:Nod_z)=z1
  z(Nod_z+1:Nod_z+12)=zaux+z1(Nod_z)

!  x=(1.0/hx)*x
!  z=(1.0/hx)*z


  !4.- ALOJAMOS MEMORIA A LOS VECTORES A
  Nx=N+24
  Nz=Nod_z+12
  N=(Nx-2)*(Nz-1)
  PRINT*,'NUMERO DE NODOS=', N
  ALLOCATE(A(N,N))
  A=0.0

IF(RMM=='H')THEN
  !GENERAMOS LA MATRIZ A, CICLO PARA GENERAR LOS COEFICIENTES DE LOS NODOS EN LAS ESQUINAS 
  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(j==1)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        ra2= (x(i)-   x(13))**2 + ( z(k) -z(1))**2
        rb2= (x(i)-x(Nx-12))**2 + ( z(k) -z(1))**2
        A(j,j)= -( 0.5*sigm(k,i)*(hs/he) + &
                &  0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs) + &
                &  0.5*sigm(k,i-1)*hs* ( ABS(x(i)-x(13))/ra2 + ABS(x(i)-x(Nx-12))/rb2) )
        A(j,j+1)= 0.5*sigm(k,i)*(hs/he)
        A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)
      ELSEIF(j==(Nx-2))THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        ra2= (x(i)-   x(13))**2 + ( z(k) -z(1))**2
        rb2= (x(i)-x(Nx-12))**2 + ( z(k) -z(1))**2
        A(j,j)= -( 0.5*sigm(k,i)*(hs/he) + &
                &  0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs) + &
                &  0.5*sigm(k,i)*hs* ( ABS(x(i)-x(13))/ra2 + ABS(x(i)-x(Nx-12))/rb2) )
        A(j,j-1)= 0.5*sigm(k,i-1)*(hs/hw)
        A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)        
      ELSEIF(j==(N-Nx+3))THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        ra2= (x(i)-   x(13))**2 + ( z(k) -z(1))**2
        rb2= (x(i)-x(Nx-12))**2 + ( z(k) -z(1))**2
        A(j,j)= -( 0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn) + &
                 & 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he) + & 
                 &(0.5*sigm(k,i-1)*hw + 0.5*sigm(k,i)*he) * ( ABS(z(k)-z(1))/ra2 + ABS(z(k)-z(1))/rb2) + &
                 &(0.5*sigm(k-1,i-1)*hn + 0.5*sigm(k,i-1)*hs) * ( ABS(x(i)-x(13))/ra2 + ABS(x(i)-x(Nx-12))/rb2) )
        A(j,j+1)= 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he)
        A(j,j-(Nx-2))= 0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn)        
      ELSEIF(j==N)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        ra2= (x(i)-   x(13))**2 + ( z(k) -z(1))**2
        rb2= (x(i)-x(Nx-12))**2 + ( z(k) -z(1))**2
        A(j,j)= -( 0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn) + &
                 & 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he) + & 
                 &(0.5*sigm(k,i-1)*hw + 0.5*sigm(k,i)*he) * ( ABS(z(k)-z(1))/ra2 + ABS(z(k)-z(1))/rb2) + &
                 &(0.5*sigm(k,i)*hs + 0.5*sigm(k-1,i)*hn) * ( ABS(x(i)-x(13))/ra2 + ABS(x(i)-x(Nx-12))/rb2) )
        A(j,j-1)= 0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw)
        A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs) 
      ELSE
        CYCLE       
      END IF
    END DO
  END DO  

  !CICLO PARA GENERAR LOS COEFICIENTES EN LAS Z=0, Z=INF, X=-INF Y X=INF, EXEPTO EN LAS ESQUINAS
  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k==1 .AND.  (i/=2 .AND. i/=(Nx-1)))THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        A(j,j)= -( 0.5*sigm(k,i)*(hs/he) + &
                &  0.5*sigm(k,i-1)*(hs/hw) + &
                &  0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs) )
        A(j,j+1)= 0.5*sigm(k,i)*(hs/he) 
        A(j,j-1)= 0.5*sigm(k,i-1)*(hs/hw)
        A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)
      ELSEIF(k==(Nz-1) .AND.  (i/=2 .AND. i/=(Nx-1)))THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        ra2= (x(i)-   x(13))**2 + ( z(k) -z(1))**2
        rb2= (x(i)-x(Nx-12))**2 + ( z(k) -z(1))**2
        A(j,j)= -( 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he) + &
                &  0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw) + &
                &  0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn) + &
                & (0.5*sigm(k,i-1)*hw + 0.5*sigm(k,i)*he)* ( ABS(z(k)-z(1))/ra2 + ABS(z(k)-z(1))/rb2) )
        A(j,j+1)= 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he)
        A(j,j-1)= 0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw)
        A(j,j-(Nx-2))= 0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn)
      ELSEIF(i==2 .AND.  (k/=1 .AND. k/=(Nz-1)))THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        ra2= (x(i)-   x(13))**2 + ( z(k) -z(1))**2
        rb2= (x(i)-x(Nx-12))**2 + ( z(k) -z(1))**2
        A(j,j)= -( 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he) + &
                &  0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn) + &
                &  0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)   + &
                & (0.5*sigm(k-1,i-1)*hn + 0.5*sigm(k,i-1)*hs) * (ABS(x(i)-x(13))/ra2 +ABS(x(i)-x(Nx-12))/rb2) )
        A(j,j+1)= 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he)
        A(j,j-(Nx-2))=  0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn)
        A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)
      ELSEIF(i==(Nx-1) .AND.  (k/=1 .AND. k/=(Nz-1)))THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        ra2= (x(i)-   x(13))**2 + ( z(k) -z(1))**2
        rb2= (x(i)-x(Nx-12))**2 + ( z(k) -z(1))**2
        A(j,j)= -( 0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw) + &
                &  0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn) + &
                &  0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)   + &
                & (0.5*sigm(k,i)*hs + 0.5*sigm(k-1,i)*hn) * (ABS(x(i)-x(13))/ra2 +ABS(x(i)-x(Nx-12))/rb2) )
        A(j,j-1)= 0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw)
        A(j,j-(Nx-2))=  0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn)
        A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)
      ELSE
        CYCLE
      END IF
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k>1 .AND. k<(Nz-1) .AND. i>2 .AND. i<(Nx-1))THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        A(j,j)= -( 0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw) + &
                 & 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he) + &
                 & 0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn) + &
                 & 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs) )
        A(j,j-1)= 0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw)
        A(j,j+1)= 0.5*sigm(k,i)*(hs/he) + 0.5*sigm(k-1,i)*(hn/he)
        A(j,j-(Nx-2))= 0.5*sigm(k-1,i)*(he/hn) + 0.5*sigm(k-1,i-1)*(hw/hn)
        A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs) 
      ELSE
        CYCLE
      END IF     
    END DO
  END DO

ELSEIF(RMM=='L')THEN
 j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k==1)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        A(j,j)= -( 0.5*sigm(k,i)*(hs/he) + &
              &    0.5*sigm(k,i-1)*(hs/hw) + &
              &    0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs) )
      ELSE
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        A(j,j)= -( 0.5*sigm(k-1,i)*(hn/he) + 0.5*sigm(k,i)*(hs/he)  + &
              &    0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw)  + &
              &    0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)  + &
              &    0.5*sigm(k-i,i-1)*(hw/hn) + 0.5*sigm(k-1,i)*(he/hn) )
      END IF
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(i==(Nx-1))CYCLE
      IF(k==1)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        A(j,j+1)= 0.5*sigm(k,i)*(hs/he)
      ELSE
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        A(j,j+1)= 0.5*sigm(k-1,i)*(hn/he) + 0.5*sigm(k,i)*(hs/he)
      END IF
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(i==2)CYCLE
      IF(k==1)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        A(j,j-1)= 0.5*sigm(k,i-1)*(hs/hw)
      ELSE
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        A(j,j-1)= 0.5*sigm(k-1,i-1)*(hn/hw) + 0.5*sigm(k,i-1)*(hs/hw)
      END IF
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k==1)CYCLE
      he=x(i+1)-x(i)
      hw=x(i)-x(i-1)
      hs=z(k+1)-z(k)
      hn=z(k)-z(k-1)
      A(j,j-(Nx-2))= 0.5*sigm(k-1,i-1)*(hw/hn) + 0.5*sigm(k-1,i)*(he/hn)
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k==(Nz-1))CYCLE
      he=x(i+1)-x(i)
      hw=x(i)-x(i-1)
      hs=z(k+1)-z(k)
      hn=z(k)-z(k-1)
      A(j,j+(Nx-2))= 0.5*sigm(k,i-1)*(hw/hs) + 0.5*sigm(k,i)*(he/hs)
    END DO
  END DO
END IF 
  PRINT*, 'SE CREO LA MATRIZ A'
  RETURN
  END SUBROUTINE



!*****************************************************************
!  SUBRUTINA PARA EL CONTRUIR LA MATRIZ DEL POTENCIAL PRIMARIO
!*****************************************************************
 SUBROUTINE MATRIX_M(COND)
  ! PARAMETROS DE ENTRADA 
  ! COND : CONDUCTIVIDAD DEL 
  REAL(KIND=4) :: COND(:,:)
  M=0.0
  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k==1)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        M(j,j)= -( 0.5*COND(k,i)*(hs/he) + &
              &    0.5*COND(k,i-1)*(hs/hw) + &
              &    0.5*COND(k,i-1)*(hw/hs) + 0.5*COND(k,i)*(he/hs) )
      ELSE
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        M(j,j)= -( 0.5*COND(k-1,i)*(hn/he) + 0.5*COND(k,i)*(hs/he)  + &
              &    0.5*COND(k-1,i-1)*(hn/hw) + 0.5*COND(k,i-1)*(hs/hw)  + &
              &    0.5*COND(k,i-1)*(hw/hs) + 0.5*COND(k,i)*(he/hs)  + &
              &    0.5*COND(k-i,i-1)*(hw/hn) + 0.5*COND(k-1,i)*(he/hn) )
      END IF
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(i==(Nx-1))CYCLE
      IF(k==1)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        M(j,j+1)= 0.5*COND(k,i)*(hs/he)
      ELSE
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        M(j,j+1)= 0.5*COND(k-1,i)*(hn/he) + 0.5*COND(k,i)*(hs/he)
      END IF
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(i==2)CYCLE
      IF(k==1)THEN
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        M(j,j-1)= 0.5*COND(k,i-1)*(hs/hw)
      ELSE
        he=x(i+1)-x(i)
        hw=x(i)-x(i-1)
        hs=z(k+1)-z(k)
        hn=z(k)-z(k-1)
        M(j,j-1)= 0.5*COND(k-1,i-1)*(hn/hw) + 0.5*COND(k,i-1)*(hs/hw)
      END IF
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k==1)CYCLE
      he=x(i+1)-x(i)
      hw=x(i)-x(i-1)
      hs=z(k+1)-z(k)
      hn=z(k)-z(k-1)
      M(j,j-(Nx-2))= 0.5*COND(k-1,i-1)*(hw/hn) + 0.5*COND(k-1,i)*(he/hn)
    END DO
  END DO

  j=0
  DO k=1, Nz-1
    DO i=2, Nx-1
      j=j+1
      IF(k==(Nz-1))CYCLE
      he=x(i+1)-x(i)
      hw=x(i)-x(i-1)
      hs=z(k+1)-z(k)
      hn=z(k)-z(k-1)
      M(j,j+(Nx-2))= 0.5*COND(k,i-1)*(hw/hs) + 0.5*COND(k,i)*(he/hs)
    END DO
  END DO
  RETURN
 END SUBROUTINE

!*****************************************************************************
!  SUBRUTINA PARA TRANFORMAR LA MATRIZ A EN FORMATO DIAGONAL BANDED DE LAPACK
!*****************************************************************************
 SUBROUTINE SPARSE
   INTEGER :: imin, imax, ii, jj
   REAL(KIND=4), ALLOCATABLE :: AB(:,:)
   ALLOCATE(AB(N,N))
   DO jj=1, N
     imax=MAXVAL((/1,jj-KU/)) 
     imin=MINVAL((/N,jj+KL/)) 
     DO ii=imax, imin
       AB(KL+KU+1+ii-jj,jj)=A(ii,jj)
     END DO
   END DO
   A=AB
   !DEALLOCATE(AB)
   RETURN
 END SUBROUTINE

FUNCTION U_LAYER(R,PROF,ESP,CORR,RHO1,RHO2)
IMPLICIT NONE
!   PROGRAMER
! GEORGE SANTIZ       19-MARZO-2016   
! DESCRIPCION - SUBRUTINA PARA EL CALCULO DEL POTENCIAL EN CUALQUIER PUNTO DE UN MEDIO BIDIMENSIONAL
!               COMPUESTO DE UN ESTRATO Y UN SEMIESPACIO, CON RESISTIVIDADES RHO1 Y RHO2 RESPECTIVAMENTE
! PARAMETROS DE ENTRADA
! X    -DISTANCIA EN [METROS] SOBRE EL EJE X RESPECTO DE LA UBICACION DE LA FUENTE
! Z    -DISTANCIA EN [METROS] Z RESPECTO DE LA FUENTE
! E    -ESPESOR DEL ESTRATO [METROS]
! IA   -CORRIENTE INYECTADA AL SUBSUELO EN [AMPERS]
! RHO1 -RESISTIVIDAD DEL ESTRATO EN [OHM-M]
! RHO2 -RESISTIVIDAD DEL SEMIESPACIO EN [OHM-M]
!
! PARAMETROS DE SALIDA
! U_LAYER  -FUNCION QUE RETORNA EL VALOR DEL POTENCIAL 
INTEGER :: NF
REAL :: U_LAYER
REAL :: R, PROF, ESP, CORR, RHO1, RHO2, P12, PI=ACOS(-1.0), RESUL
REAL :: ZZ, EE, K12
COMMON  ZZ, EE, K12

P12=(RHO2-RHO1)/(RHO2+RHO1)
ZZ=PROF
EE=ESP
K12=P12

RESUL=CABS(ZHANKS(0,R,FUN,0.0,NF,1))
U_LAYER=(CORR*RHO1*RESUL)/(2.0*PI)

RETURN
END FUNCTION

COMPLEX FUNCTION FUN(EQ)
REAL :: EQ
REAL ::  ZZ, EE, K12
COMMON ZZ, EE, K12

!PRIMERA CAPA
IF(ZZ<=EE)THEN
  IF((2.0*EQ*EE)>30.0)THEN
    IF((EQ*ZZ)>30.0)THEN
      FUN=0.0*(0.0,0.0)  
    ELSE
      FUN=EXP(-EQ*ZZ)*(1.0,0.0)
    END IF
  ELSE
    IF((EQ*ZZ)>25.0)THEN
      FUN= 0.0*(0.0,0.0)
    ELSE
      FUN=(EXP(-EQ*ZZ) + (K12/(EXP(2.0*EQ*EE)-K12))*(EXP(-EQ*ZZ)+EXP(EQ*ZZ)))*(1.0,0.0)
    END IF
  END IF
!SEMIESPACIO
ELSE
  IF((2.0*EQ*EE)>30.0)THEN
    IF((EQ*ZZ)>30.0)THEN
      FUN=0.0*(0.0,0.0)
    ELSE
      FUN=EXP(-EQ*ZZ)*(1.0,0.0)
    END IF
  ELSE
    IF((EQ*ZZ)>30.0)THEN
      FUN=0.0*(0.0,0.0)
    ELSE
      FUN=EXP(-EQ*ZZ)*( 1.0 + (K12/(EXP(2.0*EQ*EE)-K12))*(1.0+EXP(2.0*EQ*EE)) )*(1.0,0.0)
    END IF
  END IF
END IF

RETURN
END FUNCTION

COMPLEX FUNCTION ZHANKS(ORD,MULT,FUN,TOL,NF,NEW)                               
! C=======================================================================        
! C  COMPLEX HANKEL TRANSFORMS OF ORDER 0 OR 1 FOR RELATED (SAVED) KERNELS        
! C  AND FIXED TRANSFORM ARGUMENT B.GT.0.                                         
! C                                                                               
! C--REF: ANDERSON, W.L., 1979, GEOPHYSICS, VOL. 44, NO. 7, P. 1287-1305.         
! C                                                                               
! C--SUBPROGRAM ZHANKS EVALUATES THE INTEGRAL FROM 0 TO INFINITY OF               
! C  FUN(G) *JN(G *B) *DG, DEFINED AS THE COMPLEX HANKEL TRANSFORM OF                
! C  ORDER N (=0 OR 1) AND TRANSFORM ARGUMENT B.GT.0.  THE METHOD IS BY           
! C  ADAPTIVE DIGITAL FILTERING OF THE COMPLEX KERNEL FUNCTION FUN,               
! C  USING DIRECT AND/OR PREVIOUSLY SAVED KERNEL FUNCTION VALUES.                 
! C                                                                               
! C--PARAMETERS (ALL INPUT, EXCEPT NF)                                            
! C                                                                               
! C     N     = ORDER (=0 OR 1) OF THE HANKEL TRANSFORM TO BE EVALUATED.          
! C     B     = REAL TRANSFORM ARGUMENT B.GT.0.0 OF THE HANKEL TRANSFORM.         
! C             IF NEW=0, B IS ASSUMED EQUAL TO THE LAST B USED WHEN NEW=1        
! C             (SEE PARAMETER NEW AND SUBPROGRAM USAGE BELOW).                   
! C     FUN(G)= EXTERNAL DECLARED COMPLEX FUNCTION NAME (USER SUPPLIED)           
! C             OF A REAL ARGUMENT G.GT.0. THIS REFERENCE MUST BE SUPPLIED        
! C             EVEN WHEN NEW=0, SINCE THE ADAPTIVE CONVOLUTION                   
! C             MAY NEED SOME DIRECT FUNCTION CALLS (E.G. IF TOL REDUCED).        
! C             IF PARAMETERS OTHER THAN G ARE REQUIRED IN FUN, USE COMMON        
! C             IN THE CALLING PROGRAM AND IN SUBPROGRAM FUN.  BOTH               
! C             REAL AND IMAGINARY PARTS OF THE COMPLEX FUNCTION FUN(G)           
! C             MUST BE CONTINUOUS BOUNDED FUNCTIONS FOR G.GT.0.0. FOR A          
! C             REAL FUNCTION F1(G), FUN=CMPLX(F1(G),0.0) MAY BE USED.            
! C             TWO INDEPENDENT REAL-FUNCTIONS F1(G),F2(G) MAY BE                 
! C             INTEGRATED IN PARALLEL BY WRITING FUN=CMPLX(F1(G),F2(G)).         
! C     TOL   = REQUESTED REAL TRUNCATION TOLERANCE ACCEPTED AT THE FILTER        
! C             TAILS FOR ADAPTIVE FILTERING.  A TRUNCATION CRITERION IS          
! C             DEFINED DURING CONVOLUTION IN A FIXED ABSCISSA RANGE AS           
! C             THE MAX. ABSOLUTE CONVOLVED PRODUCT TIMES TOL.  TYPICALLY,        
! C             TOL.LE.0.00001 WOULD GIVE ABOUT .01 PER CENT ACCURACY             
! C             FOR WELL-BEHAVED KERNELS AND MODERATE VALUES OF B.  FOR           
! C             VERY LARGE OR SMALL B, A VERY SMALL TOL SHOULD BE USED.           
! C             IN GENERAL, DECREASING THE TOLERANCE WOULD PRODUCE HIGHER         
! C             ACCURACY IN THE CONVOLUTION SINCE MORE FILTER WEIGHTS ARE         
! C             USED (UNLESS EXPONENT UNDERFLOWS OCCUR IN THE KERNEL              
! C             EVALUATION -- SEE NOTE (1) BELOW).                                
! C             FOR MAXIMUM ACCURACY POSSIBLE, TOL=0.0 MAY BE USED.               
! C     NF    = TOTAL NUMBER OF DIRECT FUN CALLS USED DURING CONVOLUTION          
! C             FOR ANY VALUE OF NEW (NF IS AN OUTPUT PARAMETER).                 
! C             NF IS IN THE RANGE 21.LE.NF.LE.283 WHEN NEW=1.  USUALLY,          
! C             NF IS MUCH LESS THAN 283 (OR 0) WHEN NEW=0.                       
! C     NEW  =1 IS REQUIRED FOR THE VERY FIRST CALL TO ZHANKS, OR IF              
! C             FORCING DIRECT FUNCTION FUN(G) CALLS, E.G., IF USING              
! C             ZHANKS FOR UNRELATED KERNELS.                                     
! C             NEW=1 INITIALIZES COMMON/SAVE/FSAVE(283),GSAVE(283),NSAVE         
! C             FOR NSAVE COMPLEX KERNEL VALUES IN FSAVE AND CORRESPONDING        
! C             REAL ARGUMENTS IN GSAVE FOR THE GIVEN PARAMETER B.                
! C     NEW  =0 TO USE RELATED KERNELS (MODIFIED BY USER) CURRENTLY STORED        
! C             IN COMMON/SAVE/. FUN IS CALLED ONLY IF REQUIRED                   
! C             DURING THE CONVOLUTION.  ADDITIONAL FUNCTION VALUES WHEN          
! C             NEEDED ARE AUTOMATICALLY ADDED TO THE COMMON/SAVE/ BLOCK.         
! C                                                                               
! C      * * * * * * * NOTE THAT IT IS THE USERS RESPONSIBILITY TO MODIFY THE            
! C             COMMON FSAVE() VALUES FOR NEW=0 CALLS, EXTERNALLY IN              
! C             THE USERS CALLING PROGRAM (SEE SUBPROGRAM USAGE BELOW).           
! C                                                                               
! C=======================================================================        
! C--SUBPROGRAM USAGE-- ZHANKS IS CALLED AS FOLLOWS                               
! C     ...                                                                       
! C     COMPLEX Z1,Z2,ZHANKS,FSAVE                                                
! C     COMMON/SAVE/FSAVE(283),GSAVE(283),NSAVE                                   
! C     EXTERNAL ZF1,ZF2                                                          
! C     ...                                                                       
! C     Z1=ZHANKS(N1,B,ZF1,TOL,NF1,1)                                             
! C     DO 1 I=1,NSAVE                                                            
! C  C--MODIFY FSAVE IN COMMON/SAVE/ TO OBTAIN RELATED ZF2 FROM ZF1.              
! C  C--E.G. FSAVE(I)=GSAVE(I) *FSAVE(I) -- FOR RELATION ZF2(G)=G *ZF1(G)           
! C   1 CONTINUE                                                                  
! C     Z2=ZHANKS(N2,B,ZF2,TOL,NF2,0)                                             
! C     ...                                                                       
! C     END                                                                       
! C     COMPLEX FUNCTION ZF1(G)                                                   
! C     ...USER SUPPLIED CODE FOR DIRECT EVALUATION OF ZF1(G), G.GT.0.            
! C     END                                                                       
! C     COMPLEX FUNCTION ZF2(G)                                                   
! C     ...USER SUPPLIED CODE FOR DIRECT EVALUATION OF ZF2(G), G.GT.0.            
! C     END                                                                       
! C=======================================================================        
! C--NOTES                                                                        
! C       (1).  EXP-UNDERFLOW MAY OCCUR IN EXECUTING THIS SUBPROGRAM.             
! C             THIS IS OK PROVIDED THE MACHINE SYSTEM CONDITIONALLY SETS         
! C             EXP-UNDERFLOW TO 0.0.                                             
! C       (2).  ANSI FORTRAN (AMERICAN STANDARD X3.9-1966) IS USED, EXCEPT        
! C             DATA STATEMENTS MAY NEED TO BE CHANGED FOR SOME COMPILERS.        
! C             TO CONVERT ZHANKS TO THE NEW AMERICAN STANDARD FORTRAN            
! C             (X3.9-1978), ADD THE FOLLOWING DECLARATION TO THIS ROUTINE        
! C             SAVE Y1,ISAVE                                                     
! C       (3).  THE FILTER ABSCISSA CORRESPONDING TO EACH FILTER WEIGHT           
! C             IS GENERATED IN DOUBLE-PRECISION (TO REDUCE ROUND-OFF),           
! C             BUT IS USED IN SINGLE-PRECISION IN FUNCTION FUN.                  
! C       (4).  NO CHECKS ARE MADE ON CALLING PARAMETERS (TO SAVE TIME),          
! C             HENCE UNPREDICTABLE RESULTS COULD OCCUR IF ZHANKS                 
! C             IS CALLED INCORRECTLY (OR IF FUN OR COMMON IS IN ERROR).          
! C=======================================================================        
! C {VAX-11/780 VERSION FORTRAN-77 (X3.9-1978); SEE NOTE(2) ABOVE.}:              
!       SAVE Y1,ISAVE   
       INTEGER :: ORD, NF, NEW, MM, NSAVE, NONE, ISAVE, ISAVE0, IT
       REAL :: MULT, TOL, G                                                         
       COMPLEX FUN,C,CMAX,FSAVE, GSAVE   
                                               
       COMMON/SAVE/FSAVE(283),GSAVE(283),NSAVE                                   
       DOUBLE PRECISION E1,ER,Y1,Y                                                
       REAL :: T(2),TMAX(2)                                                    
       REAL :: WT0(283),WA0(76),WB0(76),WC0(76),WD0(55),&
              &  WT1(283),WA1(76),WB1(76),WC1(76),WD1(55)                                 
       EQUIVALENCE (WT0(1),WA0(1)),(WT0(77),WB0(1)),(WT0(153),WC0(1)), &          
              &    (WT0(229),WD0(1)),(WT1(1),WA1(1)),(WT1(77),WB1(1)), &                     
              &    (WT1(153),WC1(1)),(WT1(229),WD1(1))                                      
       EQUIVALENCE (C,T(1)),(CMAX,TMAX(1))                                       
 
!-----E=DEXP(.2D0), ER=1.0D0/E                                                  
       DATA E1/1.221402758160169834D0/,ER/.818730753077981859D0/                
                
!--J0-TRANSFORM FILTER WEIGHT ARRAYS (EQUIVALENT TO WT0 ARRAY)                  
       DATA WA0/&                                                                 
       2.1969101E-11, 4.1201161E-09,-6.1322980E-09, 7.2479291E-09,&              
      -7.9821627E-09, 8.5778983E-09,-9.1157294E-09, 9.6615250E-09,&             
      -1.0207546E-08, 1.0796633E-08,-1.1393033E-08, 1.2049873E-08,&              
      -1.2708789E-08, 1.3446466E-08,-1.4174300E-08, 1.5005577E-08,&             
      -1.5807160E-08, 1.6747136E-08,-1.7625961E-08, 1.8693427E-08,&             
      -1.9650840E-08, 2.0869789E-08,-2.1903555E-08, 2.3305308E-08,&             
      -2.4407377E-08, 2.6033678E-08,-2.7186773E-08, 2.9094334E-08,&              
      -3.0266804E-08, 3.2534013E-08,-3.3672072E-08, 3.6408936E-08,&              
      -3.7425022E-08, 4.0787921E-08,-4.1543242E-08, 4.5756842E-08,&              
      -4.6035233E-08, 5.1425075E-08,-5.0893896E-08, 5.7934897E-08,&              
      -5.6086570E-08, 6.5475248E-08,-6.1539913E-08, 7.4301996E-08,&              
      -6.7117043E-08, 8.4767837E-08,-7.2583120E-08, 9.7366568E-08,&              
      -7.7553611E-08, 1.1279873E-07,-8.1416723E-08, 1.3206914E-07,&              
      -8.3217217E-08, 1.5663185E-07,-8.1482581E-08, 1.8860593E-07,&              
      -7.3963141E-08, 2.3109673E-07,-5.7243707E-08, 2.8867452E-07,&              
      -2.6163525E-08, 3.6808773E-07, 2.7049871E-08, 4.7932617E-07,&              
       1.1407365E-07, 6.3720626E-07, 2.5241961E-07, 8.6373487E-07,&              
       4.6831433E-07, 1.1916346E-06, 8.0099716E-07, 1.6696015E-06,&              
       1.3091334E-06, 2.3701475E-06, 2.0803829E-06, 3.4012978E-06/              
       DATA WB0/&                                                                 
       3.2456774E-06, 4.9240402E-06, 5.0005198E-06, 7.1783540E-06,&              
       7.6367633E-06, 1.0522038E-05, 1.1590021E-05, 1.5488635E-05,&              
       1.7510398E-05, 2.2873836E-05, 2.6368006E-05, 3.3864387E-05,&              
       3.9610390E-05, 5.0230379E-05, 5.9397373E-05, 7.4612122E-05,&              
       8.8951409E-05, 1.1094809E-04, 1.3308026E-04, 1.6511335E-04,&              
       1.9895671E-04, 2.4587195E-04, 2.9728181E-04, 3.6629770E-04,&              
       4.4402013E-04, 5.4589361E-04, 6.6298832E-04, 8.1375348E-04,&              
       9.8971624E-04, 1.2132772E-03, 1.4772052E-03, 1.8092022E-03,&              
       2.2045122E-03, 2.6980811E-03, 3.2895354E-03, 4.0238764E-03,&              
       4.9080203E-03, 6.0010999E-03, 7.3216878E-03, 8.9489225E-03,&              
       1.0919448E-02, 1.3340696E-02, 1.6276399E-02, 1.9873311E-02,&              
       2.4233627E-02, 2.9555699E-02, 3.5990069E-02, 4.3791529E-02,&              
       5.3150319E-02, 6.4341372E-02, 7.7506720E-02, 9.2749987E-02,&              
       1.0980561E-01, 1.2791555E-01, 1.4525830E-01, 1.5820085E-01,&              
       1.6058576E-01, 1.4196085E-01, 8.9781222E-02,-1.0238278E-02,&              
      -1.5083434E-01,-2.9059573E-01,-2.9105437E-01,-3.7973244E-02,&              
      3.8273717E-01, 2.2014118E-01,-4.7342635E-01, 1.9331133E-01,&              
       5.3839527E-02,-1.1909845E-01, 9.9317051E-02,-6.6152628E-02,&              
       4.0703241E-02,-2.4358316E-02, 1.4476533E-02,-8.6198067E-03/              
       DATA WC0/&                                                                
      5.1597053E-03,-3.1074602E-03, 1.8822342E-03,-1.1456545E-03,&              
       7.0004347E-04,-4.2904226E-04, 2.6354444E-04,-1.6215439E-04,&              
       9.9891279E-05,-6.1589037E-05, 3.7996921E-05,-2.3452250E-05,&              
       1.4479572E-05,-8.9417427E-06, 5.5227518E-06,-3.4114252E-06,&              
       2.1074101E-06,-1.3019229E-06, 8.0433617E-07,-4.9693681E-07,&              
       3.0702417E-07,-1.8969219E-07, 1.1720069E-07,-7.2412496E-08,&              
       4.4740283E-08,-2.7643004E-08, 1.7079403E-08,-1.0552634E-08,&              
       6.5200311E-09,-4.0284597E-09, 2.4890232E-09,-1.5378695E-09,&              
       9.5019040E-10,-5.8708696E-10, 3.6273937E-10,-2.2412348E-10,&              
       1.3847792E-10,-8.5560821E-11, 5.2865474E-11,-3.2664392E-11,&              
       2.0182948E-11,-1.2470979E-11, 7.7057678E-12,-4.7611713E-12,&              
       2.9415274E-12,-1.8170081E-12, 1.1221034E-12,-6.9271067E-13,&              
       4.2739744E-13,-2.6344388E-13, 1.6197105E-13,-9.9147443E-14,&              
       6.0487998E-14,-3.6973097E-14, 2.2817964E-14,-1.4315547E-14,&              
       9.1574735E-15,-5.9567236E-15, 3.9209969E-15,-2.5911739E-15,&              
       1.6406939E-15,-8.8248590E-16, 3.0195409E-16, 2.2622634E-17,&              
      -8.0942556E-17,-3.7172363E-17, 1.9299542E-16,-3.3388160E-16,&              
       4.6174116E-16,-5.8627358E-16, 7.2227767E-16,-8.7972941E-16,&              
       1.0211793E-15,-1.0940039E-15, 1.0789555E-15,-9.7089714E-16/              
       DATA WD0/&                                                                 
       7.4110927E-16,-4.1700094E-16, 8.5977184E-17, 1.3396469E-16,&              
      -1.7838410E-16, 4.8975421E-17, 1.9398153E-16,-5.0046989E-16,&              
       8.3280985E-16,-1.1544640E-15, 1.4401527E-15,-1.6637066E-15,&              
       1.7777129E-15,-1.7322187E-15, 1.5247247E-15,-1.1771155E-15,&              
       6.9747910E-16,-1.2088956E-16,-4.8382957E-16, 1.0408292E-15,&              
      -1.5220450E-15, 1.9541597E-15,-2.4107448E-15, 2.9241438E-15,&              
      -3.5176475E-15, 4.2276125E-15,-5.0977851E-15, 6.1428456E-15,&              
      -7.3949962E-15, 8.8597601E-15,-1.0515959E-14, 1.2264584E-14,&              
      -1.3949870E-14, 1.5332490E-14,-1.6146782E-14, 1.6084121E-14,&              
      -1.4962523E-14, 1.2794804E-14,-9.9286701E-15, 6.8825809E-15,&              
      -4.0056107E-15, 1.5965079E-15,-7.2732961E-18,-4.0433218E-16,&              
      -6.5679655E-16, 3.3011866E-15,-7.3545910E-15, 1.2394851E-14,&              
      -1.7947697E-14, 2.3774303E-14,-3.0279168E-14, 3.9252831E-14,&              
      -5.5510504E-14, 9.0505371E-14,-1.7064873E-13/                             
 ! C--END OF J0 FILTER WEIGHTS                                                     
                                                                              
 ! C--J1-TRANSFORM FILTER WEIGHT ARRAYS (EQUIVALENT TO WT1 ARRAY)                  
       DATA WA1/&                                                                 
      -4.2129715E-16, 5.3667031E-15,-7.1183962E-15, 8.9478500E-15,&              
      -1.0767891E-14, 1.2362265E-14,-1.3371129E-14, 1.3284178E-14,&              
      -1.1714302E-14, 8.4134738E-15,-3.7726725E-15,-1.4263879E-15,&              
       6.1279163E-15,-9.1102765E-15, 9.9696405E-15,-9.3649955E-15,&              
       8.6009018E-15,-8.9749846E-15, 1.1153987E-14,-1.4914821E-14,&              
       1.9314024E-14,-2.3172388E-14, 2.5605477E-14,-2.6217555E-14,&              
       2.5057768E-14,-2.2485539E-14, 1.9022752E-14,-1.5198084E-14,&              
       1.1422464E-14,-7.9323958E-15, 4.8421406E-15,-2.1875032E-15,&              
      -3.2177842E-17, 1.8637565E-15,-3.3683643E-15, 4.6132219E-15,&              
      -5.6209538E-15, 6.4192841E-15,-6.8959928E-15, 6.9895792E-15,&              
      -6.5355935E-15, 5.6125163E-15,-4.1453931E-15, 2.6358827E-15,&              
      -9.5104370E-16, 1.4600474E-16, 5.6166519E-16, 8.2899246E-17,&              
       5.0032100E-16, 4.3752205E-16, 2.1052293E-15,-9.5451973E-16,&              
       6.4004437E-15,-2.1926177E-15, 1.1651003E-14, 5.8415433E-16,&              
       1.8044664E-14, 1.0755745E-14, 3.0159022E-14, 3.3506138E-14,&              
       5.8709354E-14, 8.1475200E-14, 1.2530006E-13, 1.8519112E-13,&              
       2.7641786E-13, 4.1330823E-13, 6.1506209E-13, 9.1921659E-13,&              
       1.3698462E-12, 2.0447427E-12, 3.0494477E-12, 4.5501001E-12,&              
       6.7870250E-12, 1.0126237E-11, 1.5104976E-11, 2.2536053E-11/              
       DATA WB1/&                                                                
       3.3617368E-11, 5.0153839E-11, 7.4818173E-11, 1.1161804E-10,&              
       1.6651222E-10, 2.4840923E-10, 3.7058109E-10, 5.5284353E-10,&              
       8.2474468E-10, 1.2303750E-09, 1.8355034E-09, 2.7382502E-09,&              
       4.0849867E-09, 6.0940898E-09, 9.0913020E-09, 1.3562651E-08,&              
       2.0233058E-08, 3.0184244E-08, 4.5029477E-08, 6.7176304E-08,&              
       1.0021488E-07, 1.4950371E-07, 2.2303208E-07, 3.3272689E-07,&              
       4.9636623E-07, 7.4049804E-07, 1.1046805E-06, 1.6480103E-06,&              
       2.4585014E-06, 3.6677163E-06, 5.4714550E-06, 8.1626422E-06,&              
       1.2176782E-05, 1.8166179E-05, 2.7099223E-05, 4.0428804E-05,&              
       6.0307294E-05, 8.9971508E-05, 1.3420195E-04, 2.0021123E-04,&              
       2.9860417E-04, 4.4545291E-04, 6.6423156E-04, 9.9073275E-04,&              
       1.4767050E-03, 2.2016806E-03, 3.2788147E-03, 4.8837292E-03,&              
       7.2596811E-03, 1.0788355E-02, 1.5973323E-02, 2.3612041E-02,&              
       3.4655327E-02, 5.0608141E-02, 7.2827752E-02, 1.0337889E-01,&              
       1.4207357E-01, 1.8821315E-01, 2.2996815E-01, 2.5088500E-01,&              
       2.0334626E-01, 6.0665451E-02,-2.0275683E-01,-3.5772336E-01,&              
      -1.8280529E-01, 4.7014634E-01, 7.2991233E-03,-3.0614594E-01,&              
       2.4781735E-01,-1.1149185E-01, 2.5985386E-02, 1.0850279E-02,&              
      -2.2830217E-02, 2.4644647E-02,-2.2895284E-02, 2.0197032E-02/              
       DATA WC1/&                                                                
      -1.7488968E-02, 1.5057670E-02,-1.2953923E-02, 1.1153254E-02,&              
      -9.6138436E-03, 8.2952090E-03,-7.1628361E-03, 6.1882910E-03,&              
      -5.3482055E-03, 4.6232056E-03,-3.9970542E-03, 3.4560118E-03,&              
      -2.9883670E-03, 2.5840861E-03,-2.2345428E-03, 1.9323046E-03,&              
      -1.6709583E-03, 1.4449655E-03,-1.2495408E-03, 1.0805480E-03,&              
      -9.3441130E-04, 8.0803899E-04,-6.9875784E-04, 6.0425624E-04,&              
      -5.2253532E-04, 4.5186652E-04,-3.9075515E-04, 3.3790861E-04,&              
      -2.9220916E-04, 2.5269019E-04,-2.1851585E-04, 1.8896332E-04,&              
      -1.6340753E-04, 1.4130796E-04,-1.2219719E-04, 1.0567099E-04,&              
      -9.1379828E-05, 7.9021432E-05,-6.8334412E-05, 5.9092726E-05,&              
      -5.1100905E-05, 4.4189914E-05,-3.8213580E-05, 3.3045496E-05,&              
      -2.8576356E-05, 2.4711631E-05,-2.1369580E-05, 1.8479514E-05,&              
      -1.5980307E-05, 1.3819097E-05,-1.1950174E-05, 1.0334008E-05,&              
      -8.9364160E-06, 7.7278366E-06,-6.6827083E-06, 5.7789251E-06,&              
      -4.9973715E-06, 4.3215167E-06,-3.7370660E-06, 3.2316575E-06,&              
      -2.7946015E-06, 2.4166539E-06,-2.0898207E-06, 1.8071890E-06,&              
      -1.5627811E-06, 1.3514274E-06,-1.1686576E-06, 1.0106059E-06,&              
      -8.7392952E-07, 7.5573750E-07,-6.5353002E-07, 5.6514528E-07,&              
      -4.8871388E-07, 4.2261921E-07,-3.6546333E-07, 3.1603732E-07/              
       DATA WD1/&                                                                
      -2.7329579E-07, 2.3633470E-07,-2.0437231E-07, 1.7673258E-07,&              
      -1.5283091E-07, 1.3216174E-07,-1.1428792E-07, 9.8831386E-08,&              
      -8.5465227E-08, 7.3906734E-08,-6.3911437E-08, 5.5267923E-08,&              
      -4.7793376E-08, 4.1329702E-08,-3.5740189E-08, 3.0906612E-08,&              
      -2.6726739E-08, 2.3112160E-08,-1.9986424E-08, 1.7283419E-08,&              
      -1.4945974E-08, 1.2924650E-08,-1.1176694E-08, 9.6651347E-09,&              
      -8.3580023E-09, 7.2276490E-09,-6.2501673E-09, 5.4048822E-09,&              
      -4.6739154E-09, 4.0418061E-09,-3.4951847E-09, 3.0224895E-09,&              
      -2.6137226E-09, 2.2602382E-09,-1.9545596E-09, 1.6902214E-09,&              
      -1.4616324E-09, 1.2639577E-09,-1.0930164E-09, 9.4519327E-10,&              
      -8.1736202E-10, 7.0681930E-10,-6.1122713E-10, 5.2856342E-10,&              
      -4.5707937E-10, 3.9526267E-10,-3.4180569E-10, 2.9557785E-10,&              
      -2.5560176E-10, 2.2103233E-10,-1.9113891E-10, 1.6528994E-10,&              
      -1.4294012E-10, 1.2361991E-10,-8.2740936E-11/ 

                          
 !--END OF J1 FILTER WEIGHTS                                                   
       NONE=0                                                                    
       IF(NEW.EQ.0) GO TO 100                                                    
       NSAVE=0    
                                                               
 !-----INITIALIZE KERNEL ABSCISSA GENERATION FOR GIVEN B                         
       Y1=0.7358852661479794460D0/DBLE(MULT)                                        
  100  ZHANKS=(0.0,0.0)                                                          
       CMAX=(0.0,0.0)                                                            
       NF=0                                                                      
       Y=Y1       
                                                              
 !-----BEGIN RIGHT-SIDE CONVOLUTION AT WEIGHT 131 (EITHER NEW=1 OR 0)            
       MM=1  !ASSIGN 110 TO M                                                           
       IT=131                                                                     
       Y=Y *E1                                                                      
       GO TO 200                                                                 
  110  TMAX(1)=AMAX1(ABS(T(1)),TMAX(1))                         
       TMAX(2)=AMAX1(ABS(T(2)),TMAX(2))                                          
       IT=IT+1                                                                     
       Y=Y *E1                              
       IF(IT .LE. 149) GO TO 200                                                    
       IF(TMAX(1) .EQ. 0.0 .AND. TMAX(2) .EQ. 0.0) NONE=1  
                           
 !-----ESTABLISH TRUNCATION CRITERION (CMAX=CMPLX(TMAX(1),TMAX(2))               
       CMAX=TOL *CMAX                                                             
       MM=2 !ASSIGN 120 TO M                                                           
       GO TO 200                                                                 

 !-----CHECK FOR FILTER TRUNCATION AT RIGHT END                                  
  120  IF(ABS(T(1)).LE.TMAX(1).AND.ABS(T(2)).LE.TMAX(2)) GO TO 130               
       IT=IT+1                                                                     
       Y=Y *E1                                                                     
       IF(IT.LE.283) GO TO 200                                                    
  130  Y=Y1                                                                      

 !-----CONTINUE WITH LEFT-SIDE CONVOLUTION AT WEIGHT 130                         
       MM=3 !ASSIGN 140 TO M                                                           
       IT=130                                                                     
       GO TO 200                                                                 

 !-----CHECK FOR FILTER TRUNCATION AT LEFT END                                   
  140  IF(ABS(T(1)).LE.TMAX(1).AND.ABS(T(2)).LE.TMAX(2).AND.&                     
       NONE.EQ.0) GO TO 190                                                     
       IT=IT-1                                                                     
       Y=Y *ER                                                                    
       IF(IT.GT.0) GO TO 200                                                      

 !-----RETURN WITH ISAVE=1 PRESET FOR POSSIBLE NEW=0 USE.                        
  190  ISAVE=1                                                                   

 !-----NORMALIZE BY B TO ACCOUNT FOR INTEGRATION RANGE CHANGE                    
       ZHANKS=ZHANKS/MULT                                                         
       RETURN     
                                                               
 !-----SAVE/RETRIEVE PSEUDO-SUBROUTINE (CALL FUN ONLY WHEN NECESSARY)            
  200  G=SNGL(Y)                                                                
       IF(NEW) 300,210,300                                                       
  210  IF(ISAVE.GT.NSAVE) GO TO 300                                              
       ISAVE0=ISAVE                                                              
  220  IF(G.EQ.GSAVE(ISAVE)) GO TO 240                                           
       ISAVE=ISAVE+1                                                             
       IF(ISAVE.LE.NSAVE) GO TO 220                                              
       ISAVE=ISAVE0                                                              

 !-----G NOT IN COMMON/SAVE/----- EVALUATE FUN.                                  
       GO TO 300                                                                 

 !-----G FOUND IN COMMON/SAVE/----- USE FSAVE AS GIVEN.                          
  240  C=FSAVE(ISAVE)                                                            
       ISAVE=ISAVE+1 
                                                            
 !-----SWITCH ON ORDER N                                                         
  250  IF(ORD) 270,260,270                                                         
  260  C=C *WT0(IT)                                                                 
       GO TO 280                                                                 
  270  C=C *WT1(IT)                                                                
  280  ZHANKS=ZHANKS+C                                           
       GO TO (110,120,140)MM
                                                     
 !-----DIRECT FUN EVALUATION (AND ADD TO END OF COMMON/SAVE/)                    
  300  NSAVE=NSAVE+1                                                              
       C=FUN(G)                                                                  
       NF=NF+1                                                                   
       FSAVE(NSAVE)=C                                                            
       GSAVE(NSAVE)=G                                                         
       GO TO 250                                                                 
       END FUNCTION

END MODULE resistfm

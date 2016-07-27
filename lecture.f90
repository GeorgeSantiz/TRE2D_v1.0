MODULE rutine
INTEGER :: N_rho, N_data, Nod_z
REAL(KIND=8), ALLOCATABLE :: rho(:), dat(:,:), xr(:), yr(:), ar(:), z(:),  ic1(:), ic2(:), ip1(:), ip2(:)
REAL(KIND=8) :: limits(4)

 CONTAINS
 SUBROUTINE parametros
   OPEN(1, FILE='grafic.par', STATUS='OLD')
     READ(1,*) N_rho; ALLOCATE(rho(N_rho))
     READ(1,*) rho(:)
     READ(1,*) limits(:)
     READ(1,*) Nod_z; ALLOCATE(z(Nod_z))
     READ(1,*) z(:)
     READ(1,*) N_data; ALLOCATE(dat(N_data,5))
     DO i=1, N_data
       READ(1,*) dat(i,:)
     END DO
   CLOSE(1, STATUS='DELETE')

 END SUBROUTINE


 SUBROUTINE lres
   OPEN(1, FILE='aparente.res', STATUS='OLD')
     READ(1,*) N_data
     ALLOCATE(xr(N_data), yr(N_data), ar(N_data),  ic1(N_data), ic2(N_data), ip1(N_data), ip2(N_data))
     READ(1,*) limits(:)
     DO i=1, N_data
       READ(1,*) ic1(i), ic2(i), ip1(i), ip2(i), xr(i), yr(i), ar(i) 
     END DO
   CLOSE(1)
 END SUBROUTINE
END MODULE



INCLUDE 'RESFMOD.f90'

PROGRAM main
USE resistfm

  !leer los parametros del modelo directo
  CALL PARAMETROS
  
  !genera los datos para el trazo de la grafica
  CALL plot_model

  !llama al script python para generar la grafica del modelo de resistividades
  CALL SYSTEM('python graf_model.py')
  WRITE(*,*) 'GENERANDO IMAGEN..........'
  WRITE(*,10) 'VERIFIQUE QUE SEA CORRECTA Y PRESIONE ENTER'
  READ(*,*) 
  10 FORMAT (A,$)

  IF(array=='WN')THEN
    CALL WN
  ELSEIF(array=='DD')THEN
    CALL DD
  ELSEIF(array=='WS')THEN
    CALL WS
  END IF

  CALL SYSTEM('python graf_pseudoseccion.py')
END PROGRAM main

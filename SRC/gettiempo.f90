      subroutine get_tiempo(minuto)

!      This subroutine gets the computer time

      implicit none

!      Timing variables

      integer   year, mes, dia, hora, minute, seg, mseg
      real      minuto
      character date*8
      character time*10
      character tzone*5
      integer   values(8)

      call date_and_time(date,time,tzone,values)
      year    = values(1)
      mes     = values(2)
      dia     = values(3)
      hora    = values(5)
      minute  = values(6)
      seg     = values(7)
      mseg    = values(8)
      minuto  = hora * 60.0 + minute + (seg + mseg/1000.0)/60.0

      return
      end

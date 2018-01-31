      ! BC modifications
      if ( level == 1) then 
        ! K ===
        do j = 1 , njcv 
          do i = 1 , nicv 
            ap ( i , j , 1 )    = ap ( i , j , 1 ) + ab ( i , j , 1 ) 
            ap ( i , j , nkcv ) = ap ( i , j , nkcv ) + at ( i , j , nkcv ) 
          enddo 
        enddo 
        ! I ===
        do j = 1 , njcv 
          do k = 1 , nkcv 
            ap ( 1 , j , k ) = ap ( 1 , j , k ) + aw ( 1 , j , k ) 
            ap ( nicv , j , k ) = ap ( nicv , j , k ) + ae ( nicv , j , k ) 
          enddo 
        enddo 
        ! J ===
        do i = 1 , nicv 
          do k = 1 , nkcv 
            ap ( i , 1 , k ) = ap ( i , 1 , k ) + as ( i , 1 , k ) 
            ap ( i , njcv , k ) = ap ( i , njcv , k ) + an ( i , njcv , k ) 
          enddo 
        enddo 
      endif 
      aw ( 1 , : , : ) = 0 
      ae ( nicv , : , : ) = 0 
      as ( : , 1 , : ) = 0 
      an ( : , njcv , : ) = 0 
      ab ( : , : , 1 ) = 0 
      at ( : , : , nkcv ) = 0 

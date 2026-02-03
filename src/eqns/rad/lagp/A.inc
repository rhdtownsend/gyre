A(1, 1) = -invx
A(2, 1) = invx*(this%V_2*(om2 + 3) + this%V_2)
A(1, 2) = -invx*this%V_g/this%V_2
A(2, 2) = invx*(this%V_2*(this%dV_2/this%V_2 + (this%As - this%U + 3)/ &
      this%V_2) + this%V_g)

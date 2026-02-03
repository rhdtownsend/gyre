A(1, 1) = invx*(this%V_g - 1)
A(2, 1) = invx*(om2 - this%As + this%U)/om2
A(1, 2) = -invx*om2*this%V_g
A(2, 2) = invx*(-om2*(3 - this%U) + om2*(this%As - this%U + 3))/om2

A(1, 1) = invx*(-l_e - this%U + 3)
A(2, 1) = invx*(lambda - this%As*this%U - this%U*this%V_g)
A(1, 2) = invx
A(2, 2) = invx*(-l_e - this%U + 2)

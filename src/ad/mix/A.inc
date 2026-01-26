A(1, 1) = invx*(V_g - l_i - 1)
A(2, 1) = invx*(-As_iso + om2)
A(3, 1) = 0
A(4, 1) = -As*alpha_grv*invx*U
A(1, 2) = invx*(-V_g*alpha_gam + eta)
A(2, 2) = invx*(As - l_i - U + 3)
A(3, 2) = 0
A(4, 2) = -V_g*alpha_grv*invx*U
A(1, 3) = -alpha_grv*eta*invx
A(2, 3) = -alpha_grv*invx*(1 - U)
A(3, 3) = invx*(-alpha_grv*(1 - U) + alpha_grv*(-l_i - U + 3))
A(4, 3) = invx*(alpha_grv*lambda - alpha_grv*(1 - U)*(-l_i - U + 2) + ( &
      alpha_grv*(1 - U) - alpha_grv*(-l_i - U + 3))*(U - 1) - (-As - &
      V_g - U + 3)*U)
A(1, 4) = 0
A(2, 4) = alpha_grv*invx
A(3, 4) = alpha_grv*invx
A(4, 4) = invx*(-alpha_grv*(U - 1) + alpha_grv*(-l_i - U + 2))

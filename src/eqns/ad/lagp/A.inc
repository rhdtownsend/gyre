A(1, 1) = invx*(-V_g*alpha_gam + V_g + eta - l_i - 1)
A(2, 1) = invx*(-(-V_g*alpha_gam + V_g + eta - l_i - 1)*V_2 + (As - &
      As_iso - l_i + om2 - U + 3)*V_2)
A(3, 1) = 0
A(4, 1) = invx*(As*alpha_grv*U + V_g*alpha_grv*U)
A(1, 2) = invx*(-V_g*alpha_gam + eta)/V_2
A(2, 2) = invx*(V_g*alpha_gam - eta + (dV_2/V_2 + (As - l_i - U + 3)/V_2 &
      )*V_2)
A(3, 2) = 0
A(4, 2) = V_g*alpha_grv*invx*U/V_2
A(1, 3) = alpha_grv*eta*invx
A(2, 3) = -alpha_grv*eta*invx*V_2
A(3, 3) = alpha_grv*invx*(-l_i - U + 3)
A(4, 3) = alpha_grv*invx*lambda
A(1, 4) = 0
A(2, 4) = -alpha_grv*invx*V_2
A(3, 4) = alpha_grv*invx
A(4, 4) = alpha_grv*invx*(-l_i - U + 2)

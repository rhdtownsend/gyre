A(1, 1) = invx*(-V_g*alpha_gam + V_g + eta - l_i - 1)
A(2, 1) = invx*(-V_2*(-V_g*alpha_gam + V_g + eta - l_i - 1) + V_2*(As - &
      As_iso - U - l_i + om2 + 3))
A(3, 1) = 0
A(4, 1) = invx*(As*U*alpha_grv + U*V_g*alpha_grv)
A(1, 2) = invx*(-V_g*alpha_gam + eta)/V_2
A(2, 2) = invx*(V_2*(dV_2/V_2 + (As - U - l_i + 3)/V_2) + V_g*alpha_gam &
      - eta)
A(3, 2) = 0
A(4, 2) = U*V_g*alpha_grv*invx/V_2
A(1, 3) = alpha_grv*eta*invx
A(2, 3) = -V_2*alpha_grv*eta*invx
A(3, 3) = alpha_grv*invx*(-U - l_i + 3)
A(4, 3) = alpha_grv*invx*lambda
A(1, 4) = 0
A(2, 4) = -V_2*alpha_grv*invx
A(3, 4) = alpha_grv*invx
A(4, 4) = alpha_grv*invx*(-U - l_i + 2)

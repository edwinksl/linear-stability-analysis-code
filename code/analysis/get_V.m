function V = get_V(c, phi, eta, n, xi_p, rho_s)

  E_Theta = 0;

  delta_phi_eq = 1/n*log(xi_p*(c-rho_s)) + E_Theta;
  V = eta + phi + delta_phi_eq;

end

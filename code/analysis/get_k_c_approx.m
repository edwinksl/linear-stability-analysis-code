function k_c_approx = get_k_c_approx(base, t_indices)

  rho_s = base.rho_s;
  n = base.n;
  alpha = base.alpha;
  gamma = base.gamma;

  t = base.t;
  c_1 = base.c(end,:);  % at x = 1 for all t
  cx_1 = base.cx(end,:);  % at x = 1 for all t
  phix_1 = base.phix(end,:);  % at x = 1 for all t
  eta_1 = base.eta_1;

  if t(1) == 0

    beta_1 = 1 + (rho_s+abs(rho_s))/2;

    c_1(1) = beta_1;
    cx_1(1) = 0;
    phix_1(1) = 0;
    eta_1(1) = 0;

  end

  % At x = 1
  alpha_3_N = -alpha*exp(-alpha*n*eta_1) - (1-alpha)*exp((1-alpha)*n*eta_1);
  k_c_approx = sqrt(1./(alpha_3_N*gamma).*(-alpha_3_N*n.*phix_1 + exp(-alpha*n*eta_1).*cx_1./(c_1-rho_s)));

  k_c_approx = k_c_approx(t_indices);

end

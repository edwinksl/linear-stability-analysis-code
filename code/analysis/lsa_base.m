function base = lsa_base(base_dir, str, bc_index)

  % For brevity, drop subscript 0 that denotes base state

  % bc_index: 1 = constant current

  if bc_index == 1
    subdir = 'applied-current/';
  end

  base_full_dir = [base_dir, subdir];

  domain = csvread([base_full_dir, 'domain', '_', str, '.csv'], 8, 0);
  boundaries = csvread([base_full_dir, 'boundaries', '_', str, '.csv'], 8, 0);

  cathode_str = [base_full_dir, 'cathode', '_', str, '.csv'];
  if exist(cathode_str, 'file') == 2
    cathode = csvread(cathode_str, 8, 0);
    cathode(:,1) = [];  % remove column 1
    N_vars_cathode = 4;
    N_t_cathode = size(cathode,2)/N_vars_cathode;
    delta_t_on = cathode(2,2);
    gamma_dc = cathode(2,3);
    temp_cathode = (N_t_cathode-1)*N_vars_cathode;
    V_applied_current = cathode(2, 1:N_vars_cathode:1+temp_cathode);
    J_pulse = cathode(2, 4:N_vars_cathode:4+temp_cathode);
  end

  x = domain(:,1);
  N = length(x);
  dx = 1/(N-1);

  domain(:,1) = [];  % rmeove column 1
  N_vars_domain = 24;
  N_t_domain = size(domain,2)/N_vars_domain;

  z_p = domain(1,10);
  z = domain(1,11);
  D_p = domain(1,12);
  D = domain(1,13);
  rho_s = domain(1,14);
  n = domain(1,15);
  alpha = domain(1,16);
  gamma = domain(1,17);
  beta_m = domain(1,18);
  beta_D = domain(1,19);
  beta_v = domain(1,20);
  Da = domain(1,21);
  t_s = domain(1,22);
  J_applied = domain(1,23);
  xi_p = domain(1,24);

  temp_domain = (N_t_domain-1)*N_vars_domain;
  t = domain(1, 1:N_vars_domain:1+temp_domain);
  c = domain(:, 3:N_vars_domain:3+temp_domain);
  cx = domain(:, 5:N_vars_domain:5+temp_domain);
  cxx = domain(:, 6:N_vars_domain:6+temp_domain);
  phi = domain(:, 7:N_vars_domain:7+temp_domain);
  phix = domain(:, 8:N_vars_domain:8+temp_domain);
  phixx = domain(:, 9:N_vars_domain:9+temp_domain);

  boundaries(:,1) = [];  % remove column 1
  N_vars_boundaries = 3;
  N_t_boundaries = size(boundaries,2)/N_vars_boundaries;

  temp_boundaries = (N_t_boundaries-1)*N_vars_boundaries;
  j_0_0 = boundaries(1, 1:N_vars_boundaries:1+temp_boundaries);
  j_0_1 = boundaries(2, 1:N_vars_boundaries:1+temp_boundaries);
  eta_0 = boundaries(1, 2:N_vars_boundaries:2+temp_boundaries);
  eta_1 = boundaries(2, 2:N_vars_boundaries:2+temp_boundaries);
  ct_0 = boundaries(1, 3:N_vars_boundaries:3+temp_boundaries);
  ct_1 = boundaries(2, 3:N_vars_boundaries:3+temp_boundaries);

  base = struct('t', t, ...
                'x', x, ...
                'c', c, ...
                'cx', cx, ...
                'cxx', cxx, ...
                'phi', phi, ...
                'phix', phix, ...
                'phixx', phixx, ...
                'z_p', z_p, ...
                'z', z, ...
                'D_p', D_p, ...
                'D', D, ...
                'rho_s', rho_s, ...
                'n', n, ...
                'alpha', alpha, ...
                'gamma', gamma, ...
                'beta_m', beta_m, ...
                'beta_D', beta_D, ...
                'beta_v', beta_v, ...
                'Da', Da, ...
                't_s', t_s, ...
                'J_applied', J_applied, ...
                'xi_p', xi_p, ...
                'j_0_0', j_0_0, ...
                'j_0_1', j_0_1, ...
                'eta_0', eta_0, ...
                'eta_1', eta_1, ...
                'ct_0', ct_0, ...
                'ct_1', ct_1, ...
                'dx', dx);

  if exist(cathode_str, 'file') == 2
    base.delta_t_on = delta_t_on;
    base.gamma_dc = gamma_dc;
    base.V_applied_current = V_applied_current;
    base.J_pulse = J_pulse;
  end

end

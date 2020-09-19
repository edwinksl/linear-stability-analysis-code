function lsa_ps(save_flag)

  % ps = parameter sweeps

  elapsed_time = tic;

  bc_index = 1;
  solver_index = 2;
  % solver_index = 1;
  base_dir = 'datasets/parameter-sweeps/base/';
  lsa_dir = 'datasets/parameter-sweeps/lsa/';

  rho_ss = [-1, -0.75, -0.5, -0.25, -0.05];
  Das = 1;
  J_applieds = 1.5;
  N = 1001;

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);

  L_y = 100;
  L_z = 100;
  n_min = 1;
  % n_min = 11;  % n_min = 1 to n_min = 10 are difficult and tedious to converge across all parameter values and their corresponding k values are too small to be physically relevant
  delta_n = 5;
  k_cutoff = 10;
  k_f_approx = get_k_f_approx(rho_ss, Das, J_applieds, N, base_dir, bc_index);
  n_approx = get_n_approx(k_f_approx, L_y, L_z);
  n_maxs = round(n_approx*1.25);

  % Use same n_max for for all rho_s and Da for each J_applied
  n_max = max(max(n_maxs));
  n_maxs = repmat(n_max, N_rho_s, N_Da);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        n_max = n_maxs(i,j,k);
        ks = get_ks(n_min, n_max, delta_n, L_y, L_z, k_cutoff);

        str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
        base = lsa_base(base_dir, str, bc_index);

        N_t = length(base.t);
        t_indices = 1:N_t;
        % t_indices = N_t;
        % t_indices = N_t - 1;
        % t_indices = 1;
        % t_indices = 2;
        % ks = ks([1, 6]);
        ks = repmat({ks}, N_t, 1);
        lsa(base, str, bc_index, t_indices, ks, solver_index, save_flag, lsa_dir);

      end

    end

  end

  total_elapsed_time = toc(elapsed_time);
  fprintf('Total elapsed time = %g s\n', total_elapsed_time)

end

function lsa_ca(save_flag)

  % ca = convergence analysis

  elapsed_time = tic;

  bc_index = 1;
  solver_indices = 2;
  % solver_indices = [1, 2];
  base_dir = 'datasets/convergence-analysis/base/';
  lsa_dir = 'datasets/convergence-analysis/lsa/';

  rho_ss = [-0.05, 0, 0.05];
  % rho_ss = -0.05;
  % rho_ss = 0;
  % rho_ss = 0.05;
  Das = 1;
  J_applieds = 1.5;
  Ns = [251, 501, 1001, 2001, 4001];
  % Ns = [1001, 2001, 4001];
  % Ns = 501;
  % Ns = 1001;
  % Ns = 2001;
  % Ns = 4001;

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_N = length(Ns);
  N_solver_index = length(solver_indices);

  % L_y = 1;
  % L_z = 1;
  L_y = 100;
  L_z = 100;
  n_min = 1;
  % delta_n = 1;
  delta_n = 5;
  % delta_n = 10;
  k_cutoff = 10;
  k_f_approx = get_k_f_approx(rho_ss, Das, J_applieds, Ns, base_dir, bc_index);
  n_approx = get_n_approx(k_f_approx, L_y, L_z);
  n_maxs = round(n_approx*1.25);

  % Use same n_max for all N
  n_max = max(n_maxs, [], 4);
  n_maxs = repmat(n_max, 1, N_Da, N_J_applied, N_N);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        for l = 1:N_N

          N = Ns(l);

          n_max = n_maxs(i,j,k,l);
          ks = get_ks(n_min, n_max, delta_n, L_y, L_z, k_cutoff);

          str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
          base = lsa_base(base_dir, str, bc_index);

          N_t = length(base.t);
          t_index = N_t;
          % t_index = N_t - 1;
          % t_index = 1;
          % t_index = 2;
          % ks = ks([1, 2, 46, 47]);
          ks = repmat({ks}, 1, 1);

          for m = 1:N_solver_index

            solver_index = solver_indices(m);

            lsa_parfor_k(base, str, bc_index, t_index, ks, solver_index, save_flag, lsa_dir);

          end

        end

      end

    end

  end

  total_elapsed_time = toc(elapsed_time);
  fprintf('Total elapsed time = %g s\n', total_elapsed_time)

end

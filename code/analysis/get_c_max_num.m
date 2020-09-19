function varargout = get_c_max_num(case_flag, save_flag)

  elapsed_time = tic;

  bc_index = 1;
  solver_index = 2;

  if strcmp(case_flag, 'ca')

    base_dir = 'datasets/convergence-analysis/base/';
    lsa_dir = 'datasets/convergence-analysis/lsa/';

    rho_ss = [-0.05, 0, 0.05];
    Das = 1;
    J_applieds = 1.5;
    Ns = [251, 501, 1001, 2001, 4001];

  elseif strcmp(case_flag, 'a')

    base_dir = 'datasets/analysis/base/';
    lsa_dir = 'datasets/analysis/lsa/';

    rho_ss = [-0.05, 0, 0.05];
    Das = [0.1, 1, 10];
    J_applieds = [0.5, 1, 1.5];
    Ns = 1001;

  elseif strcmp(case_flag, 'ps')

    base_dir = 'datasets/parameter-sweeps/base/';
    lsa_dir = 'datasets/parameter-sweeps/lsa/';

    rho_ss = [-1, -0.75, -0.5, -0.25, -0.05];
    Das = 1;
    J_applieds = 1.5;
    Ns = 1001;

  end

  if bc_index == 1
    subdir = 'applied-current/';
  end

  lsa_full_dir = [lsa_dir, subdir];

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_N = length(Ns);

  max_approx_data = load([lsa_full_dir, 'max_approx', '.mat']);
  k_max_approx = max_approx_data.k_max_approx;

  k_c_num = cell(N_rho_s, N_Da, N_J_applied, N_N);
  k_max_num = cell(N_rho_s, N_Da, N_J_applied, N_N);
  omega_max_num = cell(N_rho_s, N_Da, N_J_applied, N_N);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        for l = 1:N_N

          N = Ns(l);

          str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
          base = lsa_base(base_dir, str, bc_index);

          t = base.t;
          N_t = length(t);
          t_indices = 1:N_t;
          k_c_0 = get_k_c_approx(base, t_indices);

          k_c = NaN(N_t, 1);
          options_k_c = optimset('TolX', 1e-6);

          k_max = NaN(N_t, 1);
          omega_max = NaN(N_t, 1);

          k_max_approx_ijkl = k_max_approx{i,j,k,l};

          parfor m = 1:N_t

            t_temp = t;

            t_index = t_indices(m);

            t_m = t_temp(t_index);

            if t_m == 0
              % k_c(m) = fzero(@(ks) obj_func(base, str, bc_index, m, ks, solver_index, save_flag, lsa_dir), [0, k_c_0(m)], options_k_c);
              k_c(m) = 0;
            else
              k_c(m) = fzero(@(ks) obj_func(base, str, bc_index, m, ks, solver_index, save_flag, lsa_dir), k_c_0(m), options_k_c);
            end

            if t_m == 0
              k_max_m = 0;
              omega_max_temp_m = -0;
            else
              [k_max_m, omega_max_temp_m] = fminbnd(@(ks) -obj_func(base, str, bc_index, m, ks, solver_index, save_flag, lsa_dir), 0, k_max_approx_ijkl(m));
            end

            k_max(m) = k_max_m;
            omega_max(m) = -omega_max_temp_m;

          end

          k_c_num{i,j,k,l} = k_c;
          k_max_num{i,j,k,l} = k_max;
          omega_max_num{i,j,k,l} = omega_max;

        end

      end

    end

  end

  varargout{1} = k_c_num;
  varargout{2} = k_max_num;
  varargout{3} = omega_max_num;

  if save_flag == 1
    filename = 'c_max_num';

    if bc_index == 1
      subdir = 'applied-current/';
    end

    lsa_full_dir = [lsa_dir, subdir];

    save([lsa_full_dir, filename, '.mat'], 'k_c_num', 'k_max_num', 'omega_max_num')
  end

  total_elapsed_time = toc(elapsed_time);
  fprintf('Total elapsed time = %g s\n', total_elapsed_time)

end

function y = obj_func(base, str, bc_index, t_indices, ks, solver_index, ~, lsa_dir)

  omega_mu = lsa(base, str, bc_index, t_indices, {ks}, solver_index, 0, lsa_dir);
  y = omega_mu{1};

end

function pulse_analysis(sign_flag, save_flag)

  close all

  bc_index = 1;
  graphics_dir = '../../graphics/';

  if strcmp(sign_flag, 'm')

    base_dir = '/home/edwinksl/Dropbox/git/linear-stability-analysis/code/analysis/datasets/pulse/base/negative/';
    lsa_dir = 'datasets/pulse/lsa/negative/';
    rho_ss = -0.05;
    Das = 1;
    J_applieds = 1.5;  % average J_applied, not peak J_applied
    N = 1001;
    delta_t_ons = 1;  % delta_t_on/t_s, not just delta_t_on
    gamma_dcs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
    tol = 5e-3;

  elseif strcmp(sign_flag, 'z')

    base_dir = '/home/edwinksl/Dropbox/git/linear-stability-analysis/code/analysis/datasets/pulse/base/zero/';
    lsa_dir = 'datasets/pulse/lsa/zero/';
    rho_ss = 0;
    Das = 1;
    J_applieds = 1;  % average J_applied, not peak J_applied
    N = 1001;
    delta_t_ons = 0.0125;  % delta_t_on/t_s, not just delta_t_on
    gamma_dcs = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
    tol = 5e-3;

  end

  if bc_index == 1
    subdir = 'applied-current/';
  end

  lsa_full_dir = [lsa_dir, subdir];

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_delta_t_on = length(delta_t_ons);
  N_gamma_dc = length(gamma_dcs);

  max_approx_data = load([lsa_full_dir, 'max_approx', '.mat']);
  k_max_approx = max_approx_data.k_max_approx;
  omega_max_approx = max_approx_data.omega_max_approx;

  omega_ints = NaN(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);
  k_ints = NaN(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);
  lambda_ints = NaN(N_rho_s, N_Da, N_J_applied, N_delta_t_on, N_gamma_dc);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        for l = 1:N_delta_t_on

          delta_t_on = delta_t_ons(l);

          for m = 1:N_gamma_dc

            gamma_dc = gamma_dcs(m);

            str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N), '_', num2str(delta_t_on), '_', num2str(gamma_dc)];
            base = lsa_base(base_dir, str, bc_index);

            t = base.t;
            t_s = base.t_s;
            rho_s = base.rho_s;
            Da = base.Da;
            J_pulse = base.J_pulse;
            V_applied_current = base.V_applied_current;

            k_max_approx_ijklm = k_max_approx{i,j,k,l,m};
            omega_max_approx_ijklm = omega_max_approx{i,j,k,l,m};

            % Set k_max below specified tolerance and corresponding omega_max to 0
            indices = k_max_approx_ijklm < tol;
            k_max_approx_ijklm(indices) = 0;
            omega_max_approx_ijklm(indices) = 0;

            title_str = ['$$\rho_\textnormal{s}$$ = ', num2str(rho_s, '%.3g'), ', ', '$$\textnormal{Da}$$ = ', num2str(Da, '%.3g'), ', ', '$$\bar{J}_\textnormal{a}$$ = ', num2str(J_applied, '%.3g'), ', ', '$$\Delta t_\textnormal{on}/t_\textnormal{s}$$ = ', num2str(delta_t_on, '%.3g'), ', ', '$$\gamma_\textnormal{dc}$$ = ', num2str(gamma_dc, '%.3g')];  % rename J_applied to J_a in title

            Q_pulse = trapz(t, J_pulse);
            J_pulse_average = Q_pulse/(t(end)-t(1));
            omega_int = trapz(t, omega_max_approx_ijklm);
            k_omega_int = trapz(t, k_max_approx_ijklm.*omega_max_approx_ijklm);
            k_int = k_omega_int/omega_int;
            lambda_int = 2*pi/k_int;
            omega_ints(i,j,k,l,m) = omega_int;
            k_ints(i,j,k,l,m) = k_int;
            lambda_ints(i,j,k,l,m) = lambda_int;

            fprintf('rho_s = %g, Da = %g, J_applied = %g, delta_t_on = %g, gamma_dc = %g\n', rho_s, Da, J_applied, delta_t_on, gamma_dc)
            fprintf('Q_pulse = %g\n', Q_pulse)
            fprintf('J_pulse_average = %g\n', J_pulse_average)
            fprintf('omega_int = %g\n', omega_int)
            fprintf('k_int = %g\n', k_int)
            fprintf('lambda_int = %g\n', lambda_int)

            if save_flag == 1
              figure('Visible', 'off')
            else
              figure
            end

            ax = gca;

            plot(t, J_pulse, 'LineWidth', 2)
            xlabel('$$t$$', 'Interpreter', 'latex')
            ylabel('$$J_\textnormal{a}$$', 'Interpreter', 'latex')
            title(title_str, 'Interpreter', 'latex')
            axis tight
            ax.FontSize = 12;

            if save_flag == 1
              J_pulse_filename = ['J_pulse', '_', str];
              saveas(gcf, [graphics_dir, J_pulse_filename, '.eps'], 'epsc')
            end

            if save_flag == 1
              figure('Visible', 'off')
            else
              figure
            end

            ax = gca;

            plot(t, V_applied_current, 'LineWidth', 2)
            xlabel('$$t$$', 'Interpreter', 'latex')
            ylabel('$$V$$', 'Interpreter', 'latex')
            title(title_str, 'Interpreter', 'latex')
            axis tight
            ax.FontSize = 12;

            if save_flag == 1
              V_applied_current_filename = ['V_applied_current', '_', str];
              saveas(gcf, [graphics_dir, V_applied_current_filename, '.eps'], 'epsc')
            end

            if save_flag == 1
              figure('Visible', 'off')
            else
              figure
            end

            ax = gca;

            plot(t, k_max_approx_ijklm, 'LineWidth', 2)
            xlabel('$$t$$', 'Interpreter', 'latex')
            ylabel('$$k_\textnormal{max}$$', 'Interpreter', 'latex')
            title(title_str, 'Interpreter', 'latex')
            axis tight
            ax.FontSize = 12;

            if save_flag == 1
              k_max_approx_filename = ['k_max_approx', '_', str];
              saveas(gcf, [graphics_dir, k_max_approx_filename, '.eps'], 'epsc')
            end

            if save_flag == 1
              figure('Visible', 'off')
            else
              figure
            end

            ax = gca;
            hold on
            ax.ColorOrderIndex = 2;

            plot(t, omega_max_approx_ijklm, 'LineWidth', 2)
            hold off
            xlabel('$$t$$', 'Interpreter', 'latex')
            ylabel('$$\omega_\textnormal{max}$$', 'Interpreter', 'latex')
            title(title_str, 'Interpreter', 'latex')
            axis tight
            if strcmp(sign_flag, 'z')
              ax.YAxis.Exponent = 0;
            end
            ax.FontSize = 12;

            if save_flag == 1
              omega_max_approx_filename = ['omega_max_approx', '_', str];
              saveas(gcf, [graphics_dir, omega_max_approx_filename, '.eps'], 'epsc')
            end

            fprintf('\n')

          end

          int_str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N), '_', num2str(delta_t_on)];

          int_title_str = ['$$\rho_\textnormal{s}$$ = ', num2str(rho_s, '%.3g'), ', ', '$$\textnormal{Da}$$ = ', num2str(Da, '%.3g'), ', ', '$$\bar{J}_\textnormal{a}$$ = ', num2str(J_applied, '%.3g'), ', ', '$$\Delta t_\textnormal{on}/t_\textnormal{s}$$ = ', num2str(delta_t_on, '%.3g')];  % rename J_applied to J_a in title

          if save_flag == 1
            figure('Visible', 'off')
          else
            figure
          end

          ax = gca;
          hold on
          ax.ColorOrderIndex = 2;

          omega_int_ijkl = squeeze(omega_ints(i,j,k,l,:));
          plot(gamma_dcs, omega_int_ijkl, '-o', 'LineWidth', 2)
          hold off
          xlabel('$$\gamma_\textnormal{dc}$$', 'Interpreter', 'latex')
          ylabel('$$\Omega$$', 'Interpreter', 'latex')
          title(int_title_str, 'Interpreter', 'latex')
          axis tight
          ax.FontSize = 12;

          if save_flag == 1
            omega_int_filename = ['omega_int', '_', int_str];
            saveas(gcf, [graphics_dir, omega_int_filename, '.eps'], 'epsc')
          end

          if save_flag == 1
            figure('Visible', 'off')
          else
            figure
          end

          ax = gca;

          lambda_int_ijkl = squeeze(lambda_ints(i,j,k,l,:));
          plot(gamma_dcs, lambda_int_ijkl, '-o', 'LineWidth', 2)
          xlabel('$$\gamma_\textnormal{dc}$$', 'Interpreter', 'latex')
          ylabel('$$\bar{\lambda}_\textnormal{max}$$', 'Interpreter', 'latex')
          title(int_title_str, 'Interpreter', 'latex')
          axis tight
          ax.FontSize = 12;

          if save_flag == 1
            lambda_int_filename = ['lambda_int', '_', int_str];
            saveas(gcf, [graphics_dir, lambda_int_filename, '.eps'], 'epsc')
          end

        end

      end

    end

  end

end

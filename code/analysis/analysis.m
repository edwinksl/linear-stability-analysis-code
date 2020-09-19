function analysis(save_flag)

  close all

  % Define custom colormap
  addpath('DrosteEffect-BrewerMap-9e40245')
  map_1 = brewermap(9, 'Blues');
  map_2 = brewermap(9, 'Greens');
  map_3 = brewermap(9, 'Reds');
  map = [map_1([5,8],:); map_2([5,8],:); map_3(7,:)];

  bc_index = 1;
  base_dir = 'datasets/analysis/base/';
  lsa_dir = 'datasets/analysis/lsa/';
  graphics_dir = '../../graphics/';

  if bc_index == 1
    subdir = 'applied-current/';
  end

  lsa_full_dir = [lsa_dir, subdir];

  rho_ss = [-0.05, 0, 0.05];
  Das = [0.1, 1, 10];
  J_applieds = [0.5, 1, 1.5];
  N = 1001;

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);

  max_approx_data = load([lsa_full_dir, 'max_approx', '.mat']);
  k_max_approx = max_approx_data.k_max_approx;
  omega_max_approx = max_approx_data.omega_max_approx;

  c_max_num_data = load([lsa_full_dir, 'c_max_num', '.mat']);
  k_c_num = c_max_num_data.k_c_num;
  k_max_num = c_max_num_data.k_max_num;
  omega_max_num = c_max_num_data.omega_max_num;

  [axis_bounds, base_bounds] = get_axis_bounds(rho_ss, Das, J_applieds, N, lsa_full_dir);

  for i = 1:N_rho_s

    rho_s = rho_ss(i);

    for j = 1:N_Da

      Da = Das(j);

      for k = 1:N_J_applied

        J_applied = J_applieds(k);

        str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
        filename = ['lsa', '_', str];
        lsa_data = load([lsa_full_dir, filename, '.mat']);

        base = lsa_data.base;
        ks = lsa_data.ks;
        omega_mu = lsa_data.omega_mu;
        t_indices = lsa_data.t_indices;

        t = base.t;
        t_s = base.t_s;
        rho_s = base.rho_s;
        Da = base.Da;
        J_applied = base.J_applied;

        t_rs = t/t_s;  % rs = rescaled
        t_rs_num = t_rs(t_indices);
        log10_ks = cellfun(@log10, ks, 'UniformOutput', false);

        k_c_approx = get_k_c_approx(base, t_indices);
        k_c_num_ijk = k_c_num{i,j,k};

        k_max_approx_ijk = k_max_approx{i,j,k};
        omega_max_approx_ijk = omega_max_approx{i,j,k};
        k_max_num_ijk = k_max_num{i,j,k};
        omega_max_num_ijk = omega_max_num{i,j,k};

        N_t = length(t_rs_num);
        legends = cell(N_t, 1);
        for t_index = 1:N_t
          legends{t_index} = ['$$t$$ = ', num2str(t_rs_num(t_index), '%.3g'), '$$t_\textnormal{s}$$'];
        end

        title_str = ['$$\rho_\textnormal{s}$$ = ', num2str(rho_s, '%.3g'), ', ', '$$\textnormal{Da}$$ = ', num2str(Da, '%.3g'), ', ', '$$J_\textnormal{a}$$ = ', num2str(J_applied, '%.3g')];  % rename J_applied to J_a in title

        axis_bounds_k = axis_bounds(k,:);
        base_bounds_k = base_bounds(k,:);

        % Extract subsets
        if rho_s < 0 || (rho_s >= 0 && J_applied < 1)
          % subset = [1, 3, 5, 6, 7, 8, 14];
          subset = [3, 4, 6, 8, 14];
        elseif rho_s == 0 && J_applied >= 1
          % subset = [1, 3, 5, 6, 7, 8];
          subset = [3, 4, 6, 8];
        elseif rho_s > 0 && J_applied >= 1
          % subset = [1, 3, 5, 6];
          subset = [3, 4, 6];
        end
        ks_subset = ks(subset);
        t_rs_num_subset = t_rs_num(subset);
        log10_ks_subset = log10_ks(subset);
        omega_mu_subset = omega_mu(subset);
        legends_subset = legends(subset);
        N_t_subset = length(t_rs_num_subset);

        % Plot omega vs. k for all t
        % if save_flag == 1
        %   figure('Visible', 'off')
        % else
        %   figure
        % end

        % for t_index = 1:N_t
        %   hold on
        %   plot(log10_ks{t_index}, omega_mu{t_index}, 'LineWidth', 2)
        %   hold off
        % end
        % xlabel('$$\log_{10}k$$', 'Interpreter', 'latex')
        % ylabel('$$\Re\left(\omega\right)$$', 'Interpreter', 'latex')
        % legend(legends, 'Location', 'best', 'Interpreter', 'latex')
        % title(title_str, 'Interpreter', 'latex')
        % axis tight

        % if save_flag == 1
        %   omega_mu_filename = ['omega_mu', '_', str];
        %   saveas(gcf, [graphics_dir, omega_mu_filename, '.eps'], 'epsc')
        % end

        % Plot omega vs. k for subset of t
        set(groot, 'defaultAxesColorOrder', map)

        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

        k_min = 0;
        k_max = axis_bounds_k(2);
        omega_min = axis_bounds_k(3);
        omega_max = axis_bounds_k(4);

        c_min = 0;
        c_max = base_bounds_k(2);
        phi_min = base_bounds_k(3);
        phi_max = 0;
        E_min = 0;
        E_max = base_bounds_k(6);

        patch([k_min, k_max, k_max, k_min], [0, 0, omega_max, omega_max], 'red', 'FaceAlpha', 0.075, 'LineStyle', 'none', 'HandleVisibility', 'off')
        patch([k_min, k_max, k_max, k_min], [0, 0, omega_min, omega_min], 'blue', 'FaceAlpha', 0.075, 'LineStyle', 'none', 'HandleVisibility', 'off')

        ax = gca;
        for t_subset_index = 1:N_t_subset
          hold on
          % plot(log10_ks_subset{t_subset_index}, omega_mu_subset{t_subset_index}, 'LineWidth', 2)
          plot(ks_subset{t_subset_index}, omega_mu_subset{t_subset_index}, 'LineWidth', 2)
          hold off
        end
        % xlabel('$$\log_{10}k$$', 'Interpreter', 'latex')
        xlabel('$$k$$', 'Interpreter', 'latex')
        ylabel('$$\Re\left(\omega\right)$$', 'Interpreter', 'latex')
        % legend(legends_subset, 'Location', 'best', 'Interpreter', 'latex')
        legend(legends_subset, 'Location', 'southwest', 'Interpreter', 'latex')
        title(title_str, 'Interpreter', 'latex')
        axis([k_min, k_max, omega_min, omega_max])
        ax.FontSize = 14;
        ax.TickDir = 'out';  % tick marks directed inward are not visible in EPS files for unknown reasons

        set(groot, 'defaultAxesColorOrder', 'remove')

        if save_flag == 1
          omega_mu_subset_filename = ['omega_mu_subset', '_', str];
          saveas(gcf, [graphics_dir, omega_mu_subset_filename, '.eps'], 'epsc')
        end

        % Plot k_max vs. t and omega_max vs. t
        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

        ax = gca;
        yyaxis left
        hold on
        plot(t_rs_num, k_max_num_ijk, 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        plot(t_rs_num, k_max_approx_ijk, 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2)
        xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
        ylabel('$$k_\textnormal{max}$$', 'Interpreter', 'latex')
        title(title_str, 'Interpreter', 'latex')

        yyaxis right
        plot(t_rs_num, omega_max_num_ijk, 'LineStyle', '-', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 2)
        plot(t_rs_num, omega_max_approx_ijk, 'LineStyle', '--', 'Marker', 's', 'MarkerSize', 12, 'LineWidth', 2)
        hold off
        ylabel('$$\omega_\textnormal{max}$$', 'Interpreter', 'latex')
        if rho_s >= 0 && J_applied > 1
          max_legend_location = 'northwest';
        else
          max_legend_location = 'southeast';
        end
        legend({'$$k_\textnormal{max}$$ (num.)', '$$k_\textnormal{max}$$ (approx.)', '$$\omega_\textnormal{max}$$ (num.)', '$$\omega_\textnormal{max}$$ (approx.)'}, 'Location', max_legend_location, 'Interpreter', 'latex')
        ax.FontSize = 14;
        ax.Clipping = 'off';

        if save_flag == 1
          max_filename = ['max', '_', str];
          saveas(gcf, [graphics_dir, max_filename, '.eps'], 'epsc')
        end

        % Plot k_c vs. t
        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

        ax = gca;
        hold on
        plot(t_rs_num, k_c_num_ijk, 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        ax.ColorOrderIndex = 1;
        plot(t_rs_num, k_c_approx, 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2)
        hold off
        xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
        ylabel('$$k_\textnormal{c}$$', 'Interpreter', 'latex')
        legend({'$$k_\textnormal{c}$$ (num.)', '$$k_\textnormal{c}$$ (approx.)'}, 'Location', 'southeast', 'Interpreter', 'latex')
        title(title_str, 'Interpreter', 'latex')
        ax.FontSize = 14;

        if save_flag == 1
          k_c_filename = ['k_c', '_', str];
          saveas(gcf, [graphics_dir, k_c_filename, '.eps'], 'epsc')
        end

        % Plot base state for subset of t
        x = base.x;
        c = base.c;
        phi = base.phi;
        phix = base.phix;

        c_subset = c(:, subset);
        phi_subset = phi(:, subset);
        phix_subset = phix(:, subset);
        E_subset = -phix_subset;

        base_is_plotted = (Da == 1 && J_applied == 1.5);  % logical expression to check for parameter values corresponding to base state plots that are actually shown and used in paper

        % Print parameters to screen
        if base_is_plotted

          z_p = base.z_p;
          z = base.z;
          D_p = base.D_p;
          D = base.D;
          n = base.n;
          alpha = base.alpha;
          gamma = base.gamma;
          beta_m = base.beta_m;
          beta_D = base.beta_D;
          beta_v = base.beta_v;
          xi_p = base.xi_p;

          fprintf('z_p = %d\n', z_p)
          fprintf('z = %d\n', z)
          fprintf('D_p = %g\n', D_p)
          fprintf('D = %g\n', D)
          fprintf('n = %d\n', n)
          fprintf('alpha = %g\n', alpha)
          fprintf('gamma = %.3g\n', gamma)
          % fprintf('beta_m = %g\n', beta_m)
          fprintf('beta_m = %.2e\n', beta_m)
          fprintf('beta_D = %g\n', beta_D)
          % fprintf('beta_v = %g\n', beta_v)
          fprintf('beta_v = %.2e\n', beta_v)
          fprintf('xi_p = %g\n', xi_p)
          fprintf('\n')

        end

        set(groot, 'defaultAxesColorOrder', map)

        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

        ax = gca;
        plot(x, c_subset, 'LineWidth', 2)
        xlabel('$$x$$', 'Interpreter', 'latex')
        % ylabel('$$c^{\left(0\right)}$$', 'Interpreter', 'latex')
        ylabel('$$c_0$$', 'Interpreter', 'latex')
        legend(legends_subset, 'Location', 'southwest', 'Interpreter', 'latex')
        title(title_str, 'Interpreter', 'latex')
        if base_is_plotted
          ylim([c_min, c_max])
        end
        ax.FontSize = 14;

        if save_flag == 1
          c_subset_filename = ['c_subset', '_', str];
          saveas(gcf, [graphics_dir, c_subset_filename, '.eps'], 'epsc')
        end

        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

        ax = gca;
        plot(x, phi_subset, 'LineWidth', 2)
        xlabel('$$x$$', 'Interpreter', 'latex')
        % ylabel('$$\phi^{\left(0\right)}$$', 'Interpreter', 'latex')
        ylabel('$$\phi_0$$', 'Interpreter', 'latex')
        legend(legends_subset, 'Location', 'southwest', 'Interpreter', 'latex')
        title(title_str, 'Interpreter', 'latex')
        if base_is_plotted
          ylim([phi_min, phi_max])
        end
        ax.FontSize = 14;

        if save_flag == 1
          phi_subset_filename = ['phi_subset', '_', str];
          saveas(gcf, [graphics_dir, phi_subset_filename, '.eps'], 'epsc')
        end

        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

        ax = gca;
        plot(x, E_subset, 'LineWidth', 2)
        xlabel('$$x$$', 'Interpreter', 'latex')
        % ylabel('$$E^{\left(0\right)}$$', 'Interpreter', 'latex')
        ylabel('$$E_0$$', 'Interpreter', 'latex')
        legend(legends_subset, 'Location', 'northwest', 'Interpreter', 'latex')
        title(title_str, 'Interpreter', 'latex')
        if base_is_plotted
          ylim([E_min, E_max])
        end
        ax.FontSize = 14;

        if save_flag == 1
          E_subset_filename = ['E_subset', '_', str];
          saveas(gcf, [graphics_dir, E_subset_filename, '.eps'], 'epsc')
        end

        set(groot, 'defaultAxesColorOrder', 'remove')

      end

    end

  end

  % Collapse plots for numerical k_max, omega_max and k_c across rho_s values

  % Define custom colormaps
  map_1 = brewermap(9, 'Blues');
  map_2 = brewermap(9, 'Reds');
  map_k = map_1([4,6,8],:);
  map_omega = map_2([4,6,8],:);

  % line_styles = {'-', '--', ':'};

  for j = 1:N_Da

    Da = Das(j);

    for k = 1:N_J_applied

      J_applied = J_applieds(k);

      k_c_num_jk = k_c_num(:,j,k);
      k_max_num_jk = k_max_num(:,j,k);
      omega_max_num_jk = omega_max_num(:,j,k);

      % Plot k_max vs. t and omega_max vs. t
      if save_flag == 1
        f_max = figure('Visible', 'off');
      else
        f_max = figure;
      end

      % Plot k_c vs. t
      if save_flag == 1
        f_k_c = figure('Visible', 'off');
      else
        f_k_c = figure;
      end

      rho_s_k_c_legend = cell(N_rho_s, 1);
      rho_s_max_legend = cell(2*N_rho_s, 1);

      title_str = ['$$\textnormal{Da}$$ = ', num2str(Da, '%.3g'), ', ', '$$J_\textnormal{a}$$ = ', num2str(J_applied, '%.3g')];  % rename J_applied to J_a in title

      filename_str = [num2str(Da), '_', num2str(J_applied), '_', num2str(N)];

      for i = 1:N_rho_s

        rho_s = rho_ss(i);

        str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
        filename = ['lsa', '_', str];
        lsa_data = load([lsa_full_dir, filename, '.mat']);

        base = lsa_data.base;
        t_indices = lsa_data.t_indices;

        t = base.t;
        t_s = base.t_s;
        rho_s = base.rho_s;
        Da = base.Da;
        J_applied = base.J_applied;

        t_rs = t/t_s;  % rs = rescaled
        t_rs_num = t_rs(t_indices);

        rho_s_k_c_legend{i} = ['$$k_\textnormal{c}$$, $$\rho_\textnormal{s}$$ = ', num2str(rho_s)];
        rho_s_max_legend{i} = ['$$k_\textnormal{max}$$, $$\rho_\textnormal{s}$$ = ', num2str(rho_s)];
        rho_s_max_legend{i+N_rho_s} = ['$$\omega_\textnormal{max}$$, $$\rho_\textnormal{s}$$ = ', num2str(rho_s)];

        % Plot k_max vs. t and omega_max vs. t
        % figure(f_max)
        set(0, 'CurrentFigure', f_max)
        set(f_max, 'defaultAxesColorOrder', [map_k(2,:); map_omega(2,:)])
        yyaxis left
        hold on
        % plot(t_rs_num, k_max_num_jk{i}, 'LineStyle', line_styles{i}, 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        plot(t_rs_num, k_max_num_jk{i}, 'Color', map_k(i,:), 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        hold off
        yyaxis right
        hold on
        % plot(t_rs_num, omega_max_num_jk{i}, 'LineStyle', line_styles{i}, 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 2)
        plot(t_rs_num, omega_max_num_jk{i}, 'Color', map_omega(i,:), 'LineStyle', '-', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 2)
        hold off

        % Plot k_c vs. t
        % figure(f_k_c)
        set(0, 'CurrentFigure', f_k_c)
        ax_k_c = gca;
        ax_k_c.ColorOrderIndex = 1;
        hold on
        % plot(t_rs_num, k_c_num_jk{i}, 'LineStyle', line_styles{i}, 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        plot(t_rs_num, k_c_num_jk{i}, 'Color', map_k(i,:), 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        hold off

      end

      % figure(f_max)
      set(0, 'CurrentFigure', f_max)
      ax_max = gca;
      xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
      yyaxis left
      ylabel('$$k_\textnormal{max}$$', 'Interpreter', 'latex')
      yyaxis right
      ylabel('$$\omega_\textnormal{max}$$', 'Interpreter', 'latex')
      if J_applied > 1
        max_legend_location = 'northeast';
      else
        max_legend_location = 'southeast';
      end
      legend(rho_s_max_legend, 'Location', max_legend_location, 'Interpreter', 'latex')
      title(title_str, 'Interpreter', 'latex')
      ax_max.FontSize = 14;
      ax_max.Clipping = 'off';

      if save_flag == 1
        rho_s_max_filename = ['rho_s_max', '_', filename_str];
        saveas(gcf, [graphics_dir, rho_s_max_filename, '.eps'], 'epsc')
      end

      % figure(f_k_c)
      set(0, 'CurrentFigure', f_k_c)
      ax_k_c = gca;
      xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
      ylabel('$$k_\textnormal{c}$$', 'Interpreter', 'latex')
      if J_applied > 1
        k_c_legend_location = 'northeast';
      else
        k_c_legend_location = 'southeast';
      end
      legend(rho_s_k_c_legend, 'Location', k_c_legend_location, 'Interpreter', 'latex')
      title(title_str, 'Interpreter', 'latex')
      ax_k_c.FontSize = 14;

      if save_flag == 1
        rho_s_k_c_filename = ['rho_s_k_c', '_', filename_str];
        saveas(gcf, [graphics_dir, rho_s_k_c_filename, '.eps'], 'epsc')
      end

    end

  end

end

function [axis_bounds, base_bounds] = get_axis_bounds(rho_ss, Das, J_applieds, N, lsa_full_dir)

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);

  axis_bounds = NaN(N_J_applied, 4);
  base_bounds = NaN(N_J_applied, 6);

  for i = 1:N_J_applied

    J_applied = J_applieds(i);

    axis_bounds_temp = NaN(N_rho_s, N_Da, 4);
    base_bounds_temp = NaN(N_rho_s, N_Da, 6);

    for j = 1:N_rho_s

      rho_s = rho_ss(j);

      for k = 1:N_Da

        Da = Das(k);

        str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
        filename = ['lsa', '_', str];
        lsa_data = load([lsa_full_dir, filename, '.mat']);

        base = lsa_data.base;
        ks = lsa_data.ks;
        omega_mu = lsa_data.omega_mu;

        ks_ijk = horzcat(ks{:});
        omega_real_ijk = real(vertcat(omega_mu{:}));
        k_min = min(ks_ijk);
        k_max = max(ks_ijk);
        omega_min = min(omega_real_ijk);
        omega_max = max(omega_real_ijk);

        % Manually set k_max, omega_min and omega_max
        if J_applied == 0.5
          k_max = 250;
          omega_min = -6.4e-4;
          omega_max = 1.6e-3;
        elseif J_applied == 1
          k_max = 500;
          omega_min = -6e-3;
          omega_max = 0.015;
        elseif J_applied == 1.5
          k_max = 1000;
          omega_min = -0.055;
          omega_max = 0.11;
        end

        axis_bounds_temp(j,k,:) = [k_min; k_max; omega_min; omega_max];

        c = base.c;
        phi = base.phi;
        phix = base.phix;
        E = -phix;

        c_linear = c(:);
        phi_linear = phi(:);
        E_linear = E(:);

        c_min = min(c_linear);
        c_max = max(c_linear);
        phi_min = min(phi_linear);
        phi_max = max(phi_linear);
        E_min = min(E_linear);
        E_max = max(E_linear);

        base_bounds_temp(j,k,:) = [c_min; c_max; phi_min; phi_max; E_min; E_max];

      end

    end

    k_min_omega_min = min(min(axis_bounds_temp(:, :, [1,3])));
    k_max_omega_max = max(max(axis_bounds_temp(:, :, [2,4])));
    axis_bounds(i, [1,3]) = k_min_omega_min;
    axis_bounds(i, [2,4]) = k_max_omega_max;

    c_min_phi_min_E_min = min(min(base_bounds_temp(:, :, [1,3,5])));
    c_max_phi_max_E_max = max(max(base_bounds_temp(:, :, [2,4,6])));
    base_bounds(i, [1,3,5]) = c_min_phi_min_E_min;
    base_bounds(i, [2,4,6]) = c_max_phi_max_E_max;

  end

end

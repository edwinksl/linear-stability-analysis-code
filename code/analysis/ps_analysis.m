function ps_analysis(save_flag)

  close all

  % Define custom colormap
  addpath('DrosteEffect-BrewerMap-9e40245')
  map_1 = brewermap(9, 'Blues');
  map_2 = brewermap(9, 'Greens');
  map_3 = brewermap(9, 'Reds');
  map = [map_1([5,8],:); map_2([5,8],:); map_3(7,:)];

  bc_index = 1;
  base_dir = 'datasets/parameter-sweeps/base/';
  lsa_dir = 'datasets/parameter-sweeps/lsa/';
  graphics_dir = '../../graphics/';

  if bc_index == 1
    subdir = 'applied-current/';
  end

  lsa_full_dir = [lsa_dir, subdir];

  rho_ss = [-1, -0.75, -0.5, -0.25, -0.05];
  Das = 1;
  J_applieds = 1.5;
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

  axis_bounds = get_axis_bounds(rho_ss, Das, J_applieds, N, lsa_full_dir);

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
        %   omega_mu_filename = ['ps', '_', 'omega_mu', '_', str];
        %   saveas(gcf, [graphics_dir, omega_mu_filename, '.eps'], 'epsc')
        % end

        % Plot omega vs. k for subset of t
        set(groot, 'defaultAxesColorOrder', map)

        if save_flag == 1
          figure('Visible', 'off')
        else
          figure
        end

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
        axis([0, axis_bounds_k(2:end)])
        ax.FontSize = 14;

        set(groot, 'defaultAxesColorOrder', 'remove')

        % if save_flag == 1
        %   omega_mu_subset_filename = ['ps', '_', 'omega_mu_subset', '_', str];
        %   saveas(gcf, [graphics_dir, omega_mu_subset_filename, '.eps'], 'epsc')
        % end

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

        % if save_flag == 1
        %   max_filename = ['ps', '_', 'max', '_', str];
        %   saveas(gcf, [graphics_dir, max_filename, '.eps'], 'epsc')
        % end

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

        % if save_flag == 1
        %   k_c_filename = ['ps', '_', 'k_c', '_', str];
        %   saveas(gcf, [graphics_dir, k_c_filename, '.eps'], 'epsc')
        % end

      end

    end

  end

  % Collapse plots for numerical k_max, omega_max and k_c across rho_s values

  % Define custom colormaps
  map_1 = brewermap(9, 'Blues');
  map_2 = brewermap(9, 'Reds');
  map_k = map_1([5,6,7,8,9],:);
  map_omega = map_2([5,6,7,8,9],:);

  for j = 1:N_Da

    Da = Das(j);

    for k = 1:N_J_applied

      J_applied = J_applieds(k);

      k_c_num_jk = k_c_num(:,j,k);
      k_max_num_jk = k_max_num(:,j,k);
      omega_max_num_jk = omega_max_num(:,j,k);

      % Plot k_max vs. t and omega_max vs. t
      if save_flag == 1
        f_k_max = figure('Visible', 'off');
        f_omega_max = figure('Visible', 'off');
      else
        f_k_max = figure;
        f_omega_max = figure;
      end

      % Plot k_c vs. t
      if save_flag == 1
        f_k_c = figure('Visible', 'off');
      else
        f_k_c = figure;
      end

      rho_s_legend = cell(N_rho_s, 1);

      title_str = ['$$\textnormal{Da}$$ = ', num2str(Da, '%.3g'), ', ', '$$J_\textnormal{a}$$ = ', num2str(J_applied, '%.3g')];  % rename J_applied to J_a in title

      filename_str = [num2str(Da), '_', num2str(J_applied), '_', num2str(N)];

      set(f_k_max, 'defaultAxesColorOrder', map_k)
      set(f_omega_max, 'defaultAxesColorOrder', map_omega)
      set(f_k_c, 'defaultAxesColorOrder', map_k)

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

        rho_s_legend{i} = ['$$\rho_\textnormal{s}$$ = ', num2str(rho_s)];

        % Plot k_max vs. t and omega_max vs. t
        set(0, 'CurrentFigure', f_k_max)
        hold on
        plot(t_rs_num, k_max_num_jk{i}, 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        hold off

        set(0, 'CurrentFigure', f_omega_max)
        hold on
        plot(t_rs_num, omega_max_num_jk{i}, 'LineStyle', '-', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 2)
        hold off

        % Plot k_c vs. t
        set(0, 'CurrentFigure', f_k_c)
        hold on
        plot(t_rs_num, k_c_num_jk{i}, 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
        hold off

      end

      set(0, 'CurrentFigure', f_k_max)
      ax_k_max = gca;
      xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
      ylabel('$$k_\textnormal{max}$$', 'Interpreter', 'latex')
      legend(rho_s_legend, 'Location', 'northwest', 'Interpreter', 'latex')
      title(title_str, 'Interpreter', 'latex')
      ax_k_max.FontSize = 14;

      if save_flag == 1
        rho_s_k_max_filename = ['ps', '_', 'rho_s_k_max', '_', filename_str];
        saveas(gcf, [graphics_dir, rho_s_k_max_filename, '.eps'], 'epsc')
      end

      set(0, 'CurrentFigure', f_omega_max)
      ax_omega_max = gca;
      xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
      ylabel('$$\omega_\textnormal{max}$$', 'Interpreter', 'latex')
      legend(rho_s_legend, 'Location', 'northwest', 'Interpreter', 'latex')
      title(title_str, 'Interpreter', 'latex')
      ax_omega_max.FontSize = 14;

      if save_flag == 1
        rho_s_omega_max_filename = ['ps', '_', 'rho_s_omega_max', '_', filename_str];
        saveas(gcf, [graphics_dir, rho_s_omega_max_filename, '.eps'], 'epsc')
      end

      set(0, 'CurrentFigure', f_k_c)
      ax_k_c = gca;
      xlabel('$$t/t_\textnormal{s}$$', 'Interpreter', 'latex')
      ylabel('$$k_\textnormal{c}$$', 'Interpreter', 'latex')
      legend(rho_s_legend, 'Location', 'northwest', 'Interpreter', 'latex')
      title(title_str, 'Interpreter', 'latex')
      ax_k_c.FontSize = 14;

      if save_flag == 1
        rho_s_k_c_filename = ['ps', '_', 'rho_s_k_c', '_', filename_str];
        saveas(gcf, [graphics_dir, rho_s_k_c_filename, '.eps'], 'epsc')
      end

      % Plot numerical k_max, omega_max and k_c as functions of rho_s at last t/t_s value

      t_rs_num_last = t_rs_num(end);  % use t_rs_num from last iteration of previous for loop
      k_c_num_last = NaN(N_rho_s, 1);
      k_max_num_last = NaN(N_rho_s, 1);
      omega_max_num_last = NaN(N_rho_s, 1);

      title_last_str = ['$$t/t_\textnormal{s}$$ = ', num2str(t_rs_num_last, '%.3g'), ', ', title_str];

      for i = 1:N_rho_s

        k_c_num_last(i) = k_c_num_jk{i}(end);
        k_max_num_last(i) = k_max_num_jk{i}(end);
        omega_max_num_last(i) = omega_max_num_jk{i}(end);

      end

      if save_flag == 1
        figure('Visible', 'off')
      else
        figure
      end

      ax = gca;
      yyaxis left
      plot(rho_ss, k_max_num_last, 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
      xlabel('$$\rho_\textnormal{s}$$', 'Interpreter', 'latex')
      ylabel('$$k_\textnormal{max}$$', 'Interpreter', 'latex')
      title(title_last_str, 'Interpreter', 'latex')

      yyaxis right
      plot(rho_ss, omega_max_num_last, 'LineStyle', '-', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 2)
      ylabel('$$\omega_\textnormal{max}$$', 'Interpreter', 'latex')
      ax.FontSize = 14;
      ax.Clipping = 'off';

      if save_flag == 1
        figure('Visible', 'off')
      else
        figure
      end

      ax = gca;
      plot(rho_ss, k_c_num_last, 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 16, 'LineWidth', 2)
      xlabel('$$\rho_\textnormal{s}$$', 'Interpreter', 'latex')
      ylabel('$$k_\textnormal{c}$$', 'Interpreter', 'latex')
      title(title_last_str, 'Interpreter', 'latex')
      ax.FontSize = 14;

    end

  end

end

function axis_bounds = get_axis_bounds(rho_ss, Das, J_applieds, N, lsa_full_dir)

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);

  axis_bounds = NaN(N_J_applied, 4);

  for i = 1:N_J_applied

    J_applied = J_applieds(i);

    axis_bounds_temp = NaN(N_rho_s, N_Da, 4);

    for j = 1:N_rho_s

      rho_s = rho_ss(j);

      for k = 1:N_Da

        Da = Das(k);

        str = [num2str(rho_s), '_', num2str(Da), '_', num2str(J_applied), '_', num2str(N)];
        filename = ['lsa', '_', str];
        lsa_data = load([lsa_full_dir, filename, '.mat']);

        ks = lsa_data.ks;
        omega_mu = lsa_data.omega_mu;

        ks_ijk = horzcat(ks{:});
        omega_real_ijk = real(vertcat(omega_mu{:}));
        k_min = min(ks_ijk);
        k_max = max(ks_ijk);
        omega_min = min(omega_real_ijk);
        omega_max = max(omega_real_ijk);

        axis_bounds_temp(j,k,:) = [k_min; k_max; omega_min; omega_max];

      end

    end

    k_min_omega_min = min(axis_bounds_temp(:, :, [1,3]));
    k_max_omega_max = max(axis_bounds_temp(:, :, [2,4]));
    axis_bounds(i, [1,3]) = k_min_omega_min;
    axis_bounds(i, [2,4]) = k_max_omega_max;

  end

end

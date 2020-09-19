function k_f_approx = get_k_f_approx(rho_ss, Das, J_applieds, Ns, base_dir, bc_index)

  N_rho_s = length(rho_ss);
  N_Da = length(Das);
  N_J_applied = length(J_applieds);
  N_N = length(Ns);

  k_f_approx = NaN(N_rho_s, N_Da, N_J_applied, N_N);

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

          N_t = length(base.t);
          t_indices = 1:N_t;
          k_c_approx = get_k_c_approx(base, t_indices);
          k_f_approx(i,j,k,l) = max(k_c_approx);

        end

      end

    end

  end

end

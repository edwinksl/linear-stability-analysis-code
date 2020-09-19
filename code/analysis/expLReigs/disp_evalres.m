function disp_evalres(ritzval,relres)

% display computed Ritz values and corresponding eigenresidual norms

    for ii = 1 : length(ritzval)
        if isreal(ritzval(ii))
            if ritzval(ii) > 0
                fprintf('  +%.4e\t\t\t%.2e\n',ritzval(ii),relres(ii));
            else
                fprintf('  %.4e\t\t\t%.2e\n',ritzval(ii),relres(ii));
            end
        else
            if real(ritzval(ii)) >= 0 && imag(ritzval(ii)) > 0
                fprintf('  +%.4e  +%.4ei\t%.2e\n',real(ritzval(ii)),imag(ritzval(ii)),relres(ii));
            elseif real(ritzval(ii)) >= 0 && imag(ritzval(ii)) < 0
                fprintf('  +%.4e  %.4ei\t%.2e\n',real(ritzval(ii)),imag(ritzval(ii)),relres(ii));
            elseif real(ritzval(ii)) < 0 && imag(ritzval(ii)) > 0
                fprintf('  %.4e  +%.4ei\t%.2e\n',real(ritzval(ii)),imag(ritzval(ii)),relres(ii));
            elseif real(ritzval(ii)) < 0 && imag(ritzval(ii)) < 0
                fprintf('  %.4e  %.4ei\t%.2e\n',real(ritzval(ii)),imag(ritzval(ii)),relres(ii));
            end
        end
    end
end


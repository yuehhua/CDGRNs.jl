function check_confint(coefs, ses, truth)
    checks = [(c - 2s ≤ t ≤ c + 2s) for (c, s, t) in zip(coefs, ses, truth)]
    all(checks)
end

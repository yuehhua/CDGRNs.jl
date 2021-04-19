function check_confint(coefs, ses, truth; z=3)
    checks = [(c - z*s ≤ t ≤ c + z*s) for (c, s, t) in zip(coefs, ses, truth)]
    all(checks)
end

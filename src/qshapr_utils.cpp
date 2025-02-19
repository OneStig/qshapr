#include "qshapr_utils.h"

Eigen::VectorXd inv_binom_coef(int d) {
    Eigen::VectorXd coef(d + 1);
    coef(0) = 1;

    for (int i = 1; i < d / 2 + 1; i++) {
        coef(i) = (coef(i - 1) * (d - i + 1)) / i;
    }

    for (int i = d / 2 + 1; i < d + 1; i++) {
        coef(i) = coef(d - i);
    }

    return coef.cwiseInverse();
}

Eigen::MatrixXcd complex_v_invc_degree(int d) {
    Eigen::VectorXcd omega_inv(d);
    for (int k = 0; k < d; k++) {
        double theta = -2 * M_PI * k / d;
        omega_inv(k) = std::polar(1.0, theta); // e^{i * theta}
    }

    Eigen::MatrixXcd v_omega_inv(d, d);

    for (int i = 0; i < d; i++) {
        std::complex<double> w = omega_inv(i);
        std::complex<double> val(1.0, 0.0);
        for (int j = 0; j < d; j++) {
            v_omega_inv(i, j) = val;
            val *= w;
        }
    }

    Eigen::MatrixXcd v_inv_omega_theo = v_omega_inv / d;
    Eigen::VectorXd inv_binom = inv_binom_coef(d - 1);
    Eigen::MatrixXcd res = (v_inv_omega_theo * inv_binom) / d;

    return res;
}
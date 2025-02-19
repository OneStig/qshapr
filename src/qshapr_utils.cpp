#include "qshapr_utils.h"

Eigen::VectorXd inv_binom_coef(int d) {
    Eigen::VectorXd coef(d + 1);
    coef[0] = 1;

    for (int i = 1; i < d / 2 + 1; i++) {
        coef[i] = (coef[i - 1] * (d - i + 1)) / i;
    }

    for (int i = d / 2 + 1; i < d + 1; i++) {
        coef[i] = coef[d - i];
    }

    return coef.cwiseInverse();
}
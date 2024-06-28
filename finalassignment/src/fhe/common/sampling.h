/**
 * @file sampling.h
 * @brief Generate random polynomials.
 */

#include "rns.h"

namespace hehub {

/**
 * @brief Get a random ternary RnsPolynomial object with specified dimensions,
 * which represents a polynomial with coefficients uniform from {-1, 0, 1}. The
 * result will be in NTT form.
 * @param params Polynomial dimensions for initializing the RnsPolynomial.
 * @return RnsPolynomial
 */
RnsPolynomial get_rand_ternary_poly(const RnsPolyParams &params);

/**
 * @brief Get a random uniform RnsPolynomial object with specified dimensions,
 * where the k-th component polynomial's coefficients will be uniform from {0,
 * ..., moduli[k] - 1}.
 * @param params Polynomial dimensions for initializing the RnsPolynomial.
 * @param form Required representation form of the resulting polynomial.
 * @return RnsPolynomial
 */
RnsPolynomial get_rand_uniform_poly(const RnsPolyParams &params,
                                    PolyRepForm form = PolyRepForm::coeff);

/**
 * @brief Get a random gaussian RnsPolynomial object with specified dimensions,
 * which represents a polynomial with coefficients sampled (and rounded) from
 * a gaussian distribution. The result will be in NTT form.
 * @param params Polynomial dimensions for initializing the RnsPolynomial.
 * @param std_dev Standard deviation of the gaussian distribution.
 * @return RnsPolynomial
 */
RnsPolynomial get_rand_gaussian_poly(const RnsPolyParams &params,
                                     double std_dev = 3.2);

/**
 * @brief TODO
 *
 * @param params
 * @param form
 * @return RnsPolynomial
 */
RnsPolynomial get_zero_poly(const RnsPolyParams &params,
                            PolyRepForm form = PolyRepForm::value);

} // namespace hehub

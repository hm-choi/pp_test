import math

import heaan as hn

from engine.HEengine import HEengine
from engine.HEdata import Ciphertext

from coeffs.invSqrt_coeffs import _INV_SQRT_COEFFS
from coeffs.sign_coeffs import _SIGN_PPTEST_DATA


class HEApprox:

    def __init__(
        self,
        engine: HEengine,
    ):

        self.engine = engine
        self._bootstrap_count = 0

    def bootstrap_count(
        self,
    ):
        return self._bootstrap_count

    def _bootstrap(
        self,
        ctxt: Ciphertext,
    ):

        # Bootstrap is performed in place.
        self.engine.bootstrap(ctxt)
        self._bootstrap_count += 1

    def invSqrt(
        self,
        x_half: Ciphertext,
        x_cheb: Ciphertext,
        log_degree=6,
        iteration=7,
        output_scale=1.0,
    ):

        # Use a degree of the form 2^k - 1.
        degree = 2 ** log_degree - 1

        if degree not in _INV_SQRT_COEFFS:
            raise ValueError(f"Unsupported degree: {degree}")

        # Convert the stored coefficient list into the HEaaN
        # Chebyshev coefficient representation.
        coeffs = self.engine._make_cheb_coeffs(_INV_SQRT_COEFFS[degree])

        # Obtain the initial inverse square root approximation
        # from the Chebyshev-mapped ciphertext.
        y = self.cheb_invSqrt(x_cheb, coeffs, output_scale)

        # Refine the initial approximation using Newton iterations.
        return self.he_newton(x_half, y, iteration)

    def he_newton(
        self,
        x_half: Ciphertext,
        y: Ciphertext,
        iteration=10,
    ):

        # Newton iteration for inverse square root:
        #
        # y <- 1.5y - x_half * y^3,
        #
        # where x_half = x / 2.
        if x_half.level() <= 4:
            self._bootstrap(x_half)

        for _ in range(iteration):
            
            if y.level() <= 4:
                self._bootstrap(y)

            tmp_a = self.engine.mult(y, 1.5)

            tmp_b = self.engine.mult(x_half, y)
            y_squared = self.engine.mult(y, y)
            tmp_b = self.engine.mult(tmp_b, y_squared)

            y = self.engine.sub(tmp_a, tmp_b)

        return y

    def cheb_invSqrt(
        self,
        ctxt: Ciphertext,
        coeffs: hn.math.approx.ChebyshevCoefficients,
        output_scale=1.0,
    ):

        x = Ciphertext(ctxt)

        # Ensure that the ciphertext has enough remaining levels
        # before evaluating the Chebyshev polynomial.
        degree = len(coeffs.coeffs) - 1
        self.engine._ensure_cheb_level(x, degree)

        # Evaluate the inverse square root approximation
        # on the Chebyshev-domain input.
        ret = self.engine.evaluate_chebyshev(x, coeffs)

        # Apply an optional public output scaling factor.
        if output_scale != 1.0:
            ret = self.engine.mult(ret, output_scale)

        return ret

    def sign(
        self,
        ctxt: Ciphertext,
    ):

        # Convert each coefficient set used in the composite sign
        # approximation into the HEaaN Chebyshev representation.
        coeffs_list = [
            self.engine._make_cheb_coeffs(coeffs)
            for coeffs in _SIGN_PPTEST_DATA
        ]

        ret = Ciphertext(ctxt)

        for i, coeffs in enumerate(coeffs_list):
            # Before each intermediate polynomial evaluation, ensure
            # that enough ciphertext levels remain for the evaluation.
            #
            # The final polynomial does not need an additional level
            # check because no subsequent polynomial evaluation follows.
            if i != len(coeffs_list) - 1:
                degree = len(coeffs.coeffs) - 1
                self.engine._ensure_cheb_level(ret, degree)

            ret = self.engine.evaluate_chebyshev(ret, coeffs)

        return ret
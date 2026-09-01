import numpy as np

from engine.engine import HEEngine
from hedata.data import HEData
from operators.operator import HEOperator


class SecureBioinformaticsBaseline:
    """
    Secure Bioinformatics-style baseline.

    Paper-style idea:
        - Heavy linear computation is performed homomorphically.
        - Nonlinear division / lookup is handled by MPC in the original paper.

    This implementation:
        - Implements the HE linear part using HEaaN.
        - Provides an HE-only HDiv-normalized adaptation for comparison.
        - Computes Welch T^2 and df, which correspond to the core secure t-test outputs.
    """

    def __init__(self, engine: HEEngine):
        self._engine = engine
        self.__ho = HEOperator(engine)

    # ============================================================
    # 1. Normalized HDiv
    # ============================================================

    def HE_HDiv_Normalized(
        self,
        a: HEData,
        b: HEData,
        b_l: float,
        b_u: float,
        q_max: float,
        num_iter: int = 6,
        bootstrap_level: int = 11,
        bootstrap_threshold: int = 4,
        return_denormalized: bool = True,
    ):
        """
        Normalized Goldschmidt-style homomorphic division.

        Computes:
            q = a / b

        Internally computes:
            q_norm = q / q_max

        so that bootstrapping is applied to small-range ciphertexts.

        b_norm = b / b_u
        a_norm = a / (b_u * q_max)

        a_norm / b_norm
            = [a / (b_u*q_max)] / [b/b_u]
            = a / (b*q_max)
            = q_norm
        """

        ho = self.__ho

        if b_l <= 0 or b_u <= 0:
            raise ValueError("HDiv requires positive denominator bounds.")

        if b_l >= b_u:
            raise ValueError("HDiv requires b_l < b_u.")

        if q_max <= 0:
            raise ValueError("q_max must be positive.")

        b_norm = ho.mult_const(b, 1.0 / b_u)
        a_norm = ho.mult_const(a, 1.0 / (b_u * q_max))

        b_l_norm = b_l / b_u
        b_u_norm = 1.0

        denom = (
            b_l_norm ** 2
            + 6.0 * b_l_norm * b_u_norm
            + b_u_norm ** 2
        )

        # w0 approximation for 1/b_norm
        w = ho.mult_const(b_norm, -8.0 / denom)
        w = ho.add_const(w, 8.0 * (b_l_norm + b_u_norm) / denom)

        a_cur = a_norm
        b_cur = b_norm

        for i in range(num_iter):
            print(
                f"[SecureBio-HDivNorm] iter={i}, "
                f"a_level={a_cur.level()}, "
                f"b_level={b_cur.level()}, "
                f"w_level={w.level()}"
            )

            a_cur = ho.mult(a_cur, w)
            b_cur = ho.mult(b_cur, w)

            # w = 2 - b_cur
            w = ho.mult_const(b_cur, -1.0)
            w = ho.add_const(w, 2.0)

            if a_cur.level() < bootstrap_threshold:
                print("[SecureBio-HDivNorm] bootstrapping a_cur")
                a_cur = ho.do_bootstrapping(a_cur, bootstrap_level)

            if b_cur.level() < bootstrap_threshold:
                print("[SecureBio-HDivNorm] bootstrapping b_cur")
                b_cur = ho.do_bootstrapping(b_cur, bootstrap_level)

            if w.level() < bootstrap_threshold:
                print("[SecureBio-HDivNorm] bootstrapping w")
                w = ho.do_bootstrapping(w, bootstrap_level)

        if return_denormalized:
            return ho.mult_const(a_cur, q_max)

        return a_cur

    # ============================================================
    # 2. Secure Bioinformatics-style Welch t-test
    # ============================================================

    def HE_Welch_T2_DF(
        self,
        x1: HEData,
        x2: HEData,
        V_bounds,
        D_bounds,
        T2_max: float,
        DF_max: float,
        num_iter: int = 6,
        return_debug: bool = False,
    ):
        """
        Secure Bioinformatics-style Welch t-test baseline.

        Original paper-style output:
            t-statistic and degrees of freedom,
            followed by p-value lookup.

        This HE-only adapted version computes:
            T2 = T^2 = (mean1 - mean2)^2 / V
            DF = V^2 / D

        where:
            A = s1^2 / n1
            B = s2^2 / n2
            V = A + B
            D = A^2/(n1-1) + B^2/(n2-1)
        """

        ho = self.__ho

        n1 = x1.size()
        n2 = x2.size()

        print("n1, n2:", n1, n2)

        if n1 <= 1 or n2 <= 1:
            raise ValueError("Welch test requires n1 > 1 and n2 > 1.")

        V_l, V_u = V_bounds
        D_l, D_u = D_bounds

        if V_l <= 0 or D_l <= 0:
            raise ValueError("V_l and D_l must be positive.")

        # ------------------------------------------------------------
        # 1. HE linear statistics
        # ------------------------------------------------------------

        S1 = ho.sum(x1)
        S2 = ho.sum(x2)

        x1_sq = ho.mult(x1, x1)
        x2_sq = ho.mult(x2, x2)

        Q1 = ho.sum(x1_sq)
        Q2 = ho.sum(x2_sq)

        # ------------------------------------------------------------
        # 2. Means
        # ------------------------------------------------------------

        mean1 = ho.mult_const(S1, 1.0 / n1)
        mean2 = ho.mult_const(S2, 1.0 / n2)

        mean_diff = ho.sub(mean1, mean2)
        mean_diff_sq = ho.mult(mean_diff, mean_diff)

        # ------------------------------------------------------------
        # 3. Sample variances
        #
        # s^2 = (Q - S^2/n) / (n-1)
        # ------------------------------------------------------------

        S1_sq = ho.mult(S1, S1)
        S2_sq = ho.mult(S2, S2)

        var1_num = ho.sub(Q1, ho.mult_const(S1_sq, 1.0 / n1))
        var2_num = ho.sub(Q2, ho.mult_const(S2_sq, 1.0 / n2))

        var1 = ho.mult_const(var1_num, 1.0 / (n1 - 1))
        var2 = ho.mult_const(var2_num, 1.0 / (n2 - 1))

        # ------------------------------------------------------------
        # 4. Welch V and D
        # ------------------------------------------------------------

        A = ho.mult_const(var1, 1.0 / n1)
        B = ho.mult_const(var2, 1.0 / n2)

        V = ho.add(A, B)

        A_sq = ho.mult(A, A)
        B_sq = ho.mult(B, B)

        D1 = ho.mult_const(A_sq, 1.0 / (n1 - 1))
        D2 = ho.mult_const(B_sq, 1.0 / (n2 - 1))

        D = ho.add(D1, D2)

        V_sq = ho.mult(V, V)

        # ------------------------------------------------------------
        # 5. Nonlinear part: HDiv adaptation
        # ------------------------------------------------------------

        print("[SecureBio] Compute T2 = mean_diff_sq / V")
        T2 = self.HE_HDiv_Normalized(
            mean_diff_sq,
            V,
            b_l=V_l,
            b_u=V_u,
            q_max=T2_max,
            num_iter=num_iter,
            return_denormalized=True,
        )

        print("[SecureBio] Compute DF = V_sq / D")
        DF = self.HE_HDiv_Normalized(
            V_sq,
            D,
            b_l=D_l,
            b_u=D_u,
            q_max=DF_max,
            num_iter=num_iter,
            return_denormalized=True,
        )

        if return_debug:
            return {
                "T2": T2,
                "DF": DF,
                "mean1": mean1,
                "mean2": mean2,
                "mean_diff": mean_diff,
                "var1": var1,
                "var2": var2,
                "V": V,
                "D": D,
            }

        return T2, DF
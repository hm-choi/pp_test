from engine.engine import HEEngine
import heaan as hn
import numpy as np
from hedata.data import HEData
from operators.inv_sqrt import HEStatistics
from operators.operator import HEOperator
import math, json

class HEHypothesisTestingBaseline:
    def __init__(self, engine: HEEngine):
        self._engine = engine 
        self.__ho = HEOperator(engine)
        self.__hs = HEStatistics(engine)

    def HE_HDiv(
        self,
        a: HEData,
        b: HEData,
        b_l: float,
        b_u: float,
        num_iter: int = 6,
        bootstrap_level: int = 11,
        bootstrap_threshold: int = 4,
    ):
        """
        Homomorphic division baseline using Goldschmidt iteration.

        Approximate:
            a / b

        Inputs:
            a   : encrypted numerator
            b   : encrypted denominator
            b_l : public lower bound of b
            b_u : public upper bound of b

        Notes:
            - For debugging, b_l and b_u can be chosen from plaintext reference.
            - In a real protocol, b_l and b_u must be public bounds.
        """

        ho = self.__ho

        if b_l <= 0 or b_u <= 0:
            raise ValueError("HDiv requires 0 < b_l <= b <= b_u.")

        if b_l >= b_u:
            raise ValueError("HDiv requires b_l < b_u.")

        # ------------------------------------------------------------
        # Initial approximation
        #
        # w0 = -8b / (b_l^2 + 6b_l b_u + b_u^2)
        #      + 8(b_l + b_u) / (b_l^2 + 6b_l b_u + b_u^2)
        # ------------------------------------------------------------

        denom = b_l**2 + 6.0 * b_l * b_u + b_u**2

        w = ho.mult_const(b, -8.0 / denom)
        w = ho.add_const(w, 8.0 * (b_l + b_u) / denom)

        a_cur = a
        b_cur = b

        # ------------------------------------------------------------
        # Goldschmidt iteration
        #
        # a_{i+1} = a_i * w_i
        # b_{i+1} = b_i * w_i
        # w_{i+1} = 2 - b_{i+1}
        # ------------------------------------------------------------

        for i in range(num_iter):
            print(f"[HDiv] iter={i}, a_level={a_cur.level()}, b_level={b_cur.level()}, w_level={w.level()}")

            a_cur = ho.mult(a_cur, w)
            b_cur = ho.mult(b_cur, w)

            # w = 2 - b_cur
            w = ho.mult_const(b_cur, -1.0)
            w = ho.add_const(w, 2.0)

            if a_cur.level() < bootstrap_threshold:
                print("[HDiv] bootstrapping a_cur")
                a_cur = ho.do_bootstrapping(a_cur, bootstrap_level)

            if b_cur.level() < bootstrap_threshold:
                print("[HDiv] bootstrapping b_cur")
                b_cur = ho.do_bootstrapping(b_cur, bootstrap_level)

            if w.level() < bootstrap_threshold:
                print("[HDiv] bootstrapping w")
                w = ho.do_bootstrapping(w, bootstrap_level)

        return a_cur

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
        Normalized Goldschmidt HDiv.

        Computes:
            q = a / b

        But internally computes:
            q_norm = q / q_max

        so that bootstrapping is applied to small-range ciphertexts.

        Denominator is also normalized:
            b_norm = b / b_u in [b_l/b_u, 1]
            a_norm = a / (b_u * q_max)

        Then:
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

        # ------------------------------------------------------------
        # Normalize denominator and numerator
        # ------------------------------------------------------------

        b_norm = ho.mult_const(b, 1.0 / b_u)
        a_norm = ho.mult_const(a, 1.0 / (b_u * q_max))

        b_l_norm = b_l / b_u
        b_u_norm = 1.0

        # ------------------------------------------------------------
        # Initial approximation for 1 / b_norm
        # ------------------------------------------------------------

        denom = (
            b_l_norm ** 2
            + 6.0 * b_l_norm * b_u_norm
            + b_u_norm ** 2
        )

        w = ho.mult_const(b_norm, -8.0 / denom)
        w = ho.add_const(w, 8.0 * (b_l_norm + b_u_norm) / denom)

        a_cur = a_norm
        b_cur = b_norm

        # ------------------------------------------------------------
        # Goldschmidt iteration
        # ------------------------------------------------------------

        for i in range(num_iter):
            print(
                f"[HDivNorm] iter={i}, "
                f"a_level={a_cur.level()}, "
                f"b_level={b_cur.level()}, "
                f"w_level={w.level()}"
            )

            a_cur = ho.mult(a_cur, w)
            b_cur = ho.mult(b_cur, w)

            # w = 2 - b_cur
            w = ho.mult_const(b_cur, -1.0)
            w = ho.add_const(w, 2.0)

            # 이제 a_cur는 q_norm이므로 값이 0~1 근처여야 함
            if a_cur.level() < bootstrap_threshold:
                print("[HDivNorm] bootstrapping a_cur")
                a_cur = ho.do_bootstrapping(a_cur, bootstrap_level)

            if b_cur.level() < bootstrap_threshold:
                print("[HDivNorm] bootstrapping b_cur")
                b_cur = ho.do_bootstrapping(b_cur, bootstrap_level)

            if w.level() < bootstrap_threshold:
                print("[HDivNorm] bootstrapping w")
                w = ho.do_bootstrapping(w, bootstrap_level)

        if return_denormalized:
            return ho.mult_const(a_cur, q_max)

        return a_cur
    def HE_Welch_T2_DF_HDiv_Baseline(
        self,
        x1: HEData,
        x2: HEData,
        V_bounds,
        D_bounds,
        T2_max: float = None,
        DF_max: float = None,
        R: float = None,
        num_iter: int = 6,
        return_debug: bool = False,
    ):
        """
        Normalized-HDiv-based Welch baseline.

        Computes:
            T2 = T^2 = (mean1 - mean2)^2 / V
            DF = V^2 / D

        where:
            V = s1^2/n1 + s2^2/n2
            D = (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1)

        Important:
            This uses HE_HDiv_Normalized, not HE_HDiv.
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

        # ============================================================
        # 1. Sufficient statistics
        # ============================================================

        S1 = ho.sum(x1)
        S2 = ho.sum(x2)

        x1_sq = ho.mult(x1, x1)
        x2_sq = ho.mult(x2, x2)

        Q1 = ho.sum(x1_sq)
        Q2 = ho.sum(x2_sq)

        # ============================================================
        # 2. Mean and mean_diff
        # ============================================================

        mean1 = ho.mult_const(S1, 1.0 / n1)
        mean2 = ho.mult_const(S2, 1.0 / n2)

        mean_diff = ho.sub(mean1, mean2)
        mean_diff_sq = ho.mult(mean_diff, mean_diff)

        # ============================================================
        # 3. Sample variance
        #
        # s² = (Q - S²/n) / (n-1)
        # ============================================================

        S1_sq = ho.mult(S1, S1)
        S2_sq = ho.mult(S2, S2)

        var1_num = ho.sub(Q1, ho.mult_const(S1_sq, 1.0 / n1))
        var2_num = ho.sub(Q2, ho.mult_const(S2_sq, 1.0 / n2))

        var1 = ho.mult_const(var1_num, 1.0 / (n1 - 1))
        var2 = ho.mult_const(var2_num, 1.0 / (n2 - 1))

        # ============================================================
        # 4. V and D
        # ============================================================

        A = ho.mult_const(var1, 1.0 / n1)
        B = ho.mult_const(var2, 1.0 / n2)

        V = ho.add(A, B)

        A_sq = ho.mult(A, A)
        B_sq = ho.mult(B, B)

        D1 = ho.mult_const(A_sq, 1.0 / (n1 - 1))
        D2 = ho.mult_const(B_sq, 1.0 / (n2 - 1))

        D = ho.add(D1, D2)

        V_sq = ho.mult(V, V)

        # ============================================================
        # 5. Quotient upper bounds
        # ============================================================

        if T2_max is None:
            if R is None:
                raise ValueError("Either T2_max or R must be provided.")
            # conservative public bound: mean_diff^2 <= R^2
            T2_max = (R ** 2) / V_l

        if DF_max is None:
            # public bound: DF = V^2/D <= V_u^2/D_l
            DF_max = (V_u ** 2) / D_l

        print("[HDivNorm baseline] T2_max:", T2_max)
        print("[HDivNorm baseline] DF_max:", DF_max)

        # ============================================================
        # 6. Normalized HDiv
        # ============================================================

        print("[HDivNorm baseline] Compute T2 = mean_diff_sq / V")
        T2 = self.HE_HDiv_Normalized(
            mean_diff_sq,
            V,
            b_l=V_l,
            b_u=V_u,
            q_max=T2_max,
            num_iter=num_iter,
            return_denormalized=True,
        )

        print("[HDivNorm baseline] Compute DF = V_sq / D")
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
                "T2_max": T2_max,
                "DF_max": DF_max,
            }

        return T2, DF
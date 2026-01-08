import numpy as np
from scipy.stats import shapiro, ttest_1samp, wilcoxon, rankdata

def _rank_biserial_from_diffs(diffs, zero_method="wilcox"):
    """
    Rank-biserial correlation (RBC) for Wilcoxon signed-rank test, computed from diffs.
    RBC = (W_plus - W_minus) / (W_plus + W_minus) = (2*W_plus / total_rank_sum) - 1

    Notes:
    - We rank |diffs| with average ranks for ties.
    - Handling zeros:
        * "wilcox": drop zeros (matches common Wilcoxon default behavior)
        * "pratt": include zeros in ranking but give them zero sign contribution (common interpretation)
        * "zsplit": split zero ranks half to + and half to - (rare; included)
    """
    diffs = np.asarray(diffs, dtype=float)

    is_zero = diffs == 0
    if zero_method == "wilcox":
        diffs_eff = diffs[~is_zero]
        if diffs_eff.size == 0:
            return np.nan
        abs_vals = np.abs(diffs_eff)
        ranks = rankdata(abs_vals, method="average")
        signs = np.sign(diffs_eff)
        w_plus = np.sum(ranks[signs > 0])
        w_minus = np.sum(ranks[signs < 0])
        denom = w_plus + w_minus
        return np.nan if denom == 0 else (w_plus - w_minus) / denom

    # For pratt and zsplit we rank all abs(diffs), including zeros
    abs_vals = np.abs(diffs)
    ranks = rankdata(abs_vals, method="average")

    signs = np.sign(diffs)
    w_plus = np.sum(ranks[signs > 0])
    w_minus = np.sum(ranks[signs < 0])

    if zero_method == "pratt":
        # zeros contribute nothing (neither + nor -)
        pass
    elif zero_method == "zsplit":
        # split the zero ranks evenly between + and -
        w_zero = np.sum(ranks[is_zero])
        w_plus += 0.5 * w_zero
        w_minus += 0.5 * w_zero
    else:
        raise ValueError("zero_method must be one of: 'wilcox', 'pratt', 'zsplit'.")

    denom = w_plus + w_minus
    return np.nan if denom == 0 else (w_plus - w_minus) / denom


def paired_test_auto(
    diffs,
    *,
    alpha_normality=0.05,
    alternative="two-sided",
    zero_method="wilcox",
    nan_policy="omit",
):
    """
    Automatic paired test from a 1D array of paired differences (x - y):
      1) Shapiro–Wilk normality test on diffs
      2) If normal at alpha_normality -> paired t-test (one-sample t-test on diffs)
         else -> Wilcoxon signed-rank test
      3) Report effect size:
         - t-test: Cohen's d (paired) = mean(diffs) / sd(diffs, ddof=1)
         - Wilcoxon: rank-biserial correlation (RBC)

    Returns a results dict (easy to log or serialize).
    """
    diffs = np.asarray(diffs, dtype=float)

    if nan_policy == "omit":
        diffs = diffs[~np.isnan(diffs)]
    elif nan_policy != "raise" and np.isnan(diffs).any():
        raise ValueError("nan_policy must be 'omit' or 'raise'.")

    n = diffs.size
    if n < 3:
        raise ValueError("Need at least 3 paired differences (Shapiro–Wilk requires >=3).")

    # Shapiro–Wilk
    sh_w, sh_p = shapiro(diffs)
    is_normal = sh_p >= alpha_normality

    results = {
        "n": int(n),
        "diff_mean": float(np.mean(diffs)),
        "diff_sd": float(np.std(diffs, ddof=1)) if n > 1 else np.nan,
        "diff_median": float(np.median(diffs)),
        "diff_iqr": float(np.subtract(*np.percentile(diffs, [75, 25]))),
        "shapiro_W": float(sh_w),
        "shapiro_p": float(sh_p),
        "normal_at_alpha": bool(is_normal),
        "alpha_normality": float(alpha_normality),
        "chosen_test": None,
        "statistic": None,
        "p_value": None,
        "effect_size_name": None,
        "effect_size": None,
        "alternative": alternative,
        "zero_method": zero_method,
    }

    if is_normal:
        # Paired t-test == one-sample t-test on differences
        t_stat, p = ttest_1samp(diffs, popmean=0.0, alternative=alternative)

        sd = np.std(diffs, ddof=1)
        d = np.nan if sd == 0 else np.mean(diffs) / sd  # Cohen's d (paired)

        results.update({
            "chosen_test": "paired_t_test (ttest_1samp on diffs)",
            "statistic": float(t_stat),
            "p_value": float(p),
            "effect_size_name": "Cohen_d_paired",
            "effect_size": float(d) if np.isfinite(d) else np.nan,
        })
    else:
        w_stat, p = wilcoxon(
            diffs,
            alternative=alternative,
            zero_method=zero_method,
        )
        rbc = _rank_biserial_from_diffs(diffs, zero_method=zero_method)

        results.update({
            "chosen_test": "wilcoxon_signed_rank",
            "statistic": float(w_stat),
            "p_value": float(p),
            "effect_size_name": "rank_biserial_correlation",
            "effect_size": float(rbc) if np.isfinite(rbc) else np.nan,
        })

    return results


def format_report(res, *, digits=4):
    """Human-readable, single-block report."""
    n = res["n"]
    alt = res["alternative"]
    lines = []
    lines.append(f"n = {n}")
    lines.append(
        f"Differences: mean = {res['diff_mean']:.{digits}g}, "
        f"SD = {res['diff_sd']:.{digits}g}, "
        f"median = {res['diff_median']:.{digits}g}, "
        f"IQR = {res['diff_iqr']:.{digits}g}"
    )
    lines.append(
        f"Shapiro–Wilk normality: W = {res['shapiro_W']:.{digits}g}, "
        f"p = {res['shapiro_p']:.{digits}g} "
        f"(α = {res['alpha_normality']}) -> "
        f"{'normal' if res['normal_at_alpha'] else 'non-normal'}"
    )
    lines.append(f"Chosen test: {res['chosen_test']} (alternative = '{alt}')")
    lines.append(
        f"Test result: statistic = {res['statistic']:.{digits}g}, "
        f"p = {res['p_value']:.{digits}g}"
    )
    lines.append(
        f"Effect size ({res['effect_size_name']}): {res['effect_size']:.{digits}g}"
    )
    return "\n".join(lines)


if __name__ == "__main__":
    # tier 1 - knee compression then shear
    # print('\nKnee compression:'); diffs = np.array([-0.102287073, 1.781602499, 0.450474826, 0.717245925, 1.292819939, 1.890431661, 0.681113247])  # example (x - y)
    # print('\nKnee shear:'); diffs = np.array([-0.153084042, 0.203847533, -0.132418695, 0.01031966, -0.323410273, 0.267793971, -0.014709844])  # example (x - y)
    # # tier 2 - quads
    # print('\nQuads compression:'); diffs = np.array([-0.150932521, 1.416347522, 0.578245078, 0.615183912, 1.835710767, 1.113863595, 0.372006246])  # example (x - y)
    # # tier 3 - exploration (hip compression/shear and then ankle compression/shear)
    # print('\nHip compression:'); diffs = np.array([-0.589784077, 1.270201452, 0.102214341, 0.53888374, 0.714631287, 0.276242841, -1.02234926])  # example (x - y)
    # print('\nHip shear:'); diffs = np.array([-0.089662651, 0.300522123, 0.183639582, 0.344884016, 0.474511264, 0.088238692, -0.536236139])  # example (x - y)
    # print('\nAnkle compression:'); diffs = np.array([0.15269778, 1.854828272, -0.00101458, 0.616307171, -0.526644559, 1.644364595, 1.132647017])  # example (x - y)
    print('\nAnkle shear:'); diffs = np.array([0.162786588, 0.931908201, 0.078858327, 0.104042828, -0.125004089, 1.360212357, 0.657593848])  # example (x - y)

    res = paired_test_auto(
        diffs,
        alpha_normality=0.05,
        alternative="two-sided",
        zero_method="wilcox",   # or "pratt" if you want to keep zeros
        nan_policy="omit",
    )

    print(format_report(res))


# current results N=7 Tiered approach, and ordered. (note bonferroni correction not applied here)
"""
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel> python .\stats_both.py

Knee compression:
n = 7
Differences: mean = 0.9588, SD = 0.7276, median = 0.7172, IQR = 0.9714
Shapiro–Wilk normality: W = 0.9438, p = 0.6727 (α = 0.05) -> normal
Chosen test: paired_t_test (ttest_1samp on diffs) (alternative = 'two-sided')
Test result: statistic = 3.486, p = 0.01304
Effect size (Cohen_d_paired): 1.318
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel> python .\stats_both.py

Knee shear:
n = 7
Differences: mean = -0.02024, SD = 0.2067, median = -0.01471, IQR = 0.2498
Shapiro–Wilk normality: W = 0.9626, p = 0.8405 (α = 0.05) -> normal
Chosen test: paired_t_test (ttest_1samp on diffs) (alternative = 'two-sided')
Test result: statistic = -0.2591, p = 0.8042
Effect size (Cohen_d_paired): -0.09792
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel> python .\stats_both.py

Quads compression:
n = 7
Differences: mean = 0.8258, SD = 0.673, median = 0.6152, IQR = 0.79
Shapiro–Wilk normality: W = 0.9777, p = 0.9478 (α = 0.05) -> normal
Chosen test: paired_t_test (ttest_1samp on diffs) (alternative = 'two-sided')
Test result: statistic = 3.247, p = 0.01754
Effect size (Cohen_d_paired): 1.227
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel> python .\stats_both.py

Hip compression:
n = 7
Differences: mean = 0.1843, SD = 0.7807, median = 0.2762, IQR = 0.8705
Shapiro–Wilk normality: W = 0.9721, p = 0.9135 (α = 0.05) -> normal
Chosen test: paired_t_test (ttest_1samp on diffs) (alternative = 'two-sided')
Test result: statistic = 0.6246, p = 0.5552
Effect size (Cohen_d_paired): 0.2361
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel> python .\stats_both.py

Hip shear:
n = 7
Differences: mean = 0.1094, SD = 0.3386, median = 0.1836, IQR = 0.3234
Shapiro–Wilk normality: W = 0.9108, p = 0.4015 (α = 0.05) -> normal
Chosen test: paired_t_test (ttest_1samp on diffs) (alternative = 'two-sided')
Test result: statistic = 0.8551, p = 0.4253
Effect size (Cohen_d_paired): 0.3232
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel> python .\stats_both.py

Ankle compression:
n = 7
Differences: mean = 0.6962, SD = 0.8865, median = 0.6163, IQR = 1.313
Shapiro–Wilk normality: W = 0.9539, p = 0.765 (α = 0.05) -> normal
Chosen test: paired_t_test (ttest_1samp on diffs) (alternative = 'two-sided')
Test result: statistic = 2.078, p = 0.08299
Effect size (Cohen_d_paired): 0.7853
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel> python .\stats_both.py

Ankle shear:
n = 7
Differences: mean = 0.4529, SD = 0.5438, median = 0.1628, IQR = 0.7033
Shapiro–Wilk normality: W = 0.8977, p = 0.3173 (α = 0.05) -> normal
Chosen test: paired_t_test (ttest_1samp on diffs) (alternative = 'two-sided')
Test result: statistic = 2.204, p = 0.06976
Effect size (Cohen_d_paired): 0.8329
(moco451) PS C:\Users\jonstingel\code\musclemodel\muscleEnergyModel>
"""
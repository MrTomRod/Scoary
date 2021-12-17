from scipy.stats import boschloo_exact


def simple_boschloo(test_string: str) -> (float, float):
    """
    can be pickled
    """
    c1r1, c2r1, c1r2, c2r2 = (int(i) for i in test_string.split(','))
    boschloo_res = boschloo_exact([[c1r1, c2r1], [c1r2, c2r2]])
    return boschloo_res.pvalue, boschloo_res.statistic

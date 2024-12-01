import numpy as np
from kneed import KneeLocator
import matplotlib.pyplot as plt
import logging
from numba import jit

log = logging.getLogger("my_logger")

plt_log = logging.getLogger("matplotlib")
plt_log.setLevel(logging.ERROR)

logging.getLogger("PIL").setLevel(logging.WARNING)


# Numba-optimized functions
@jit(nopython=True, cache=True)
def entropy(counts):
    counts_sum = np.sum(counts)
    frequencies = counts / counts_sum
    H = 0.0
    for p in frequencies:
        if p != 0:
            H = H - p * np.log(p)
    return H


@jit(nopython=True, cache=True)
def get_even_distribution(n_obs, n_bins):
    quotient = n_obs // n_bins
    remainder = n_obs % n_bins
    values = np.empty(n_bins, dtype=np.int64)
    values[: n_bins - remainder] = quotient
    values[n_bins - remainder :] = quotient + 1
    return values


@jit(nopython=True, cache=True)
def norm_entropy(counts):
    counts_sum = np.sum(counts)
    n_bins = len(counts)

    if counts_sum == 1:
        return 1.0

    values = get_even_distribution(counts_sum, n_bins)
    max_possible_entropy = entropy(values)
    return entropy(counts) / max_possible_entropy


@jit(nopython=True, cache=True)
def cumsum(arr):
    result = np.empty(len(arr) + 1)
    result[0] = 0
    sum_val = 0
    for i in range(len(arr)):
        sum_val += arr[i]
        result[i + 1] = sum_val
    return result


@jit(nopython=True, cache=True)
def trapezoid(y, dx):
    n = len(y) - 1
    area = 0.0
    for i in range(n):
        area += (y[i] + y[i + 1]) * dx / 2
    return area


@jit(nopython=True, cache=True)
def gini_coeff(values):
    total = np.sum(values)
    if total == 0:
        return 0.0

    norm_values = values / total
    sorted_values = np.sort(norm_values)

    # Calculate cumulative distribution including 0
    cum_distr = cumsum(sorted_values)

    # Compute area under Lorenz curve using trapezoidal rule
    n_classes = len(values)
    dx = 1.0 / n_classes
    under_lorenz = trapezoid(cum_distr, dx)

    return (0.5 - under_lorenz) / 0.5


@jit(nopython=True, cache=True)
def norm_gini_coeff(values):
    n_bins = len(values)
    n_obs = np.sum(values)

    if n_obs == 0:
        return 0.0

    gini = gini_coeff(values)

    min_values = get_even_distribution(n_obs, n_bins)
    min_gini = gini_coeff(min_values)

    max_values = np.zeros(n_bins, dtype=np.int64)
    max_values[-1] = n_obs
    max_gini = gini_coeff(max_values)

    if max_gini - min_gini == 0:
        return 0.0
    return (gini - min_gini) / (max_gini - min_gini)


def find_knee(df, out_plot_name):
    """
    Find the knee point of a curve.
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with two columns: "x" and "y".
    Returns
    -------
    kneedle : float
        Knee point of the curve.
    """
    logging.info("Finding knee point for the gini and entropy values...")
    df = df[["norm_gini", "norm_entropy"]].copy()
    df["norm_gini"] = df["norm_gini"].round(1)
    df2_round = df.groupby("norm_gini", as_index=False).agg({"norm_entropy": "mean"})
    # Find knee point
    x = df2_round["norm_gini"].values
    y = df2_round["norm_entropy"].values
    kneedle = KneeLocator(
        x, y, S=1.0, curve="concave", direction="decreasing", online=True
    )

    if kneedle.knee is None:
        return None, None
    min_gini = round(kneedle.elbow, 3)
    min_entropy = round(kneedle.knee_y, 3)
    logging.info(f"Knee point found: gini={min_gini}, entropy={min_entropy}")
    if out_plot_name is not None:
        logging.info("Knee point figure saved at kneedle.png")
        fig, ax = plt.subplots(nrows=1, ncols=1)
        plt.suptitle("Knee Point")
        plt.title(f"Knee point: gini={min_gini}, entropy={min_entropy}")
        plt.plot(kneedle.x, kneedle.y, "c", label="data")
        plt.vlines(
            kneedle.knee,
            plt.ylim()[0],
            plt.ylim()[1],
            linestyles="--",
            label="knee/elbow",
        )
        plt.xlabel("Normalized Gini coefficient")
        plt.ylabel("Normalized Entropy")
        plt.legend(loc="best")  # create figure & 1 axis
        fig.savefig(out_plot_name, dpi=300)
        plt.close(fig)

    return min_gini, min_entropy

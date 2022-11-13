import numpy as np
from kneed import KneeLocator
import matplotlib.pyplot as plt
import logging

log = logging.getLogger("my_logger")

plt_log = logging.getLogger("matplotlib")
plt_log.setLevel(logging.ERROR)

logging.getLogger("PIL").setLevel(logging.WARNING)

# Code from: https://www.frontiersin.org/articles/10.3389/fmicb.2022.918015/full
def entropy(counts):
    counts_vector = np.array(counts)
    frequencies = counts_vector / counts_vector.sum()
    H = 0
    for p in frequencies:
        if p != 0:
            H = H - p * np.log(p)
    return H


def norm_entropy(counts):
    """

    Entropy divided by the maximum entropy possible with that number of counts
    and that number of bins.

    Parameters
    ----------
    counts : array-like object
        Counts associated to each class.
    Returns
    -------
    rel_possible_ent : float
        Ranges from 0, when entropy is 0, to 1, when entropy is the maximum
        possible entropy. The maximum possible entropy depends on the number of
        counts and bins, and it's achieved when the counts are distributed as
        evenly as possible among the bins. Example: with 10 bins and 12 counts,
        maximum possible entropy is the entropy of the distribution where 2
        bins contain 2 counts, and 8 bins contain 1 count.
    """

    counts_vector = np.array(counts)
    n_obs = counts_vector.sum()
    n_bins = len(counts_vector)
    if n_obs == 1:
        rel_possible_ent = 1
    else:
        # Compute max entropy possible with that number of obs and bins
        quotient = n_obs // n_bins
        remainder = n_obs % n_bins
        chunk_1 = np.repeat(quotient, n_bins - remainder)
        chunk_2 = np.repeat(quotient + 1, remainder)
        values = np.hstack((chunk_1, chunk_2))  # values distr as evenly as possible
        max_possible_entropy = entropy(values)
        # Compute relative entropy
        rel_possible_ent = entropy(counts) / max_possible_entropy
    return rel_possible_ent


def gini_coeff(values_for_each_class):
    """
    Gini coefficient measures distribution inequality.
    Parameters
    ----------
    values_for_each_class : array-like object
        Values associated to each class.
        They don't need to be already sorted and/or normalized.
    Returns
    -------
    gini_coeff : float
        Ranges from 0 (perfect equality) to 1 (maximal inequality).
    """

    import numpy as np

    values = np.array(values_for_each_class)
    norm_values = values / values.sum()  # normalize

    # Generate Lorenz curve
    norm_values.sort()
    cum_distr = np.cumsum(norm_values)
    cum_distr = list(cum_distr)
    cum_distr.insert(0, 0)

    # Get area under Lorenz curve
    n_classes = len(cum_distr) - 1
    under_lorenz = np.trapz(y=cum_distr, dx=1 / n_classes)

    # Area under Perfect Equality curve
    # It's the area of a triangle with base = 1 and height = 1
    under_PE = 0.5

    # Compute Gini coefficient
    gini_coeff = (under_PE - under_lorenz) / under_PE

    return gini_coeff


def norm_gini_coeff(values_for_each_class):
    """
    Normalized Gini coefficient.
    The minimum and maximum possible Gini coefficient with that number of
    bins and observations are computed. Then, norm_Gini_coefficient is
    defined as
    norm_Gini_coefficient := (Gini - min_Gini) / (max_Gini - min_Gini)
    Parameters
    ----------
    values_for_each_class : array-like object
        Values associated to each class.
        They don't need to be already sorted and/or normalized.
    Returns
    -------
    norm_gini_coeff : float
        Ranges from 0 (minimal inequality possible) to 1 (maximal
        inequality possible).
    """

    import numpy as np

    # Compute Gini coefficient
    nuber_of_bins = len(values_for_each_class)
    number_of_obs = np.array(values_for_each_class).sum()
    Gini = gini_coeff(values_for_each_class)

    # Compute minimum possible Gini coefficient
    quotient = number_of_obs // nuber_of_bins
    remainder = number_of_obs % nuber_of_bins
    chunk_1 = np.repeat(quotient, nuber_of_bins - remainder)
    chunk_2 = np.repeat(quotient + 1, remainder)
    vect = np.hstack((chunk_1, chunk_2))  # values distr as evenly as possible
    min_Gini = gini_coeff(vect)

    # Compute maximum possible Gini coefficient
    chunk_1 = np.repeat(0, nuber_of_bins - 1)
    chunk_2 = np.repeat(number_of_obs, 1)
    vect = np.hstack((chunk_1, chunk_2))  # values distr as unevenly as possible
    vect = [int(v) for v in vect]
    max_Gini = gini_coeff(vect)

    # Compute normalized Gini coefficient
    if max_Gini - min_Gini == 0:
        norm_gini = 0
    else:
        norm_gini = (Gini - min_Gini) / (max_Gini - min_Gini)

    return norm_gini


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
    logging.info("Knee point figure saved at kneedle.png")
    fig, ax = plt.subplots(nrows=1, ncols=1)
    plt.suptitle("Knee Point")
    plt.title(f"Knee point: gini={min_gini}, entropy={min_entropy}")
    plt.plot(kneedle.x, kneedle.y, "c", label="data")
    plt.vlines(
        kneedle.knee, plt.ylim()[0], plt.ylim()[1], linestyles="--", label="knee/elbow"
    )
    plt.xlabel("Normalized Gini coefficient")
    plt.ylabel("Normalized Entropy")
    plt.legend(loc="best")  # create figure & 1 axis
    fig.savefig(out_plot_name, dpi=300)
    plt.close(fig)

    return min_gini, min_entropy

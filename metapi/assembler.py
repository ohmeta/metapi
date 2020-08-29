#!/usr/bin/env python

import numpy as np
from metapi import tooler


def assembler_init(contigs_length_range, groups):
    global CONTIGS_LENGTH_RANGES__
    global GROUP_BY_

    temp = sorted(list(set(contigs_length_range)))
    if (temp is None) or (temp == [0]) or (len(temp) == 0):
        CONTIGS_LENGTH_RANGES__ = [(0, 1500), (1500, 2000), (2000, 2500), 2500]
    elif temp[0] == 0:
        CONTIGS_LENGTH_RANGES__ = [(i, j) for i, j in zip(temp[0:-1], temp[1:])] + [
            temp[-1]
        ]
    else:
        CONTIGS_LENGTH_RANGES__ = [(i, j) for i, j in zip([0] + temp[0:-1], temp)] + [
            temp[-1]
        ]

    GROUP_BY_ = groups


def Nx(x):
    Nx_len = np.zeros(21)
    Nx_rate = np.linspace(0, 1, 21)
    Nx_len_ = Nx_rate * sum(x)

    len_sum = 0
    for i in sorted(x)[::-1]:
        len_sum += i
        for j in range(0, len(Nx_len)):
            if (Nx_len[j] == 0.0) and (len_sum >= Nx_len_[j]):
                Nx_len[j] = i
    return list(Nx_len)


def GC_content(x):
    return (x[("#G", "sum")] + x[("#C", "sum")]) / x[("length", "sum")]


def cal_len_range(x):
    len_dict = {i: 0 for i in CONTIGS_LENGTH_RANGES__}
    for seq_len in x:
        if seq_len >= CONTIGS_LENGTH_RANGES__[-1]:
            len_dict[CONTIGS_LENGTH_RANGES__[-1]] += 1
        else:
            for i in CONTIGS_LENGTH_RANGES__[:-1]:
                if i[0] <= seq_len < i[1]:
                    len_dict[i] += 1
                    break
    return tuple([len_dict[i] for i in CONTIGS_LENGTH_RANGES__])


def cumulative_len(x):
    """
    for x, y in enumerate(len_y):
        plt.plot(x, y)
    """
    len_y = [0]
    for l in sorted(x)[::-1]:
        len_y.append(len_y[-1] + l)
    return len_y


def parse_assembly(stats_file):
    df_ = tooler.parse(stats_file)
    if df_ is not None:
        df = (
            df_.groupby(GROUP_BY_)
            .agg(
                {
                    "chr": ["count"],
                    "length": ["sum", "min", "max", "std", Nx, cal_len_range],
                    "#A": ["sum", "min", "max", "std"],
                    "#C": ["sum", "min", "max", "std"],
                    "#G": ["sum", "min", "max", "std"],
                    "#T": ["sum", "min", "max", "std"],
                    "#CpG": ["sum", "min", "max", "std"],
                }
            )
            .reset_index()
        )
        df["GC_content"] = df.apply(lambda x: GC_content(x), axis=1)

        N_ = []
        LEN_LIST = []

        for N in range(0, 105, 5):
            N_.append("N" + str(N))

        for i in range(0, 21):
            df[("length", N_[i])] = df.apply(lambda x: x[("length", "Nx")][i], axis=1)


        for i in range(0, len(CONTIGS_LENGTH_RANGES__) - 1):
            len_tuple = CONTIGS_LENGTH_RANGES__[i]
            LEN_LIST.append(len_tuple[0])
            LEN_LIST.append(len_tuple[1])
            len_range_str = "[%d, %d)" % (len_tuple[0], len_tuple[1])
            df[("length", len_range_str)] = df.apply(
                lambda x: x[("length", "cal_len_range")][i], axis=1
            )

        LEN_LIST.append(CONTIGS_LENGTH_RANGES__[-1])

        df[("length", "[%d, )" % CONTIGS_LENGTH_RANGES__[-1])] = df.apply(
            lambda x: x[("length", "cal_len_range")][-1], axis=1
        )

        LEN_LIST = sorted(list(set(LEN_LIST)))
        for i in range(0, len(LEN_LIST)):
            df[("length", ">=%d" % LEN_LIST[i])] = df.apply(
                lambda x: sum(x[("length", "cal_len_range")][i:]),
                axis=1
            )
        return df
    else:
        return None

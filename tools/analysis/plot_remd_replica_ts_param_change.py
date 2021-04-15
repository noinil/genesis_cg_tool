#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(remd_par_on_rep_fname, fig_all_in_one, rep_indices):
    remd_rep_ts_data = np.loadtxt(remd_par_on_rep_fname)

    num_timesteps = remd_rep_ts_data[-1, 0]
    num_replicas  = np.shape(remd_rep_ts_data)[1] - 1

    remd_timesteps = remd_rep_ts_data[:, 0]
    remd_parameters = remd_rep_ts_data[:, 1:]

    if len(rep_indices) == 0:
        rep_indices = [i for i in range(num_replicas)]

    # ==================
    # axis tick settings
    # ==================
    ax_tk_dx = num_timesteps / 5
    ax_tk_power = int( np.log10(num_timesteps) )
    ax_tk_ticks = [i * ax_tk_dx for i in range(6)]
    ax_tk_tick_lbs = [round( tk / 10**ax_tk_power, 3 ) for tk in ax_tk_ticks]

    # =========
    # Plotting!
    # =========
    if fig_all_in_one:
        fig, ax = plt.subplots(1, 1, figsize=(9, 5), constrained_layout=True)
        for i, j in enumerate( rep_indices ):
            ax.plot(remd_timesteps, remd_parameters[:, j], c=[i / len(rep_indices), 0, 1 - i / len(rep_indices)])
        ax.set_xticks(ax_tk_ticks)
        ax.set_xticklabels(ax_tk_tick_lbs, fontsize=12)
        ax.set_xlim(0, num_timesteps)
        ax.set_xlabel(r"MD steps ($10^{0}$)".format(ax_tk_power), fontsize=16)
        ax.set_ylim(1, num_replicas)

        figname = remd_par_on_rep_fname[:-4] + "_all_in_one.svg"
        plt.savefig(figname)
    else:
        fig, axes = plt.subplots(len(rep_indices), 1, figsize=(9, 2 * len(rep_indices)**0.5), constrained_layout=True, sharex=True, sharey=False)
        for i, j in enumerate( rep_indices ):
            axes[i].plot(remd_timesteps, remd_parameters[:, j], c=[i / len(rep_indices), 0, 1 - i / len(rep_indices)])
            axes[i].set_xticks(ax_tk_ticks)
            axes[i].set_xlim(0, num_timesteps)
            axes[i].set_ylim(1, num_replicas)
        axes[-1].set_xticklabels(ax_tk_tick_lbs, fontsize=12)
        axes[-1].set_xlabel(r"MD steps ($10^{0}$)".format(ax_tk_power), fontsize=16)

        figname = remd_par_on_rep_fname[:-4] + "_share_x.svg"
        plt.savefig(figname)

if __name__ == '__main__':
    import argparse

    def parse_arguments():
        parser = argparse.ArgumentParser(description='Plot replica/parameter walk in param/rep space.')
        parser.add_argument('filename', type=str, help="file name of remd log")
        parser.add_argument('-A', '--all-in-one', help="plot the all-in-one figure or not. (Default: no)", action='store_true')
        parser.add_argument('-n', '--rep-ID', metavar='N', type=int, nargs='+', help='indices of replicas to plot.', default=[])
        return parser.parse_args()

    args = parse_arguments()

    main(args.filename, args.all_in_one, args.rep_ID)

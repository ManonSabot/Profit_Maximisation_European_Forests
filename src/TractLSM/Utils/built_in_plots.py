# -*- coding: utf-8 -*-

"""
functions used to plot preset diagnostic plots (comparisons of diurnal
cycles)

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""
__title__ = "Plotting functions of TractLSM"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (10.02.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import os  # check for paths
import sys  # check for system version
import pandas as pd  # to read/write csv
import numpy as np  # array manipulations, math operators

# plotting modules
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.lines as mlines  # custom legends
from matplotlib import transforms, patches  # handle text boxes
from matplotlib.backends.backend_pdf import PdfPages  # metadata in pdf

try:
    from pdfminer.pdfparser import PDFParser  # extract data from a pdf
    from pdfminer.pdfdocument import PDFDocument  # read metadata in pdf

except (ImportError, ModuleNotFoundError):
    from pdfminer.pdfparser import PDFParser, PDFDocument


# ======================================================================

class FigInfo(object):

    """
    Metadata information stored in the diagnostic figures.

    """

    def __init__(self, df):

        # read params & met data
        p = df.iloc[0]  # fixed params & first doy

        self.doy = str(int(p['doy']) + 1)
        self.doy2 = str(int(df['doy'].max()) + 1)  # last doy
        self.minT = str(df['Tair'].min())
        self.maxT = str(df['Tair'].max())
        self.minVPD = str(df['VPD'].min())
        self.maxVPD = str(df['VPD'].max())
        self.PPT = str(df['precip'].sum())
        self.PPFD = str(df['PPFD'].mean())
        self.u = str(p['u'])
        self.maxleaf = str(p['max_leaf_width'])
        self.Ps = str(p['Ps'])
        self.P50 = str(p['P50'])
        self.P88 = str(p['P88'])
        self.kmax = str(p['kmax'])

        return


def check_fig_exists(md, fig_dir, name_str, force_write):

    """
    Checks whether the exact same metadata is not already present in
    another figure.

    Arguments:
    ----------
    md: class or pandas object
        parameters the model has been run with

    fig_dir: string
        directory (path) where the other comparable figures are stored

    name_str: string
        specific file name

    force_write: string
        if 'yes' then same figure gets written over

    Returns:
    --------
    num: string
        number extension to give the figure in case the same name
        already exists but not with the exact same metadata

    """

    # does output exists for same P50, P88, kmax, with which tag #?
    highest_num = []

    try:

        for file in os.listdir(fig_dir):

            if file.endswith('.pdf'):
                if name_str in file:
                    if file[0] == name_str[0]:
                        highest = file.rsplit('.pdf', 1)[0]
                        highest_num += highest.rsplit('_', 1)[1]

        highest_num = [int(s) for s in highest_num]
        num = max(highest_num) + 1

    except ValueError:  # no file with the same parameters
        num = 0

    # read metadata of existing file to avoid recreating the same plot
    if num != 0:

        for i in range(num):

            fpname = os.path.join(fig_dir, name_str + '_' + str(i) + '.pdf')
            fp = open(fpname, 'rb')
            parser = PDFParser(fp)

            if (sys.version_info < (3, 0)):
                doc = PDFDocument(parser)

            else:
                doc = PDFDocument()
                parser.set_document(doc)
                doc.set_parser(parser)
                doc.initialize()

            fp.close()
            metadata = doc.info  # parameters used in previous runs

            # exit if the model has been run with same values before
            if (np.isclose(metadata[0]['PPFD'], md.PAR) and
               np.isclose(metadata[0]['minT'], md.minT) and
               np.isclose(metadata[0]['maxT'], md.maxT) and
               np.isclose(metadata[0]['PPT'], md.PPT) and
               np.isclose(metadata[0]['minVPD'], md.minVPD) and
               np.isclose(metadata[0]['maxVPD'], md.maxVPD) and
               np.isclose(metadata[0]['u'], md.u) and
               np.isclose(metadata[0]['maxleaf'], md.maxleaf) and
               np.isclose(metadata[0]['Ps'], md.Ps)):

                if force_write == 'no':
                    msg = ('This figure has already been generated on %s. '
                           'Please, check %s'
                           % (metadata[0]['CreationDate'], fpname))

                    if (sys.version_info < (3, 0)):
                        sys.tracebacklimit = 0

                    else:
                        sys.tracebacklimit = None

                    raise Exception(msg)

                if force_write == 'yes':
                    num = i
                    os.remove(fpname)
                    break

    return str(num)


def txt_info_box(md, fig, axis):

    """
    Writes a text box containing basic parameter information.

    Arguments:
    ----------
    md: class or pandas object
        parameters the model has been run with

    fig: matplotlib object
        the figure to add the text box to

    axis: matplotlib object
        axis on which to plot

    Returns:
    --------
    Plots the relevant text box on the axis.

    """

    txtstr = ["Parameters",
              "P$_{\mathrm{50}}$ = " + str(round(float(md.P50), 2)) + " MPa",
              "P$_{\mathrm{88}}$ = " + str(round(float(md.P88), 2)) + " MPa",
              r"$\frac{k_{\mathrm{max}}}{LAI}$ = " + str(round(float(md.kmax),
                                                               2)) +
              " mmol m$^{-2}$ s$^{-1}$ MPa$^{-1}$", "$PPFD$ = " +
              str(round(float(md.PPFD), 2)) + r" $\mu$mol m$^{-2}$ s$^{-1}$",
              "$VPD_{min}$ = " + str(round(float(md.minVPD), 2)) + " kPa",
              "$VPD_{max}$ = " + str(round(float(md.maxVPD), 2)) + " kPa"]

    t = axis.transData
    total_height = 0.

    for s in reversed(txtstr):

        if s != txtstr[0]:
            text = axis.text(0.05, 0.05, " " + s + " ", fontsize=8,
                             bbox={'facecolor': 'none', 'edgecolor': 'none',
                                   'pad': 0.3}, transform=t)
            text.draw(fig.canvas.get_renderer())

            if s == txtstr[len(txtstr)-1]:
                ex = text.get_window_extent()
                line_height = (ex.height) * 0.8

            if s == txtstr[1]:
                line_height = (ex.height) * 1.

            total_height += line_height

            t = transforms.offset_copy(text._transform, y=line_height,
                                       units='dots')

        if s == txtstr[0]:
            text = axis.text(0.05, 0.05, " " + s + " ", fontsize=8,
                             bbox={'facecolor': 'none', 'edgecolor': 'none'},
                             transform=t)
            text.draw(fig.canvas.get_renderer())
            ex = text.get_window_extent()
            t = transforms.offset_copy(text._transform, units='dots')

    # the position of the text box depends on the number of axes
    if len(fig.axes) == 4:  # one subplot
        ax_loc = axis.get_position()
        new_loc = total_height * (ax_loc.y1 - ax_loc.y0) / 10.
        ax_loc.y0 = new_loc
        ax_loc.y1 = 2. * new_loc
        ax_loc.x1 = ax_loc.x1 + 0.0015

    if len(fig.axes) == 6:  # two subplots
        ax_loc = axis.get_position()
        new_loc = total_height * (ax_loc.y0 - ax_loc.y1) / 10.
        ax_loc.y0 = 2. * new_loc - (ax_loc.y0 - ax_loc.y1)
        ax_loc.y1 = new_loc - 2.1 * (ax_loc.y0 - ax_loc.y1)
        ax_loc.x1 = ax_loc.x1 + 0.0015

    ax = fig.add_axes(ax_loc)
    ax.set_zorder(-1)
    ax.axis('off')

    r = patches.Rectangle((0., 0.), .001, .001, fill=False, edgecolor='none',
                          visible=False)
    handles = [r for e in txtstr[:len(txtstr)-1]]
    leg = ax.legend(handles, txtstr[:len(txtstr)-1], fontsize=6.8)
    plt.setp(leg.get_texts(), color='none')

    return fig


def set_fig_info(infodic, md):

    """
    Writes a text box containing basic parameter information.

    Arguments:
    ----------
    infodic: pdf object
        metadata info dictionary that goes with the figure

    md: class or pandas object
        parameters the model has been run with

    Returns:
    --------
    Sets the relevant metadata.

    """

    infodic['Title'] = __title__
    infodic['Version'] = __version__
    infodic['Author'] = __author__
    infodic['Contact'] = __email__
    infodic['Subject'] = 'Parameters used to run the models'
    infodic['P50'] = md.P50
    infodic['P88'] = md.P88
    infodic['kmax:LAI'] = md.kmax
    infodic['Ps'] = md.Ps
    infodic['PPFD'] = md.PPFD
    infodic['minT'] = md.minT
    infodic['maxT'] = md.maxT
    infodic['PPT'] = md.PPT
    infodic['minVPD'] = md.minVPD
    infodic['maxVPD'] = md.maxVPD
    infodic['u'] = md.u
    infodic['maxleaf'] = md.maxleaf
    infodic['CreationDate'] = datetime.now().strftime('%Y-%m-%d')

    return


def plt_intra_std_n_opt(fpname, df1, df2, psi_case):

    """
    Plots the carbon and water fluxes for the ProfitMax compared to the
    Control.

    Arguments:
    ----------
    fpname: string
        name (path) of the figure to create

    df1: pandas dataframe
        dataframe containing the input data

    df2: pandas dataframe
        dataframe containing the output data

    psi_case: int
        which ProfitMax case

    Returns:
    --------
    Plots the relevant figure.

    """

    # read params & met data
    md = FigInfo(df1)

    # fill the nans
    df1.fillna(value=0.)
    df2.fillna(value=0.)

    vplot = ['A(std)', 'A(psi%d)' % (psi_case), 'E(std)',
             'E(psi%d)' % (psi_case), 'Tair', 'precip']
    colours = ['#e5f5f9', '#99d8c9', '#2ca25f', '#fee8c8', '#fdbb84',
               '#e34a33', 'darkgrey', 'silver']
    labels = ['$A_{standard}$', r'A$_{\psi}$', '$E_{standard}$',
              r'E$_{\psi}$', '$T_{air}$', '$precip$']

    try:
        with PdfPages(fpname) as pdf:

            fig = plt.figure()

            ax1 = plt.subplot(211)
            ax1.tick_params(labelsize=7)

            ax2 = plt.subplot(212, sharex=ax1)
            ax2.tick_params(labelsize=7)

            # add right axes
            ax3 = fig.add_subplot(211, sharex=ax1, frameon=False)
            ax3.yaxis.tick_right()
            ax3.yaxis.set_label_position("right")
            ax3.tick_params(labelsize=7)

            ax4 = fig.add_subplot(212, sharex=ax2, frameon=False)
            ax4.yaxis.tick_right()
            ax4.yaxis.set_label_position("right")
            ax4.tick_params(labelsize=7)

            # add x axis for doy display, below ax2
            ax5 = ax2.twiny()
            ax5.xaxis.set_ticks_position("bottom")
            ax5.xaxis.set_label_position("bottom")
            ax5.spines['bottom'].set_position(("axes", -0.05))
            ax5.spines['bottom'].set_visible(False)
            ax5.tick_params(colors='w', labelcolor='k', labelsize=7)

            # text box with input parameters, below figure
            ax6 = fig.add_axes([0.12, -0.22, 0.3, 0.02])
            ax6.axis('off')

            # plot the data
            df2[vplot[:2]].plot(linewidth=1.5, color=colours[:2],
                                label=labels[:2], ax=ax1)
            df2[vplot[2:4]].plot(linewidth=1.5, color=colours[2:4],
                                 label=labels[2:4], ax=ax2)
            df1[vplot[4]].plot(linestyle=':', color=colours[6],
                               label=labels[4], ax=ax3)
            df1[vplot[5]].plot(kind='bar', width=0.8, color=colours[7],
                               label=labels[5], ax=ax4)

            # legend must display labels from both left & right axes
            handles1, __ = ax1.get_legend_handles_labels()
            handles3, __ = ax3.get_legend_handles_labels()
            handles1.insert(4, handles3[0])
            labels1 = labels[:2]
            labels1.insert(4, labels[4])
            ax1.legend(handles1, labels1, fontsize='small',
                       bbox_to_anchor=(0.2 + ax6.get_position().x1, -1.05 +
                                       ax6.get_position().y1 +
                                       ax6.get_position().y0))

            handles2, __ = ax2.get_legend_handles_labels()
            handles4, __ = ax4.get_legend_handles_labels()
            handles2.insert(4, handles4[0])
            labels2 = labels[2:4]
            labels2.insert(4, labels[5])
            ax2.legend(handles2, labels2, fontsize='small',
                       bbox_to_anchor=(0.4365 + ax6.get_position().x1,
                                       ax6.get_position().y1 +
                                       ax6.get_position().y0))

            # time series
            plt.minorticks_off()
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax3.get_xticklabels(), visible=False)
            plt.setp(ax4.get_xticklabels(), visible=False)

            freq_tick = int(len(df1) / 18)  # max 18 xticks for time
            doys = [datetime.strptime(str(e), '%j').strftime('%d-%m') for e in
                    range(int(md.doy), int(md.doy2) + 1)]
            timestep = df1['hod'][1] - df1['hod'][0]
            timeseries = pd.date_range(doys[0] + '-2018', periods=len(df1),
                                       freq=str(timestep) + 'H')
            timeseries = timeseries[::freq_tick]
            timeseries = [str(timeseries[e]) for e in range(len(timeseries))]
            timeseries = [timeseries[e][11:16] for e in range(len(timeseries))]

            xticks = ax2.get_xticks()
            ax2.set_xticks(xticks[::freq_tick])
            ax2.set_xticklabels(timeseries, fontsize=6, rotation='horizontal')

            doys = [datetime.strptime(str(e), '%j').strftime('%d %b') for e in
                    range(int(md.doy), int(md.doy2) + 1)]
            ax5.set_xticks(range(len(doys) * 2 + 1))

            for i in range(len(doys)):

                doys.insert(2 * i, '')

            ax5.set_xticklabels(doys, rotation='horizontal')

            # axes labels
            ax5.set_xlabel('time (h, date)', fontsize=7)
            ax1.set_ylabel(r'assimilation rate [$\mu$mol s$^{-1}$ m$^{-2}$]',
                           fontsize=7)
            ax3.set_ylabel(r'temperature [$^\circ$C]', fontsize=7)
            ax2.set_ylabel('transpiration [mmol s$^{-1}$ m$^{-2}$]',
                           fontsize=7)
            ax4.set_ylabel('precip [mm d$^{-1}$]', fontsize=7)

            # textbox
            fig = txt_info_box(md, fig, ax6)

            # insure that left axis is above right axis
            ax1.set_zorder(ax3.get_zorder() + 1)
            ax1.patch.set_visible(False)

            ax2.set_zorder(ax4.get_zorder() + 1)
            ax2.patch.set_visible(False)

            # set vertical gap between subplots
            plt.subplots_adjust(hspace=0.05)

            title = 'ProfitMax optimization between ' + doys[1] + ' and ' + \
                    doys[len(doys)-1]
            fig.suptitle(title, size=12, color='dimgray', y=0.95)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # set the file's metadata via the PdfPages object
            set_fig_info(pdf.infodict(), md)

    except Exception:
        os.remove(fpname)

    return


def plt_Farq_Col(fpname, df1, df2, df3, psi_case):

    """
    Plots the carbon and water fluxes for the Farquhar photosynthetic
    model vs the Collatz photosynthetic model.

    Arguments:
    ----------
    fpname: string
        name (path) of the figure to create

    df1: pandas dataframe
        dataframe containing the input data

    df2: pandas dataframe
        dataframe containing the output data

    df3: pandas dataframe
        dataframe containing the output data

    psi_case: int
        which ProfitMax case

    Returns:
    --------
    Plots the relevant figure.

    """

    # read params & met data
    md = FigInfo(df1)

    # fill the nans
    df1.fillna(value=0.)
    df2.fillna(value=0.)
    df3.fillna(value=0.)

    vplot = ['A(std)', 'A(psi%d)' % (psi_case), 'E(std)',
             'E(psi%d)' % (psi_case)]
    colours = ['#e5f5f9', '#99d8c9', '#fee8c8', '#fdbb84']

    labels = ['$A_{std}$', r'A$_{\psi}$', '$E_{std}$', r'E$_{\psi}$',
              'Farquhar', 'Collatz']

    try:
        with PdfPages(fpname) as pdf:

            fig = plt.figure()

            ax1 = plt.subplot(211)
            ax1.tick_params(labelsize=7)

            ax2 = plt.subplot(212, sharex=ax1)
            ax2.tick_params(labelsize=7)

            # add x axis for doy display, below ax2
            ax3 = ax2.twiny()
            ax3.xaxis.set_ticks_position("bottom")
            ax3.xaxis.set_label_position("bottom")
            ax3.spines['bottom'].set_position(("axes", -0.05))
            ax3.spines['bottom'].set_visible(False)
            ax3.tick_params(colors='w', labelcolor='k', labelsize=7)

            # text box with input parameters, to the right of figure
            ax4 = fig.add_axes([0.99295, 0.074, 0.3, 0.02])
            ax4.axis('off')

            # plot the data
            df2[vplot[:2]].plot(linestyle='--', linewidth=1.5,
                                color=colours[:2], label=labels[:2], ax=ax1)
            df3[vplot[:2]].plot(linestyle=':', linewidth=1.5,
                                color=colours[:2], label=labels[:2], ax=ax1)

            df2[vplot[2:]].plot(linestyle='--', linewidth=1.5,
                                color=colours[2:], label=labels[2:4], ax=ax2)
            df3[vplot[2:]].plot(linestyle=':', linewidth=1.5,
                                color=colours[2:], label=labels[2:4], ax=ax2)

            # legend must display labels from both left & right axes
            handles = [mlines.Line2D([], [], color=c, markersize=15)
                       for c in colours]
            handles.insert(7, mlines.Line2D([], [], linestyle='--', color='k',
                           markersize='15'))
            handles.insert(8, mlines.Line2D([], [], linestyle=':', color='k',
                           markersize='15'))

            ax1.legend(handles, labels, fontsize='small', loc='lower left',
                       bbox_to_anchor=(ax4.get_position().x0 +
                                       ax1.get_position().x0, 0.05))
            ax2.legend_.remove()

            # time series
            plt.setp(ax1.get_xticklabels(), visible=False)

            freq_tick = int(len(df1) / 18)  # max 18 xticks for time
            doys = [datetime.strptime(str(e), '%j').strftime('%d-%m') for e in
                    range(int(md.doy), int(md.doy2) + 1)]
            timestep = df1['hod'][1] - df1['hod'][0]
            timeseries = pd.date_range(doys[0] + '-2018', periods=len(df1),
                                       freq=str(timestep) + 'H')
            timeseries = timeseries[::freq_tick]
            timeseries = [str(timeseries[e]) for e in range(len(timeseries))]
            timeseries = [timeseries[e][11:16] for e in range(len(timeseries))]

            xticks = list(range(len(df1)))
            ax2.set_xticks(xticks[::freq_tick])
            ax2.set_xticklabels(timeseries, fontsize=6, rotation='horizontal')

            doys = [datetime.strptime(str(e), '%j').strftime('%d %b') for e in
                    range(int(md.doy), int(md.doy2) + 1)]
            ax3.set_xticks(range(len(doys) * 2 + 1))

            for i in range(len(doys)):

                doys.insert(2 * i, '')

            # axes labels
            ax3.set_xticklabels(doys, rotation='horizontal')
            ax3.set_xlabel('time (h, date)', fontsize=7)
            ax1.set_ylabel(r'assimilation rate [$\mu$mol s$^{-1}$ m$^{-2}$]',
                           fontsize=7)
            ax2.set_ylabel('transpiration [mmol s$^{-1}$ m$^{-2}$]',
                           fontsize=7)

            # textbox
            fig = txt_info_box(md, fig, ax4)

            # set vertical gap between subplots
            plt.subplots_adjust(hspace=0.05)

            title = 'ProfitMax optimization between ' + doys[1] + ' and ' + \
                    doys[len(doys)-1]
            fig.suptitle(title, size=13, color='dimgray', y=0.95)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # set the file's metadata via the PdfPages object
            set_fig_info(pdf.infodict(), md)

    except Exception:
        os.remove(fpname)

    return


def plt_intra_psi_opt(fpname, df1, df2):

    """
    Plots the two possible ProfitMax model solver cases for the carbon
    and water fluxes.

    Arguments:
    ----------
    fpname: string
        name (path) of the figure to create

    df1: pandas dataframe
        dataframe containing the input data

    df2: pandas dataframe
        dataframe containing the output data

    Returns:
    --------
    Plots the relevant figure.

    """

    # read params & met data
    md = FigInfo(df1)

    # fill the nans
    df1.fillna(value=0.)
    df2.fillna(value=0.)

    vplot = ['A(psi1)', 'A(psi2)', 'E(psi1)', 'E(psi2)', 'Tair', 'precip']
    colours = ['#01665e', '#5ab4ac', '#c7eae5', '#f6e8c3', 'darkgrey',
               'silver']
    labels = ['$case$ 1', '$case$ 2', '$T_{air}$', '$precip$']
    try:
        with PdfPages(fpname) as pdf:

            fig = plt.figure()

            ax1 = plt.subplot(211)
            ax1.tick_params(labelsize=7)

            ax2 = plt.subplot(212, sharex=ax1)
            ax2.tick_params(labelsize=7)

            # add right axes
            ax3 = fig.add_subplot(211, sharex=ax1, frameon=False)
            ax3.yaxis.tick_right()
            ax3.yaxis.set_label_position("right")
            ax3.tick_params(labelsize=7)

            ax4 = fig.add_subplot(212, sharex=ax2, frameon=False)
            ax4.yaxis.tick_right()
            ax4.yaxis.set_label_position("right")
            ax4.tick_params(labelsize=7)

            # add x axis for doy display, below ax2
            ax5 = ax2.twiny()
            ax5.xaxis.set_ticks_position("bottom")
            ax5.xaxis.set_label_position("bottom")
            ax5.spines['bottom'].set_position(("axes", -0.05))
            ax5.spines['bottom'].set_visible(False)
            ax5.tick_params(colors='w', labelcolor='k', labelsize=7)

            # text box with input parameters, below figure
            ax6 = fig.add_axes([0.12, -0.22, 0.3, 0.02])
            ax6.axis('off')

            # plot the data
            df2[vplot[:2]].plot(linewidth=1.5, color=colours[:2],
                                label=labels[:2], ax=ax1)
            df2[vplot[2:4]].plot(linewidth=1.5, color=colours[:2],
                                 label=labels[:2], ax=ax2)

            df1[vplot[4]].plot(linestyle=':', color=colours[3],
                               label=labels[3], ax=ax3)
            df1[vplot[5]].plot(kind='bar', width=0.8, color=colours[4],
                               label=labels[4], ax=ax4)

            # legend must display labels from both left & right axes
            handles, __ = ax1.get_legend_handles_labels()
            handles3, __ = ax3.get_legend_handles_labels()
            handles4, __ = ax4.get_legend_handles_labels()
            handles.insert(7, handles3[0])
            handles.insert(8, handles4[0])
            ax1.legend(handles, labels, fontsize=6,
                       bbox_to_anchor=(0.54, -1.05 + ax6.get_position().y1 +
                                       ax6.get_position().y0 / 1.94))
            ax2.legend_.remove()

            # time series
            plt.minorticks_off()
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax3.get_xticklabels(), visible=False)
            plt.setp(ax4.get_xticklabels(), visible=False)

            freq_tick = int(len(df1) / 18)  # max 18 xticks for time
            doys = [datetime.strptime(str(e), '%j').strftime('%d-%m') for e in
                    range(int(md.doy), int(md.doy2) + 1)]
            timestep = df1['hod'][1] - df1['hod'][0]
            timeseries = pd.date_range(doys[0] + '-2018', periods=len(df1),
                                       freq=str(timestep) + 'H')
            timeseries = timeseries[::freq_tick]
            timeseries = [str(timeseries[e]) for e in range(len(timeseries))]
            timeseries = [timeseries[e][11:16] for e in range(len(timeseries))]

            xticks = ax2.get_xticks()
            ax2.set_xticks(xticks[::freq_tick])
            ax2.set_xticklabels(timeseries, fontsize=6, rotation='horizontal')

            doys = [datetime.strptime(str(e), '%j').strftime('%d %b') for e in
                    range(int(md.doy), int(md.doy2) + 1)]
            ax5.set_xticks(range(len(doys) * 2 + 1))

            for i in range(len(doys)):

                doys.insert(2 * i, '')

            ax5.set_xticklabels(doys, rotation='horizontal')

            # axes labels
            ax5.set_xlabel('time (h, date)', fontsize=7)

            ax1.set_ylabel(r'assimilation rate [$\mu$mol s$^{-1}$ m$^{-2}$]',
                           fontsize=7)
            ax3.set_ylabel(r'temperature [$^\circ$C]', fontsize=7)

            ax2.set_ylabel('transpiration [mmol s$^{-1}$ m$^{-2}$]',
                           fontsize=7)
            ax4.set_ylabel('precip [mm d$^{-1}$]', fontsize=7)

            # textbox
            fig = txt_info_box(md, fig, ax6)

            # insure that left axis is above right axis
            ax1.set_zorder(ax3.get_zorder() + 1)
            ax1.patch.set_visible(False)

            ax2.set_zorder(ax4.get_zorder() + 1)
            ax2.patch.set_visible(False)

            # set vertical gap between subplots
            plt.subplots_adjust(hspace=0.05)

            title = 'ProfitMax optimizations between ' + doys[1] + ' and ' + \
                    doys[len(doys)-1]
            fig.suptitle(title, size=12, color='dimgray', y=0.95)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # set the file's metadata via the PdfPages object
            set_fig_info(pdf.infodict(), md)

    except Exception:
        os.remove(fpname)

    return

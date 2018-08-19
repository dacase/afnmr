#!/usr/bin/env python2
from __future__ import division

from argparse import ArgumentParser
from collections import OrderedDict
import matplotlib
import numpy as np
from scipy.stats import linregress
import sys

def calc_rmse(x, y):
   return np.sqrt(((x-y)**2).mean())

parser = ArgumentParser()
group = parser.add_argument_group('Files')
group.add_argument('-i', '--input-file', metavar='FILE',default=['hshifts.dat'],
                   dest='input_file', help='''Input file with shift data.
                   Default name is %(default)s. Data should be in 4 columns:
                   <Atom name> <observed shifts> <calculated shifts> <error>''',
                   nargs='*')
group.add_argument('-o', '--output-file', metavar='IMAGE', default=None,
                    dest='output_file', help='''File to save image in. Default
                    is to just show it on the screen.''')
group = parser.add_argument_group('Filtering options')
group.add_argument('-s', '--skip', dest='skipped', metavar='<ATOM>', nargs='*',
                    help='List of atom names to skip when plotting.',
                    default=[])
group.add_argument('-a', '--avg-adjust', dest='coil', default=False,
                   action='store_true', help='''Subtract a 'random coil' shift
                   from each point calculated as the average shift from all data
                   on the X-axis for that type of nucleus.''')
group.add_argument('--only', dest='only', metavar='<ATOM>', nargs='*',
                   help='List of atom names to (exclusively) plot.', default=[])
group = parser.add_argument_group('Plot Appearance Options')
group.add_argument('-t', '--title', dest='title', metavar='TITLE',
                   default=None, help='''Title of the plot(s)''', nargs='*')
group.add_argument('-x', '--x-label', dest='xlabel', metavar='X-AXIS-TITLE',
                   default='Observed Chemical Shift (ppm)', help='''Title of the
                   X-axis. Default is [%(default)s]''')
group.add_argument('-y', '--y-label', dest='ylabel', metavar='Y-AXIS-TITLE',
                   default='Calculated Chemical Shift (ppm)', help='''Title of
                   the Y-axis. Default is [%(default)s]''')
group.add_argument('--legend-columns', default=2, metavar='INT', dest='cols',
                   help='''Number of columns to divide legend into. Default is
                   %(default)s.''', type=int)
group.add_argument('--no-legend', dest='legend', default=True,
                   action='store_false', help='''Do not show a legend separating
                   nuclei types''')
group = parser.add_argument_group('Graph Properties')
group.add_argument('-r', '--range', dest='range', metavar='<FLOAT>',
                   nargs=2, default=None, help='''Chemical shift range to plot
                   (same for both axes)''')
group.add_argument('-f', '--linear-fit', dest='fit', default=False,
                   action='store_true', help='''Fit all of the data using a
                   linear regression. Also plots the best-fit line and prints
                   the fitting statistics to the plot''')
group.add_argument('--no-stats', dest='printstats', default=True,
                   action='store_false', help='''Do not print the linear
                   regression statistics on the plot''')
group.add_argument('--no-grid', dest='grid', default=True, action='store_false',
                   help='Turn off grid in plot')
group.add_argument('--r-label', dest='rlab', default=None, nargs='*',
                   metavar='<FLOAT>', help='''Location to put a label of the
                   correlation coefficient for the best-fit line''', type=float)
group.add_argument('--rrmse', dest='rlabrmse', default=False,
                   action='store_true', help='''Include the RMSE in the label
                   printing the correlation coefficient. Default is not to''')
group = parser.add_argument_group('Other Options:')
group.add_argument('-d', '--detailed-stats', dest='details', default=False,
                   action='store_true', help='''Compute RMSE and R^2 values for
                   individual families of protons.''')
group.add_argument('--xkcd', dest='xkcd', default=False, action='store_true',
                   help='Use XKCD-style plots')
group.add_argument('--scale-factor', dest='scale_factor', default=1.0,
                   type=float, help='''The scaling factor by which to increase
                   the size of the plot and its various attributes''')
group.add_argument('--legend-scale-factor', dest='legend_scale_factor',
                   default=1.0, type=float, help='''The scaling factor by which
                   to increase the size of the legend text.''')

opt = parser.parse_args()

if opt.title is None:
   opt.title = ['AF-NMR'] * len(opt.input_file)

if opt.output_file is not None:
   matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

if opt.xkcd:
    plt.xkcd()

symbols = ['o', 'x', '^', '*', 'h', '+', 'D', '8', 'v', '<', '>']
LS = len(symbols)
colors = ['r', 'b', 'g', 'k', 'm']
LC = len(colors)

factor = opt.scale_factor
legend_factor = opt.legend_scale_factor
ninp = len(opt.input_file)
thin = False
if ninp in (1, 6):
   fig = plt.figure(1, figsize=(15*factor,10*factor))
elif ninp == 2:
   fig = plt.figure(1, figsize=(10*factor, 5*factor))
   thin = True
elif ninp == 3:
   fig = plt.figure(1, figsize=(15*factor, 5*factor))
   thin = True
elif ninp == 4:
   fig = plt.figure(1, figsize=(10*factor, 10*factor))

if opt.rlab is not None and opt.fit:
   if len(opt.rlab) % 2 != 0:
      sys.exit('Cannot have an odd number of R-label coordinates.')
   if len(opt.rlab) != 2 * ninp and len(opt.rlab) != 2:
      sys.exit('Bad number of R-label coordinates (expected 2 or %d, got %d)' %
               (2*ninp, len(opt.rlab)))
   # Give every plot the same label location
   if len(opt.rlab) == 2:
      opt.rlab *= ninp

MAXROWS = 2
MAXCOLS = 3

FITLINE_LABEL_POSITIVE = """\
$y = %.2f x + %.2f$
$R^2 = %.3f$
$RMSE = %.2f$ ppm
"""

FITLINE_LABEL_NEGATIVE = """\
$y = %.2f x - %.2f$
$R^2 = %.3f$
$RMSE = %.2f$ ppm
"""

if not len(opt.input_file) in (1,2,3,4,6):
   sys.exit('Only support 1, 2, 3, 4, or 6 plots per figure!')

# Set some layout-specific sizes
LAYOUTS = {1 : 110, 2 : 120, 3 : 130, 4 : 220, 6 : 230}
ROWS = {1 : 1, 2 : 1, 3 : 1, 4 : 2, 6 : 2}

layout = LAYOUTS[len(opt.input_file)]
rows = ROWS[len(opt.input_file)]
fontsizes = dict(title=int(26-2*(len(opt.input_file)))*factor,
                 label=int(23-2*(len(opt.input_file)))*factor)
legendfont = FontProperties()
legendfont.set_size((16-len(opt.input_file))*legend_factor)
markersize = 8*factor - 0.6 * (len(opt.input_file)-1)

if opt.range is not None:
   opt.range = [float(opt.range[0]), float(opt.range[1])]

if opt.skipped and opt.only:
   sys.exit('Cannot specify both skipped atoms and only atoms!')

# Add a "ghost" plot encompassing the entire figure to implement a set of common
# axis labels for each plot. We need to turn off all lines associated with the
# blank plot, though
ax = fig.add_subplot(111)
for key in ax.spines:
    ax.spines[key].set_color('none')
ticksize = ax.get_xaxis().get_ticklabels()[0].get_fontsize() * factor
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off',
               labelsize=ticksize)
ax.set_xlabel(opt.xlabel, fontdict=dict(fontsize=fontsizes['label']))
ax.set_ylabel(opt.ylabel, fontdict=dict(fontsize=fontsizes['label']))

for ii, input_file in enumerate(opt.input_file):
   
   atnames = []
   obs = []
   calc = []
   calcerr = []

   with open(input_file, 'r') as f:
      for line in f:
         words = line.split()
         atnames.append(words[0])
         obs.append(float(words[1]))
         calc.append(float(words[2]))
         try:
            calcerr.append(float(words[3]))
         except IndexError:
            calcerr.append(0.0)
      obs = np.asarray(obs)
      calc = np.asarray(calc)
      calcerr = np.asarray(calcerr)
   
   # Now get the atom types
   datapoints = OrderedDict()
   if opt.details:
      nucdata = OrderedDict()
   
   for i in range(len(obs)):
      try:
         datapoints[atnames[i]].append(i)
      except KeyError:
         datapoints[atnames[i]] = [i]
   
   ax = fig.add_subplot(layout + ii + 1)
   pls = []
   if opt.grid: ax.grid(lw=1)
   if opt.range is not None:
      ax.set_xlim(opt.range)
      ax.set_ylim(opt.range)
   
   used_keys = []

   gbl_minmax_x = None

   if opt.fit:
      allx = np.ndarray(0)
      ally = np.ndarray(0)

   for i, key in enumerate(datapoints):
      # Skip over those atoms we don't want to plot
      if key in opt.skipped: continue
      if opt.only and not key in opt.only: continue
      used_keys.append(key)
      xdata = np.asarray([obs[j] for j in datapoints[key]])
      if gbl_minmax_x is None:
         gbl_minmax_x = [np.min(xdata), np.max(xdata)]
      else:
         gbl_minmax_x[0] = min(np.min(xdata), gbl_minmax_x[0])
         gbl_minmax_x[1] = max(np.max(xdata), gbl_minmax_x[1])
      ydata = np.asarray([calc[j] for j in datapoints[key]])
      ybars = np.asarray([calcerr[j] for j in datapoints[key]])
      if opt.coil:
         mean = xdata.mean()
         xdata -= mean
         ydata -= mean
      if opt.fit:
         allx = np.concatenate((allx, xdata))
         ally = np.concatenate((ally, ydata))
      if opt.details:
         nucdata[key] = (xdata, ydata)
   
      if opt.legend and len(datapoints.keys()) > 1:
         marker, color = symbols[i%LS], colors[i%LC]
      else:
         marker, color = symbols[0], colors[3]

      if not any(calcerr):
         pl, = ax.plot(xdata, ydata, marker=marker, color=color, ls='None',
                       markersize=markersize)
      else:
         pl, cl, bl = ax.errorbar(xdata, ydata, yerr=ybars, marker=marker,
                              color=color, ls='None', markersize=markersize)
      pls.append(pl)

   # Plot the diagonal line
   pts = (min(ax.xaxis.get_view_interval()[0], ax.yaxis.get_view_interval()[0]),
         max(ax.xaxis.get_view_interval()[1], ax.yaxis.get_view_interval()[1]))

   ax.set_xlim(pts)
   ax.set_ylim(pts)
   ax.plot(pts, pts, 'k-', lw=2)

   # If we want a fit, do a linear regression
   if opt.fit:
      m, b, r, p, stderr = linregress(allx, ally)
      rmse = calc_rmse(allx, ally)
      xpt = gbl_minmax_x
      ypt = (m * gbl_minmax_x[0] + b, m * gbl_minmax_x[1] + b)
      ax.plot(xpt, ypt, 'k--', lw=2)
      if rows == 1 and not thin:
         xt, yt = 0.05, 0.8
      else:
         xt, yt = 0.05, 0.65

      if opt.printstats:
         if b >= 0:
            FITLINE_LABEL = FITLINE_LABEL_POSITIVE
         else:
            FITLINE_LABEL = FITLINE_LABEL_NEGATIVE
         ax.text(xt, yt, FITLINE_LABEL % (m, abs(b), r**2, rmse),
                transform=ax.transAxes, fontsize=14)
      if opt.rlab is not None:
         if opt.rlabrmse:
            lab = '$R = %.3f$\n$RMSE = %.3f ppm$' % (r, rmse)
         else:
            lab = '$R = %.3f$' % r
         ax.text(opt.rlab[2*ii], opt.rlab[2*ii+1], lab,
                 fontsize=fontsizes['label'], transform=ax.transAxes)

      print 'Plot %d' % ii
      print '--------------------------------'
      print 'Linear regression:'
      print '     Slope = %.2f' % m
      print ' Intercept = %.2f' % b
      print ' Corr.Coef = %.3f' % r
      print '        2'
      print '       R   = %.3f' % (r * r)
      print '         p = %.4f' % p
      print ' Std. Err. = %.4f' % stderr
      print '      RMSE = %.4f' % rmse
      if opt.details:
         # Now print out the R^2 and RMSE for each type of nucleus
         print ''
         print 'Nucleus |  R^2   |  RMSE'
         print '--------+--------+-------'
         for key in nucdata:
            xdata, ydata = nucdata[key]
            m, b, r, p, stderr = linregress(xdata, ydata)
            print '   %-4s | %.4f | %.4f' % (key, r*r, calc_rmse(xdata, ydata))
         print ''
      print '--------------------------------'

   if len(pls) > 1 and opt.legend:
      if opt.printstats:
         loc = 4
      else:
         loc = 2
      leg = ax.legend(pls, [key.replace('_', ' ') for key in used_keys],
                     loc=loc, ncol=opt.cols, numpoints=1, prop=legendfont)
      leg.get_title().set_size(legendfont.get_size())
   ax.set_title(opt.title[ii], fontdict=dict(fontsize=fontsizes['title']))
   # Get the current tick sizes, then multiply it by the factor and set the tick
   # sizes
   ticksize = ax.get_xaxis().get_ticklabels()[0].get_fontsize() * factor
   ax.tick_params(labelsize=ticksize)

plt.tight_layout()

if opt.output_file is None:
   plt.show()
else:
   fig.savefig(opt.output_file)

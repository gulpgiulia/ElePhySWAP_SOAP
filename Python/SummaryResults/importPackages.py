#------------------------------------------------------
#  Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
#  Version: 1.0 - June 6, 2017
#-----------------------------------------------------

# ---- importPackages

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.ticker import LinearLocator, MultipleLocator, FormatStrFormatter

import numpy as np

import scipy.io as sio 
import scipy.stats as stats
import scipy.interpolate

from operator import itemgetter, attrgetter, methodcaller

from mpl_toolkits.mplot3d import Axes3D

import math

import os

import copy

import sys

import statsmodels
from statsmodels.stats import multitest

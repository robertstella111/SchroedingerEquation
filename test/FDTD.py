from Schroedinger import Schroedinger
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np


length = 20
b = 0.4
a = 0
obj = Schroedinger.FDTD()
obj.setNumPoints(length)
obj.setDomain(a,b)
obj.solve()



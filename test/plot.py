from Schroedinger import Schroedinger
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np


length = 160
b = 1
a = 0
obj = Schroedinger.Schroedinger1D()
obj.setNumPoints(length)
obj.setDomain(a,b)
obj.solve(100000)


def f1(n):
    return obj.getEigenvektor(int(n))

ykoord = obj.getEigenvektor(2)
xkoord = obj.getXKoord()   
print(len(xkoord))
print(len(ykoord))

init_n = 0
n = 0
fig, ax = plt.subplots()
line, = ax.plot(xkoord, f1(init_n))

ax.set_ylim(-1.05*np.sqrt(2/(b-a)), 1.05*np.sqrt(2/(b-a)))

fig.subplots_adjust(left=0.15, bottom = 0.15)

axamp = fig.add_axes([0.05, 0.25, 0.0225, 0.63])

amp_slider = Slider(
    ax=axamp,
    label="n",
    valmin=0,
    valmax=length-2,
    valinit=init_n,
    orientation="vertical",
    valfmt="%i"
)


def update(val):
    line.set_ydata(f1(int(amp_slider.val)))
    fig.canvas.draw_idle()

amp_slider.on_changed(update)


resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button1 = Button(resetax, '+1', hovercolor='0.975')


def add(event):
    if amp_slider.val < length -2:
        n = amp_slider.val + 1
        amp_slider.set_val(n)
        line.set_ydata(f1(int(n)))
        fig.canvas.draw_idle()

resetax2 = fig.add_axes([0.6, 0.025, 0.1, 0.04])
button2 = Button(resetax2, '-1', hovercolor='0.975')


def sub(event):
    if amp_slider.val > 0:
        n = amp_slider.val -  1
        amp_slider.set_val(n)
        line.set_ydata(f1(int(n)))
        fig.canvas.draw_idle()

button1.on_clicked(add)
button2.on_clicked(sub)

plt.show()

print(obj.getEigenValues())
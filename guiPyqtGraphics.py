__author__ = 'ruben'

# import pyqtgraph as pg
#
# ### create data to plot ######
# from numpy import *
# pi=3.1415
# X=linspace(-10,10,100)
# X1=linspace(-1,1,10);X11=linspace(-1,1,11)
# Y=2+sin(X1);Y1=2+sin(X);
# Y2=(cos(10*Y)/(X1+0.0131415))**2
# Y3=cos(10*Y1)/(X+0.0131415)
# Y4=4+sin(X)*cos(2*X)
# ###############################
#
# ##### plot two graphs in two windows ########
# pg.plot(X,Y1)
# # another graph
# pg.plot(X,Y1,pen=(255,0,255),symbol='+')
# # show
# pg.QtGui.QApplication.exec_()
#
# ### plot two graphs on the same window ####
# H=pg.plot(X,Y1, pen='b',symbolPen="w",symbol='t')
# H.plot(X,Y1, pen=(255,0,0),symbolPen="r",symbol="s")
# H.plot(1.1*X,sqrt(Y1),symbolPen="y")
# # show
# pg.QtGui.QApplication.exec_()
#
# ##### make a window with multiple graphs ####
# win= pg.GraphicsWindow(title="subplot window") # make the window
# p=win.addPlot(title="fig 1");p.plot(X,Y1)
# p=win.addPlot(title="fig 2");p.plot(X,Y4)
# win.nextRow() # go down one row in the window
# p=win.addPlot(title="fig 3");p.plot(X,Y3)
# curve=pg.PlotCurveItem(X11,Y2,stepMode=True,fillLevel=0, brush=(0, 0, 255, 80))
# p=win.addPlot(title="fig 4");p.addItem(curve)
#
# # show
# pg.QtGui.QApplication.exec_()

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

app = QtGui.QApplication([])
win = QtGui.QWidget()

layout = QtGui.QGridLayout()
win.setLayout(layout)

iv = pg.ImageView()
pw = pg.PlotWidget()

img = (np.random.random((300, 200)) * 255).astype('uint8')

data = np.zeros(((255,) + img.shape))
for i in xrange(255):
    data[i, :, :] = np.where(img > i, img, 0)

iv.setImage(data, xvals=np.linspace(0, 255, data.shape[0]))
layout.addWidget(iv, 0, 1)

prof = np.array([img[:, x].sum() for x in xrange(data.shape[2])])


pw.setFixedWidth(100)
pw.setYLink(iv.view)
pw.plot(prof, np.arange(prof.shape[0]))
pw.invertY()
layout.addWidget(pw, 0, 0)


def update_pw():
    global pw, iv, data
    prof = np.array([data[iv.currentIndex, :, x].sum() for x in xrange(data.shape[2])])
    pw.clear()
    pw.plot(prof, np.arange(prof.shape[0]))

iv.timeLine.sigPositionChanged.connect(update_pw)


win.show()


if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
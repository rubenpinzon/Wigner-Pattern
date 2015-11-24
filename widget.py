from PyQt4 import QtCore, QtGui


class Window(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.resize(800, 600)
        self.setWindowTitle(QtGui.QApplication.translate("self", "self", None, QtGui.QApplication.UnicodeUTF8))
        self.setDockOptions(QtGui.QMainWindow.AnimatedDocks)
        self.centralwidget = QtGui.QWidget(self)
        self.centralwidget.hide()
        self.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(self)
        self.setStatusBar(self.statusbar)
        self.dock1Widget = QtGui.QDockWidget(self)
        self.dock1Widget.setFeatures(QtGui.QDockWidget.AllDockWidgetFeatures)
        self.dock1Widget.setWindowTitle(
            QtGui.QApplication.translate("self", "dock1", None, QtGui.QApplication.UnicodeUTF8))
        self.dock1WidgetContents = QtGui.QWidget()
        self.dock1Widget.setWidget(self.dock1WidgetContents)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dock1Widget)
        self.dock2Widget = QtGui.QDockWidget(self)
        self.dock2Widget.setFeatures(QtGui.QDockWidget.AllDockWidgetFeatures)
        self.dock2Widget.setWindowTitle(
            QtGui.QApplication.translate("self", "dock2", None, QtGui.QApplication.UnicodeUTF8))
        self.dock2WidgetContents = QtGui.QWidget()
        self.dock2Widget.setWidget(self.dock2WidgetContents)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dock2Widget)
        self.dock3Widget = QtGui.QDockWidget(self)
        self.dock3Widget.setFeatures(QtGui.QDockWidget.AllDockWidgetFeatures)
        self.dock3Widget.setWindowTitle(
            QtGui.QApplication.translate("self", "dock3", None, QtGui.QApplication.UnicodeUTF8))
        self.dock3WidgetContents = QtGui.QWidget()
        self.dock3Widget.setWidget(self.dock3WidgetContents)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dock3Widget)
        self.dock4Widget = QtGui.QDockWidget(self)
        self.dock4Widget.setFeatures(QtGui.QDockWidget.AllDockWidgetFeatures)
        self.dock4Widget.setWindowTitle(
            QtGui.QApplication.translate("self", "dock4", None, QtGui.QApplication.UnicodeUTF8))
        self.dock4WidgetContents = QtGui.QWidget()
        self.dock4Widget.setWidget(self.dock4WidgetContents)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dock4Widget)


if __name__ == '__main__':
    import sys
    import cv2

    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    app.exec_()

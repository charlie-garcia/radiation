# ------------------------------------------------------
# -------------------- mplwidget.py --------------------
# ------------------------------------------------------
from PyQt5.QtWidgets import*
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure
    
class widget_plot1(QWidget):
    
    def __init__(self, parent = None):

        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvas(Figure())
        self.canvas.figure.set_facecolor('#EDEDEE')

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.canvas.axes.set_facecolor('#EDEDEE')
        self.canvas.axes.tick_params(color='#C0C0C0', labelcolor='k')
        self.canvas.axes.tick_params(axis='both', which='major', labelsize=7)
        self.canvas.axes.tick_params(axis='both', which='minor', labelsize=6)
        self.setLayout(vertical_layout)
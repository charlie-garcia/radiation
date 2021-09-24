# ------------------------------------------------------
# -------------------- mplwidget.py --------------------
# ------------------------------------------------------
from PyQt5.QtWidgets import*
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure
    
class widget_radiation_modes(QWidget):
    
    def __init__(self, parent = None):

        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvas(Figure())
        self.canvas.figure.set_facecolor('#EDEDEE')
        self.canvas.figure.suptitle('Radiation Modes')
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        
        n_modes = 7
        self.canvas.axes = [0]*n_modes

        for iax in range(0,n_modes):
            
            self.canvas.axes[iax] = self.canvas.figure.add_subplot(1, n_modes, iax+1, projection='3d')
            self.canvas.axes[iax].set_facecolor('#EDEDEE')
            self.canvas.axes[iax].tick_params(color='#C0C0C0', labelcolor='k')
            self.canvas.axes[iax].tick_params(axis='both', which='major', labelsize=5)
            self.canvas.axes[iax].tick_params(axis='both', which='minor', labelsize=4)

        self.setLayout(vertical_layout)
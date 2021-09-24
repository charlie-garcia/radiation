# ------------------------------------------------------
# ---------------------- main.py -----------------------
# ------------------------------------------------------
import sys, os, random, time
from PyQt5.QtWidgets import*
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np
import scipy.interpolate as inter

from mf.radiation import radiation, mode_excitation
from mf.dcm import RadModes, MatrixR, distance3
from mf.plots import create_plate2D, create_plate2D_gui

class MatplotlibWidget(QMainWindow):
    def __init__(self):
        
        QMainWindow.__init__(   self)
        loadUi("radiation_gui.ui",self)
        self.setWindowTitle("Radiation")

        self.pushButton_calculate.clicked.connect(self.compute_ARM)

        self.doubleSpinBox_young.editingFinished.connect(self.update_graph)
        self.doubleSpinBox_density.editingFinished.connect(self.update_graph)
        self.doubleSpinBox_heigth.editingFinished.connect(self.update_graph)
        self.doubleSpinBox_damping.editingFinished.connect(self.update_graph)
        self.doubleSpinBox_force.editingFinished.connect(self.update_graph)
        self.doubleSpinBox_px.editingFinished.connect(self.update_graph)
        self.doubleSpinBox_py.editingFinished.connect(self.update_graph)

        # self.checkBox_modes.valueChanged.connect(self.plot_mode)
        self.checkBox_modes.stateChanged.connect(self.plot_mode)
        self.scroll_frequency.sliderReleased.connect(self.plot_mode)
        self.scroll_frequency.actionTriggered.connect(self.plot_mode)

        self.freq_up_button.clicked.connect(self.plot_up)
        self.freq_down_button.clicked.connect(self.plot_down)

        self.doubleSpinBox_px.setMaximum(float(self.lineEdit_length.text()))
        self.doubleSpinBox_py.setMaximum(float(self.lineEdit_width.text()))

        self.addToolBar(NavigationToolbar(self.widget_plot1.canvas, self))
        self.addToolBar(NavigationToolbar(self.widget_plot2.canvas, self))

        self.compute_ARM()


    def compute_ARM(self):
        t = time.time()
        
        self.get_spin_values()

        self.color_plot = ['steelblue'] # , 'orange', 'green']    

        M = 400
        N = int(np.sqrt(M))
        dx = self.Lx/N
        dy = self.Ly/N

        cs, self.x, self.y = create_plate2D(self.Lx,self.Ly,N, 'no')
        xx = cs[:,0]
        yy = cs[:,1]

        A = dx*dy * np.ones([N**2,1])

        #  Analytical Solution : Modes in SS Plate
        Nm = self.spinBox_Nm.value();           Nn = self.spinBox_Nn.value()
        m  =  np.r_[1:Nm+1];                    n  = np.r_[1:Nn+1]
        km = np.array([m*np.pi/self.Lx]);       kn = np.array([n*np.pi/self.Ly]);

        alpha = np.sin(km.T*xx)         # m * nx
        beta  = np.sin(kn.T*yy)          # n * ny

        self.mass = self.rho*self.Lx*self.Ly*self.h
        self.wmn  = np.sqrt(self.D/(self.rho*self.h)) * (km.T**2+kn**2)

        # Radiation
        Nf    = 450                                                          # N frequencies
        
        self.c0  = 343.5                                                       # [m/s] c air  
        self.rho0  = 1.23                                              

        self.f0 = np.min(self.wmn)/(2*np.pi)
        self.fc = self.c0**2/(2*np.pi)*np.sqrt(self.rho*self.h/self.D)
        self.f= np.logspace( 0, np.log10(self.fc*1.5), Nf )  

        omega = 2*np.pi*self.f                                                     
        k0    = omega/self.c0  

        rij   = distance3(cs, cs)

        # Pack variables
        self.N = N
        self.M = M
        self.Nf = Nf
        self.Nm = Nm
        self.Nn = Nn
        self.alpha =alpha
        self.beta = beta
        self.m = m
        self.n = n
        self.km =km
        self.kn =kn
        self.omega =omega
        self.A =A

        self.R = [0]*Nf
        # # Store ARM
        for iif in range(0,Nf):
            self.R[iif] = MatrixR(self.rho0, self.c0, k0[iif], self.A, rij)

        elapsed = time.time() - t  
        print('ARM computed, elapsed time : %.2f sec' %(elapsed))
        
        # ============================  Plot radiation modes ==================================

        x = np.linspace(0,self.Lx,self.N)
        y = np.linspace(0,self.Ly,self.N)
        self.X, self.Y = np.meshgrid(x,y)
        
        self.n_modes = 7
        
        f_indx = 0

        nn = int(f_indx)
        mode_i, efficiency_i = RadModes( self.R[nn], nn, self.n_modes)

        for iax in range(0,self.n_modes):
            
            my_col = cm.jet( mode_i[iax]/np.amax( mode_i[iax] ))
            self.widget_radiation_modes.canvas.axes[iax].clear()
            self.widget_radiation_modes.canvas.axes[iax].plot_surface(self.X, self.Y, mode_i[iax], facecolors = my_col) 
            self.widget_radiation_modes.canvas.axes[iax].set_aspect('auto', adjustable='box')
            self.widget_radiation_modes.canvas.axes[iax].set_title('$\lambda_{%.f}$ : %.2e'%(iax+1, np.real(efficiency_i[iax])), fontsize=15)

        
        self.widget_radiation_modes.canvas.draw()

        self.update_graph()

    def update_graph(self):

        t = time.time()
        
        self.get_spin_values()

        self.mass = self.rho*self.Lx*self.Ly*self.h

        self.wmn  = np.sqrt(self.D/(self.rho*self.h)) * (self.km.T**2+self.kn**2)

        # compute Radiation
        self.sigma_force, self.W, self.mean_v = radiation(self)
        
        f0 = np.min(self.wmn)/(2*np.pi)
        fc = self.c0**2/(2*np.pi)*np.sqrt(self.rho*self.h/self.D)
        self.fres = (np.sort(self.wmn.reshape(-1)/(2*np.pi))).reshape(self.Nm*self.Nn,1)

        id_f0 = (np.abs(self.fres[0]*.6 - self.f)).argmin()                         # find closest freq to resonant mode
        id_fc = (np.abs(fc*2 - self.f)).argmin()                                    # find closest freq to resonant mode

        self.ff = np.round((self.wmn[0,0])/(2*np.pi),0)
        self.scroll_frequency.setMaximum(len(self.f))
       
        # ============================  Plot 1 : Radiation Efficiency ==================================
        self.widget_plot1.canvas.axes.clear()
        
        self.widget_plot1.canvas.axes.loglog(self.f, np.real(self.sigma_force), label='$\sigma$')        
        self.widget_plot1.canvas.axes.vlines(f0, 0, 3, linestyles='dotted', colors=[.6,.6,.6], label='_nolegend_')
        self.widget_plot1.canvas.axes.text(f0-15, 10**(-4.5), '$f_0$', fontsize=15, color=[.6,.6,.6])
        self.widget_plot1.canvas.axes.vlines(fc, 0, 3, linestyles='dotted', colors=[.6,.6,.6],  label='_nolegend_')
        self.widget_plot1.canvas.axes.text(fc+15, 10**(-4.5), '$f_c$', fontsize=15, color=[.6,.6,.6])
        self.widget_plot1.canvas.axes.set_ylim([10**(-5), 3])
        self.widget_plot1.canvas.axes.set_xlabel('f (Hz)')
        self.widget_plot1.canvas.axes.set_ylabel('$\sigma$')
        self.widget_plot1.canvas.axes.set_title('Radiation effiiciency \n $\sigma$ = $\mathcal{P}_a}$/$\mathcal{P}_m$')

        self.line_plot1 = self.widget_plot1.canvas.axes.axvline(self.ff, np.min(np.real(self.sigma_force)), 3, visible=False, linestyle='dotted', color='#ef2929')

        self.widget_plot1.canvas.axes.legend()
        self.widget_plot1.canvas.draw()

        # ============================  Plot 2: Sound Power ==================================
        self.widget_plot2.canvas.axes.clear()
        
        self.cpow = 10*np.log10(np.real(self.W)/(10e-12))
        self.widget_plot2.canvas.axes.semilogx(self.f, self.cpow ,  color='#cc0099', label ='$L_w$')

        self.widget_plot2.canvas.axes.set_xlim([self.f[id_f0], self.f[id_fc]])     
        self.widget_plot2.canvas.axes.set_ylim([np.min(self.cpow[id_f0:id_fc])*0.95, np.max(self.cpow[id_f0:id_fc])*1.05])    
        self.widget_plot2.canvas.axes.set_xlabel('f (Hz)')
        self.widget_plot2.canvas.axes.set_ylabel('Lw (dB ref $10^{-12} W$)')
        self.widget_plot2.canvas.axes.vlines(fc, 0, 200, linestyles='dotted', colors=[.6,.6,.6],  label='_nolegend_')
        self.widget_plot2.canvas.axes.set_title('Sound power,  $\mathcal{P}_a = v_e^H \mathbf{R} v_e$') 
        
        self.line_plot2 = self.widget_plot2.canvas.axes.axvline(self.ff, -100, 100, visible=False, linestyle='dotted', color='#ef2929')

        self.widget_plot2.canvas.axes.legend()

        self.widget_plot2.canvas.draw()

        # ============================  Plot 3:  mean squared velocity ==================================
        self.widget_plot3.canvas.axes.clear()

        self.cmsv =10*np.log10( np.abs(self.mean_v) )
        self.widget_plot3.canvas.axes.semilogx(self.f, self.cmsv, color='#990000', label='$<v^2>$')
        self.widget_plot3.canvas.axes.set_xlim([self.f[id_f0], self.f[id_fc]])          
        self.widget_plot3.canvas.axes.set_ylim([np.min(np.real(self.cmsv[id_f0:id_fc]))*1.05,   np.max(np.real(self.cmsv[id_f0:id_fc]))*0.95])      
        self.widget_plot3.canvas.axes.set_xlabel('f (Hz)')
        self.widget_plot3.canvas.axes.set_ylabel('mean velocity (dB ref 1m/s/N)')
        self.widget_plot3.canvas.axes.vlines(fc, -100, 0, linestyles='dotted', colors=[.6,.6,.6],  label='_nolegend_')

        self.line_plot3 = self.widget_plot3.canvas.axes.axvline(self.ff, -100, 100, visible=False, linestyle='dotted', color='#ef2929')

        self.widget_plot3.canvas.axes.legend()

        self.widget_plot3.canvas.draw()

        # ============================  Plot Response to excitation ==================================
        self.auxp = 0
        self.auxm = 0

        my_map = plt.get_cmap('seismic')
        
        self.widget_excitation.canvas.axes.clear()
        plot_square_grid = 'no'

        cs, x, y = create_plate2D_gui(self.Lx,self.Ly,self.N, plot_square_grid, self.widget_excitation.canvas.axes)
        self.widget_excitation.canvas.axes.plot(self.px, self.py, 'x', c='g', ms=15, lw=15)
        self.cs = cs
        response = mode_excitation(self)
        mode_response = inter.griddata(self.cs, np.real(response), (self.X,self.Y), method = 'cubic').reshape(len(self.X), len(self.Y))
        mode_response = np.nan_to_num(mode_response)
        
        self.surf_excitation = self.widget_excitation.canvas.axes.pcolorfast(self.x, self.y, mode_response, cmap = my_map)


        self.widget_excitation.canvas.axes.set_aspect('equal', adjustable='box')
        self.widget_excitation.canvas.axes.set_title('Response to force')
        self.widget_excitation.canvas.draw()
        
        self.plot_mode()
        elapsed = time.time() - t  
        print('Radiation modes computed, elapsed time : %.2f sec' %(elapsed))


    def get_spin_values(self):
        self.E   = self.doubleSpinBox_young.value()*1e9         # =====================================
        self.rho = self.doubleSpinBox_density.value()           # =
        self.Lx  = float(self.lineEdit_length.text())           # = 
        self.Ly  = float(self.lineEdit_width.text())            # =
        self.h   = self.doubleSpinBox_heigth.value()*1e-3       # =
        self.eta = self.doubleSpinBox_damping.value()           # = verifiy class name of doubleSpinBox
        
        self.F  = self.doubleSpinBox_force.value()              # =
        self.px = self.doubleSpinBox_px.value()                 # =
        self.py = self.doubleSpinBox_py.value()                 # =====================================

        self.nu   = 0.23
        self.rho0 = 1.23
        self.D    = self.E*self.h**3/(12*(1-self.nu**2))

    def plot_mode(self):

        if self.checkBox_modes.isChecked() == True:
            self.ff = self.f[self.scroll_frequency.value()]
            response = mode_excitation(self)

            # Better resolution for plate deformation
            X_, Y_ = np.meshgrid(np.linspace(0,self.Lx,self.N), 
                                 np.linspace(0,self.Ly,self.N))

            mode_response = inter.griddata(self.cs, np.real(response), (X_,Y_), method = 'cubic').reshape(len(X_), len(Y_))
            mode_response = np.nan_to_num(mode_response)


            f_indx = np.argwhere(self.f >= self.ff)[0]
            nn = int(f_indx)
            mode_i, efficiency_i = RadModes( self.R[nn], nn, self.n_modes)

            for iax in range(0,self.n_modes):
                my_col = cm.jet( mode_i[iax]/np.amax( mode_i[iax] ))

                self.widget_radiation_modes.canvas.axes[iax].clear()
                self.widget_radiation_modes.canvas.axes[iax].plot_surface(self.X, self.Y, mode_i[iax], facecolors = my_col) 
                self.widget_radiation_modes.canvas.axes[iax].set_aspect('auto', adjustable='box')
                self.widget_radiation_modes.canvas.axes[iax].set_title('$\lambda_{%.f}$ : %.2e'%(iax+1, np.real(efficiency_i[iax])), fontsize=15)

            self.line_plot1.set_xdata(self.ff)
            self.line_plot2.set_xdata(self.ff)
            self.line_plot3.set_xdata(self.ff)

            self.surf_excitation.set_data(np.real(mode_response))
            self.surf_excitation.set_clim(np.min(mode_response), np.max(mode_response))

            self.line_plot1.set_visible(True)
            self.line_plot2.set_visible(True)
            self.line_plot3.set_visible(True)

            self.surf_excitation.set_visible(True)

            # # Plot
            self.widget_plot1.canvas.draw()
            self.widget_plot2.canvas.draw()
            self.widget_plot3.canvas.draw()

            self.widget_radiation_modes.canvas.draw()
            self.widget_excitation.canvas.draw()

            self.label_frequency.setText(QtCore.QCoreApplication.translate("MainWindow", "<html><head/><body><p><span style=\" color:#ef2929;\">\
                                                                           "+  '{:.2f}'.format(self.ff) +"</span></p></body></html>"))

        else:
            # self.renew()
            self.line_plot1.set_visible(False)
            self.line_plot2.set_visible(False)
            self.line_plot3.set_visible(False)
            self.surf_excitation.set_visible(False)

            self.widget_plot1.canvas.draw()
            self.widget_plot2.canvas.draw()
            self.widget_plot3.canvas.draw()
            self.widget_excitation.canvas.draw()

    def plot_up(self):
        self.ff = self.f[self.scroll_frequency.value()]
        id_res = (np.abs(self.fres - self.ff)).argmin()                             # find closest freq value
        idx = (np.abs(self.fres[id_res + self.auxp] - self.f)).argmin()             # find closest freq value
        self.scroll_frequency.setValue(idx)
        self.plot_mode()   
        self.auxp =  1
        

    def plot_down(self):
        self.ff = self.f[self.scroll_frequency.value()]
        id_res = (np.abs(self.fres - self.ff)).argmin()                             # find closest freq value
        idx = (np.abs(self.fres[id_res - self.auxm] - self.f)).argmin()             # find closest freq value
        
        self.scroll_frequency.setValue(idx)
        self.plot_mode()  
        self.auxm = 1
        
app = QApplication([])
window = MatplotlibWidget()
window.show()
app.exec_()
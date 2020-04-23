import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['text.usetex'] = True

class animation(object):
    '''Basic Animation Class

    Keyword arguments:
    nrows -- σ^2 is the variance
    ncols -- Mean (i.e. center)
    x,y   -- Arrays, or list/tuples of arrays
    nt    -- number of timesteps
    L     --
    '''
    def __init__(self,nrows,ncols,x,y,nt,L):
        self.nrows = nrows
        self.ncols = ncols
        self.x     = x
        self.y     = y
        self.nt    = nt
        self.L     = L

    def labels(self,xlabels,ylabels):
        # if (type(line) == list) | (type(line) == tuple):
        #     print('yerry')

        self.xlabels = xlabels
        self.ylabels = ylabels

    def animate(self):
        from matplotlib import animation, rc
        rc('animation', html='jshtml')
        plt.rcParams['font.size'] = 14

        fig, ax = plt.subplots(1,1,figsize=(12,6))
        ax.set_xlim(self.x[0], self.x[-1] )
        ax.set_ylim(0, 17.2)

        ax.set_ylabel(self.ylabels)
        ax.set_xlabel(self.xlabels)

        line1, = ax.plot([], [], 'x', lw=2, color='r',label='Crank-Nicolson')
        line2, = ax.plot([], [], lw=2, color='darkslategrey',label='Analytical')

        line = [line1,line2]
        ax.legend(loc=1)

        def animater(i):
            global x,u
            line[0].set_data(self.x, self.y[0][i,:])
            line[1].set_data(self.x, self.y[1][i,:])
            return line

        anim = animation.FuncAnimation(fig, animater,
                                       frames=range(0,self.nt), interval=40, blit=True)
        plt.close()

        return anim

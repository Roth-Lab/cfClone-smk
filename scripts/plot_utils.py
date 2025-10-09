from matplotlib import pyplot as plt


def setup_plot(font_size=None):
    """
    Set up matplotlib in a consistent way between plots.
    """
    plt.rcParams['font.family'] = ' sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'Nimbus Sans', 'Liberation Sans', 'DejaVu Sans', 'Arial']
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['grid.color'] = 'darkgrey'
    if font_size is not None:
        plt.rcParams["font.size"] = font_size


def setup_axes(ax):
    """
    Set up axes in a consistent way between plots.  Spines are at left and bottom, offset by 10pts.
    Other spines not visible.  Ticks only for left and bottom.
    """
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    ax.xaxis.grid(True, which="major", linestyle=':')
    ax.yaxis.grid(True, which="major", linestyle=':')

    return ax

import numpy as np
import pylab
pylab.ion()

def get_fig(fig_num, some_data, some_labels):

    fig = pylab.figure(fig_num,figsize=(8,8),frameon=False)
    ax = fig.add_subplot(111)
    ax.set_ylim([0.1,0.8]); ax.set_xlim([0.1, 0.8]);
    ax.set_title("Quarterly Stapler Thefts")
    ax.pie(some_data, labels=some_labels, autopct='%1.1f%%', shadow=True);
    return fig

my_labels = ("You", "Me", "Some guy", "Bob")

# To ensure first plot is always made.
do_plot = 1; num_plots = 0;

while do_plot:
    num_plots = num_plots + 1;
    data = np.random.rand(1,4).tolist()[0]

    fig = get_fig(num_plots,data,my_labels)
    fig.canvas.draw()
    pylab.draw()

    print "Close any of the previous plots? If yes, enter its number, otherwise enter 0..."
    close_plot = raw_input()

    if int(close_plot) > 0:
        pylab.close(int(close_plot))

    print "Create another random plot? 1 for yes; 0 for no."
    do_plot = raw_input();

    # Don't allow plots to go over 10.
    if num_plots > 10:
        do_plot = 0

pylab.show()

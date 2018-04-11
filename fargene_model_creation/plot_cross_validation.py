import numpy as np
import matplotlib.pyplot as plt
import matplotlib
plt.switch_backend('agg')

def plot_cross_validation(scores,sensitivity,specificity,figfile):
    plt.figure(figsize=(8,6), dpi=80)
    plt.subplot(111)

    plt.plot(scores, sensitivity, color="blue", linewidth=2.5, linestyle="-",label="Sensitivity")

    plt.plot(scores, specificity, color="green", linewidth=2.5, linestyle="-",label="1-Specificity")

    # Set x limits
    plt.xlim(min(scores),max(scores))
    plt.xlabel('Domain score')
    plt.ylabel('Identified sequences [fraction]')

    #plt.xticks(np.linspace(-4,4,9,endpoint=True))

    plt.ylim(0,1.01)

    # Set y ticks
    #plt.yticks(np.linspace(-1,1,5,endpoint=True))

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.legend(loc='upper right',frameon=False)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)

    plt.savefig(figfile)
    #plt.show()

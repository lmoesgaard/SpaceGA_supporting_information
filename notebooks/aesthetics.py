from matplotlib import pyplot as plt

def remove_spines(ax):
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    return ax

def remove_sticks(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    return ax

def invisible(ax):
    ax = remove_spines(ax)
    ax = remove_sticks(ax)
    ax.patch.set_alpha(0.0) 
    ax.set_xlim(0, 1)
    return ax

def remove_some_spines(ax):
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return ax

def remove_all_spines(ax):
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    return ax

def remove_sticks(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    return ax

def fix(ax):
    ax = remove_some_spines(ax)
    ax = remove_sticks(ax)
    return ax

def multiplot(y, x, xsize, ysize):
    # Make 6x4 subplots: top and left axis + 5 methods x 3 proteins
    fig, axes = plt.subplots(y+1, x+1, figsize=(xsize, ysize), width_ratios=[0] + [4]*x)

    # Format left axis
    for i in range(y+1):
        ax = axes[i, 0]
        ax = remove_spines(ax)
        ax = remove_sticks(ax)
        if i > 0:
            ax.spines['left'].set_visible(True)
            ax.plot([0], [1])
    
    # Format top axis
    for k in range(1, x+1):
        ax = axes[0, k]
        ax = remove_spines(ax)
        ax = remove_sticks(ax)
        ax.plot([0, 20], [0,0], c="k") # x-axis

    # Format main grid
    for i in range(1, x+1):
        for k in range(1, y+1):
            ax = axes[k, i]
            ax = remove_spines(ax)
            if k != y:
                ax = remove_sticks(ax)
            else:
                ax.spines['bottom'].set_visible(True)
                ax.set_yticks([])

    # Format bottom axis
    for i in range(1, x+1):
        k = y
        ax = axes[k, i]

    
    fig.tight_layout(pad=0.0)
    return fig, axes
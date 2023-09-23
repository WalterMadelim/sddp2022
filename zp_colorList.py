from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def plot_colortable(colors, sort_colors=True, emptycols=0):

    cell_width = 100
    cell_height = 22
    swatch_width = 48
    margin = 1

    # Sort colors by hue, saturation, value and name.
    if sort_colors is True:
        by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))),
                         name)
                        for name, color in colors.items())
        names = [name for hsv, name in by_hsv]
    else:
        names = list(colors)

    n = len(names)
    ncols = 4 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-margin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.add_patch(
            Rectangle(xy=(swatch_start_x, y-9), width=swatch_width,
                      height=18, facecolor=colors[name], edgecolor='0.7')
        )

    return fig

cl_cloR   = 138/256,33 /256,66/256
cl_lotusR = 207/256,131/256,185/256
cl_leaveG = 96 /256,124/256,56/256
cl_fenceG = 66 /256,100/256,93/256
cl_hairB  = 32 /256,25/256,25/256

color1 = (0.1647, 0.7176, 0.6157)  # 海绿
color2 = (0.2667, 0.7804, 0.9569,)  # 海蓝
color3 = (0.9059, 0.6627, 0.1255,)  # 脏橙
color4 = (0.8784, 0.7882, 0.2431,)  # 暗黄
color5 = (0.7686, 0.3529, 0.1804,)  # 血红
color6 = (0.6902, 0.6784, 0.9765)  # 淡紫色
color7 = (0.3922, 0.4078, 0.8549)  # 紫色
color8 = (0.4784, 0.8000, 0.4392,)  # 草绿色
color9 = (0.0863, 0.7843, 0.8353,)  # 电吉他青
color10 = (0.9961,0,0.5882,) # 桃红
color11 = (0.5412,0.1059,0.6431,) # 紫罗兰
color12 = (0.36862,1,0.003921) # 亮绿
color14 = (0,0.40784,0.70588) # intel蓝
color15 = (0.46274,0.7255,0) # Nvidia绿

dct = {'$P^{CHP}$':cl_cloR,'$P^{CS,d}$':color6,'$P^{pur}$':color7,'$P^{load}$':color2,'$P^{CS,c}$':color3,'$P^{PM}$':cl_lotusR,'$Q^{CHP}$':color1,'$Q^{load}$':color12,'$Q^{PM}$':color15,'$H^{CHP}$':color10,'$H^{pur}$':color11,'$H^{load}$':color14}
plot_colortable(dct, sort_colors=False)
plt.show()
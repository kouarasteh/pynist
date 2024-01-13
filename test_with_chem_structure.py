# stdlib
import datetime
import sys

# 3rd party
from pyms.Spectrum import MassSpectrum  # type: ignore
from matplotlib import rc
import numpy as np
# this package
import pyms_nist_search
import matplotlib.pyplot as plt
from mendeleev import element
from matplotlib import transforms
def color_element_text(element_list, x, y, ax):
    txts = element_list
    colors = []
    for elem in element_list:
        if elem.isdigit():
            color = colors[-1]
        elif elem == 'H':
            color = '#000000'
        else:
            color = element(elem).cpk_color
        colors.append(color)
    rainbow_text(x,y,txts,colors,ax,va='center', ha='center')


def rainbow_text(x,y,ls,lc,ax,**kw):
    """
    Take a list of strings ``ls`` and colors ``lc`` and place them next to each
    other, with text ls[i] being shown in color lc[i].

    This example shows how to do both vertical and horizontal text, and will
    pass all keyword arguments to plt.text, so you can set the font size,
    family, etc.
    """
    t = ax.transData
    fig = plt.gcf()

    text = plt.text(x,y,''.join(ls),color='white')
    text.set_bbox(dict(facecolor='white', alpha=1.0, edgecolor='white'))

    #horizontal version
    for s,c in zip(ls,lc):
        text = plt.text(x,y,s,color=c, transform=t, **kw)
        text.draw(fig.canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(text._transform, x=2*ex.width, units='dots')


# #vertical version
    # for s,c in zip(ls,lc):
    #     text = plt.text(x,y,s+" ",color=c, transform=t,
    #             rotation=90,va='bottom',ha='center',**kw)
    #     text.draw(fig.canvas.get_renderer())
    #     ex = text.get_window_extent()
    #     t = transforms.offset_copy(text._transform, y=ex.height, units='dots')


def draw_structure(stdata):
    print(stdata.keys())
    print(stdata)
    print(stdata['bond_type_list'])
    plt.figure(figsize=(10,10))
    ax = plt.gca()
    ax.set_xlim((-100,100))
    ax.set_ylim((-100,100))

    #ax.scatter(stdata['xl_list'],stdata['yl_list'])
    xl = stdata['xl_list'][:stdata['num_points']]
    yl = stdata['yl_list'][:stdata['num_points']]

    for i in range(0,stdata['num_points'],2):
        print(i, xl[i],yl[i])
        print('bond_type', stdata['bond_type_list'][i])
        ax.plot(xl[i:i+2],yl[i:i+2], color='black')

    # add aromatic circles
    if stdata['num_circs'] >0:
        for i in range(stdata['num_circs']):
            x_center,y_center = stdata['xc_list'][i],stdata['yc_list'][i]
            c = plt.Circle((x_center,y_center),10, fill=False, color='black')
            ax.add_patch(c)

    # add heteratom text
    if stdata['num_strs'] >0:
        for i in range(stdata['num_strs']):
            x_text, y_text = stdata['xs_list'][i], stdata['ys_list'][i]
            txt = stdata['str_list'][i]
            color_element_text(txt, x_text,y_text, ax)
            # t = plt.text(x_text,y_text, ''.join(txt), va='center', ha='center')
    plt.show()


if sys.platform == "win32":
    FULL_PATH_TO_MAIN_LIBRARY = r"C:\home\domdf\Python\mainlib"
    FULL_PATH_TO_USER_LIBRARY = r"C:\home\domdf\Python\01 GitHub Repos\pynist\tests\MoNA-export-GC-MS_Spectra"
    FULL_PATH_TO_WORK_DIR = r"C:\home\domdf\Python\00 Projects\pynist"
else:
    FULL_PATH_TO_MAIN_LIBRARY = "/home/domdf/Python/mainlib"
    FULL_PATH_TO_USER_LIBRARY = "/home/domdf/Python/01 GitHub Repos/pynist/tests/MoNA-export-GC-MS_Spectra"
    FULL_PATH_TO_WORK_DIR = "/home/domdf/Python/00 Projects/pynist"

search = pyms_nist_search.Engine(
        FULL_PATH_TO_MAIN_LIBRARY,
        pyms_nist_search.NISTMS_MAIN_LIB,
        FULL_PATH_TO_WORK_DIR,
        debug=True,
        )

# search = pyms_nist_search.Engine(
# 		FULL_PATH_TO_USER_LIBRARY,
# 		pyms_nist_search.NISTMS_USER_LIB,
# 		FULL_PATH_TO_WORK_DIR,
# 		debug=True,
# 		)

mz_int_pairs = [
        (27, 138),
        (28, 210),
        (32, 59),
        (37, 70),
        (38, 273),
        (39, 895),
        (40, 141),
        (41, 82),
        (50, 710),
        (51, 2151),
        (52, 434),
        (53, 49),
        (57, 41),
        (59, 121),
        (61, 73),
        (62, 229),
        (63, 703),
        (64, 490),
        (65, 1106),
        (66, 932),
        (67, 68),
        (70, 159),
        (71, 266),
        (72, 297),
        (73, 44),
        (74, 263),
        (75, 233),
        (76, 330),
        (77, 1636),
        (78, 294),
        (84, 1732),
        (87, 70),
        (88, 86),
        (89, 311),
        (90, 155),
        (91, 219),
        (92, 160),
        (93, 107),
        (101, 65),
        (102, 111),
        (103, 99),
        (104, 188),
        (113, 107),
        (114, 120),
        (115, 686),
        (116, 150),
        (117, 91),
        (126, 46),
        (127, 137),
        (128, 201),
        (129, 73),
        (130, 69),
        (139, 447),
        (140, 364),
        (141, 584),
        (142, 279),
        (143, 182),
        (152, 37),
        (153, 60),
        (154, 286),
        (166, 718),
        (167, 3770),
        (168, 6825),
        (169, 9999),
        (170, 1210),
        (171, 85),
        ]

mass_list = []
intensity_list = []
for mass, intensity in mz_int_pairs:
    mass_list.append(mass)
    intensity_list.append(intensity)

mass_spec = MassSpectrum(mass_list, intensity_list)

start_time = datetime.datetime.now()
print("Performing Full Search")

hit_list = search.full_search_with_ref_data(mass_spec)

for hit_no, (hit, ref_data) in enumerate(hit_list):
    print(f"Hit {hit_no}")
    print(hit)
    print(ref_data)
    print(ref_data.mass_spec)
    draw_structure(ref_data.structure_data)

    # reference_data = search.get_r#eference_data(hit.spec_loc)
    # print(reference_data.mass_spec == ref_data.mass_spec)
    # print(reference_data == ref_data)

print(f"Completed Full Search in {(datetime.datetime.now() - start_time).total_seconds()}")

############

start_time = datetime.datetime.now()
print("Performing Quick Search")

hit_list = search.spectrum_search(mass_spec, n_hits=10)

for hit_no, hit in enumerate(hit_list):
    print(f"Hit {hit_no}")
    print(hit)
    print()

print(f"Completed Quick Search in {(datetime.datetime.now() - start_time).total_seconds()}")

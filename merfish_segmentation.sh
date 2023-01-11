import os
import numpy as np
import cv2
import glob
import math
from scipy.spatial import KDTree
from shapely.geometry import Polygon
from copy import deepcopy
from tensorflow.keras.preprocessing import image

from cellpose import utils, io
from cellpose import models, io
from cellpose import plot

import merlin
from merlin.util import imagereader
from merlin.util import dataportal

from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import make_outline_overlay
from deepcell.utils.plot_utils import create_rgb_image

import PIL
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300



## define functions
def outlines_list(masks):
    """ get outlines of masks as a list to loop over for plotting """
    outpix=[]
    for n in np.unique(masks)[1:]:
        mn = masks==n
        if mn.sum() > 0:
            contours = cv2.findContours(mn.astype(np.uint8), mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE)
            contours = contours[-2]
            cmax = np.argmax([c.shape[0] for c in contours])
            pix = contours[cmax].astype(int).squeeze()
            if len(pix)>4:
                outpix.append(pix)
            else:
                outpix.append(np.zeros((0,2)))
    return outpix

def outlines_to_text(base, outlines):
    with open(base + '_cp_outlines.txt', 'w') as f:
        for o in outlines:
            xy = list(o.flatten())
            xy_str = ','.join(map(str, xy))
            f.write(xy_str)
            f.write('\n')


def read_coordinate(file_name):
    polygons = []
    with open(file_name, 'r') as f:
        for line in f:
            line = line.split(',')
            coordinates = [(int(line[2*i]), int(line[2*i+1])) for i in range(len(line)//2)]
            if not coordinates:
              continue
            else:
              polygons.append(coordinates)
    return polygons


def mean(l):
    return sum(l)/len(l)


def center(polygon):
    return mean([x[0] for x in polygon]), mean([x[1] for x in polygon])


def length(vector):
    return math.sqrt(vector[0]**2 + vector[1]**2)


def radius(polygon):
    """
    radius of circumcircle
    """
    center_point = center(polygon)
    return max([length((x[0]-center_point[0], x[1]-center_point[1])) for x in polygon])


def build_kd_tree(polygons, **kwargs):
    kd_tree = KDTree(np.array([center(x) for x in polygons]), **kwargs)
    return kd_tree


def combine_polygons(polygons_1, polygons_2, thres=0.5):
    new_polygons = list()
    if not polygons_1:
      new_polygons = polygons_2
      return  new_polygons
    else:
      kd_tree = build_kd_tree(polygons_1)
      r = max([radius(x) for x in polygons_1])
      new_polygons = deepcopy(polygons_1)
      for polygon in polygons_2:
          r0 = radius(polygon)
          indices = kd_tree.query_ball_point(center(polygon), r0+r)
          p0 = Polygon(polygon).convex_hull
          area = p0.area
          intersection_flag = False
          for idx in indices:
              p1 = Polygon(polygons_1[idx]).convex_hull
              intersection_area = p1.intersection(p0).area
              if intersection_area > thres * area:
                  intersection_flag = True
                  break
          if intersection_flag is False:
              new_polygons.append(polygon)
      return new_polygons


def write_polygon(polygons, filename):
    with open(filename, 'w') as f:
        for item in polygons:
            line=""
            for i, j in enumerate(item):

                line= line + ','.join(map(str, j))
                if i+1 < len(item):
                    line = line + ','
                else: 
                    break
            f.write("%s\n" % line)  
            
def combine_polygons_onlyifinDAPI(polygons_nuclei, polygons_wholecell):
    new_polygons = list()
    if not polygons_nuclei:
      return  new_polygons
    else:
      kd_tree = build_kd_tree(polygons_nuclei)
      r = max([radius(x) for x in polygons_nuclei])

      for polygon in polygons_wholecell:
          r0 = radius(polygon)
          indices = kd_tree.query_ball_point(center(polygon), r0+r)
          p0 = Polygon(polygon).convex_hull
          area = p0.area
          intersection_flag = False
          for idx in indices:
              p1 = Polygon(polygons_nuclei[idx]).convex_hull
              intersection_area = p1.intersection(p0).area
              if intersection_area > 0:
                  intersection_flag = True
                  break
          if intersection_flag is False:
              continue
          if intersection_flag is True:
              new_polygons.append(polygon)
      return new_polygons
            
    
    
## set variables
#'/storage/chenseq/MERFISH/raw_data/12022022_mouse_retina_VA45_XC/Data/' '/storage/chenseq/MERFISH/raw_data/12162022_mouse_retina_regeneration_VA45_XC/Data/'
#'/storage/chenseq/MERFISH/cell_segmentation/12022022_mouse_retina_VA45_XC/' '/storage/chenseq/MERFISH/cell_segmentation/12162022_mouse_retina_regeneration_VA45_XC/'

datadir = '/storage/chenseq/MERFISH/raw_data/111022_mouse_retina_VA45/Data/'
outdir = '/storage/chenseq/MERFISH/cell_segmentation/111022_mouse_retina_VA45/'

## create output dir for cellpose, mesmer, and filtered mesmer output files
os.makedirs(outdir+'/cellpose_cb2only_cb2dapi', exist_ok = True)
os.makedirs(outdir+'/mesmer_cb2_c20_cl2', exist_ok = True)
os.makedirs(outdir + '/mesmer_wholecell_nuclei_intersept', exist_ok = True)
os.makedirs(outdir + '/cellpose_mesmer_overlap_v2', exist_ok = True)

## loop through individual FOVs - segmentation
FilenamesList = glob.glob(datadir + 'hal-config-749z7-638z7-546z7-477z9-405z7' + '*.dax')
print(FilenamesList)

filecount = 0

for filename in FilenamesList:    
    os.chdir(outdir+'/mesmer_cb2_c20_cl2')
    ## import dax files 
    ## need to have inf + dax files in the same folder
    dataPortal = dataportal.LocalDataPortal(datadir)
    daxPortal = dataPortal.open_file(filename)
    a = imagereader.DaxReader(daxPortal)
    
    cb1_stack_list = range(14, 21)
    cb2_stack_list = range(7, 14)
    cb3_stack_list = range(0, 7)
    dapi_stack_list = range(30, 37)
    polyT_stack_list = range(21, 28)

    cb_1_image_data = []
    cb_2_image_data = []
    cb_3_image_data = []
    dapi_image_data = []
    polyT_image_data = []

    ## append Z stacks to each cell boundary stainging + DAPI
    for j in cb1_stack_list:
        this_frm = a.load_frame(j)
        cb_1_image_data.append(this_frm)
    for j in reversed(cb2_stack_list):
        this_frm = a.load_frame(j)
        cb_2_image_data.append(this_frm)
    for j in cb3_stack_list:
        this_frm = a.load_frame(j)
        cb_3_image_data.append(this_frm)
    for j in dapi_stack_list:
        this_frm = a.load_frame(j)
        dapi_image_data.append(this_frm)
    for j in reversed(polyT_stack_list):
        this_frm = a.load_frame(j)
        polyT_image_data.append(this_frm)

        
    # Calculate blended CB2 image and contrast for cellpose
    cb_2_dst = cb_2_image_data[0]

    for i in range(0,2):
        if i == 0:
            pass
        else:
            alpha = 1.0/(i + 1)
            beta = 1.0 - alpha
            cb_2_dst = cv2.addWeighted(cb_2_image_data[i], alpha, cb_2_dst, beta, 0.0)
 
    clahec = cv2.createCLAHE(clipLimit = 30)
    cb2_blended = clahec.apply(cb_2_dst)

    ## Contrast CB2 and DAPI for mesmer
    clahe = cv2.createCLAHE(clipLimit = 20, tileGridSize=(8, 8))
    clahed = cv2.createCLAHE(clipLimit = 2, tileGridSize=(8, 8))
    
    cb2_final_img = clahe.apply(cb_2_image_data[1])
    #dapi_final_img = clahed.apply(dapi_image_data[0])
    dapi_final_img = clahed.apply(dapi_image_data[4])
    
    ## convert 16 bits to 8 bits
    cb2_final_img8 = (cb2_final_img/256).astype('uint8')
    dapi_final_img8 = (dapi_final_img/256).astype('uint8')

    im = np.stack((dapi_final_img8, cb2_final_img8), axis=-1)
    im = np.expand_dims(im,0)
    
       
    ## run cellpose using blended CB2
    model = models.Cellpose(gpu=True, model_type='cyto2')
    channels = [0,0]
    masks, flows, styles, diams = model.eval(cb2_blended, diameter=None, channels=channels)
    #file1 = str(os.path.splitext(filename.split("/")[-1])[0] + '_cb2' + '_cL30_Z01blended')
    file1 = str(outdir+'/cellpose_cb2only_cb2dapi/' + os.path.splitext(filename.split("/")[-1])[0] + '_cb2_cL30_Z01blended')
    outlines = utils.outlines_list(masks)
    io.outlines_to_text(file1, outlines)
    
    print('\n Cellpose segmenting ' + str(filecount) + ' out of ' + str(len(FilenamesList)) + '\n') 
        
        
    ## run mesmer using CB2 and DAPI images
    app = Mesmer()
    segmentation_predictions = app.predict(im, image_mpp=0.1667)
    segmentation_predictions_nu = app.predict(im, image_mpp=0.1667, compartment='nuclear')

    outlines = outlines_list(segmentation_predictions[0,...,0])
    outlines_nu = outlines_list(segmentation_predictions_nu[0,...,0])
    
    ## write mesmer outlies to files
    file = str(outdir+'/mesmer_cb2_c20_cl2/' + os.path.splitext(filename.split("/")[-1])[0] + '_cb2_Z1_cL20_cL2_mesmer')
    file_nu = str(outdir+'/mesmer_cb2_c20_cl2/' + os.path.splitext(filename.split("/")[-1])[0] + '_cb2_Z1_cL20_cL2_mesmer_nuclei')  
    outlines_to_text(file, outlines)    
    outlines_to_text(file_nu, outlines_nu)

    print('Mesmer segmenting ' + str(filecount) + ' out of ' + str(len(FilenamesList)) + '\n') 

    filecount += 1

    ## can be used to plot mesmer prediction images
    '''
    rgb_images = create_rgb_image(im, channel_colors=['blue', 'green'])
    overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions)
    overlay_data_nu = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions_nu)
    
    fig, ax = plt.subplots(1, 2, figsize=(15, 15))
    ax[0].imshow(rgb_images[0, ...])
    ax[1].imshow(overlay_data[0, ...])
    ax[0].set_title('Raw data')
    ax[1].set_title('Whole-cell Predictions')
    plt.savefig(file + '_wholecell_prediction_plots.png', dpi = 300)
    plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(15, 15))
    ax[0].imshow(rgb_images[0, ...])
    ax[1].imshow(overlay_data_nu[0, ...])
    ax[0].set_title('Raw data')
    ax[1].set_title('Nuclei Predictions')
    plt.savefig(file_nu + '_nuclei_prediction_plots.png', dpi = 300)
    plt.close()
    '''

    ### filter out mesmer segmentation that does not have nuclei segmentation
    poly1 = read_coordinate(file_nu + '_cp_outlines.txt')    
    poly2 = read_coordinate(file + '_cp_outlines.txt')
    combined_mask = list()
    combined_mask = combine_polygons_onlyifinDAPI(poly1, poly2)    
    write_polygon(combined_mask, str(outdir + '/mesmer_wholecell_nuclei_intersept/' + os.path.splitext(filename.split("/")[-1])[0] + 
                                     "_cb2_Z1_cL20_wholecell_nuclei_intercept.txt"))

    
    ### use cellpose semgentation as base and fill in missing areas with mesmer segmentation
    poly_cellpose = read_coordinate(str(file1 + '_cp_outlines.txt'))
    poly_mesmer_filtered = read_coordinate(str(outdir + '/mesmer_wholecell_nuclei_intersept/' + os.path.splitext(filename.split("/")[-1])[0] + 
                                     "_cb2_Z1_cL20_wholecell_nuclei_intercept.txt"))
    combined_mask = list()
    combined_mask = combine_polygons(poly_cellpose, poly_mesmer_filtered, thres=0.2)
    write_polygon(combined_mask, str(outdir + '/cellpose_mesmer_overlap_v2/' + os.path.splitext(filename.split("/")[-1])[0] + 
                                      "_cb2_cellpose_mesmer_wholecell_nuclei_intercept.txt"))


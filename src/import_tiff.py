
"""
Load TIFF stack (exported from lif microscopy file by imagej)
"""

import sys, os
import numpy as np
from xml import etree as et
import matplotlib.pyplot as plt

from tifffile import TiffWriter, TiffFile, imsave

def load_tiff(fn, seq=0):
    with TiffFile(fn) as tiff:

        data = tiff.asarray()#colormapped=False)

#        print 'len(tiff)', len(tiff)

        data = np.squeeze(data)

#        print tiff.imagej_metadata

        try:
            unit = tiff.pages[0].tags['resolution_unit'].value
            print('unit', unit)
            u_sc = 0.00001 if unit==3 else 1.0
        except KeyError:
            u_sc = 1.0

        try:
            x_sp = tiff.pages[0].tags['XResolution'].value
        except KeyError:
            try:
                x_sp = tiff.pages[0].tags['x_resolution'].value            
            except KeyError:
                x_sp = (1, 1)


        try:
            y_sp = tiff.pages[0].tags['YResolution'].value
        except KeyError:
            try:
                y_sp = tiff.pages[0].tags['y_resolution'].value            
            except KeyError:
                y_sp = (1, 1)

        try:
            z_sp = tiff.imagej_metadata.get('spacing', 1.0)
        except AttributeError:
            z_sp = 1.0

    print(x_sp, y_sp, z_sp)
    if x_sp[0]==0:
        x_sp=(1.0,x_sp[1])
    if y_sp[0]==0:
        y_sp=(1.0,y_sp[1])

    x_sp = float(x_sp[1])/x_sp[0]*u_sc
    y_sp = float(y_sp[1])/y_sp[0]*u_sc
#    return np.transpose(data, (1,2,0)), (x_sp, y_sp, z_sp)
    print(data.shape, data.dtype)
    if len(data.shape)==2:
        data = data[np.newaxis,:,:]

#    if len(data.shape)==4:
#        data = np.max(data, axis=3)
    return data, (x_sp, y_sp, z_sp)

def write_tiff(fn, A, spacing):
#    A = np.transpose(A, (2, 0, 1))
    imsave(fn, A[np.newaxis, :, np.newaxis, :, :, np.newaxis], imagej=True,
                    resolution=(1.0/spacing[0], 1.0/spacing[1]),
                metadata={'spacing': spacing[2], 'unit': 'micron'})  



def load_tiff_stack(filename, im_num=None, verbose=True):
    """Load 3D TIFF from file

    Args:
        filename (string): Path of image file
        im_num (int): Im_Num of stack in lif file.

    Returns:
        array (numpy.ndarray): Array of image data
        spacing (tuple): Voxel spacing (in microns) in the i, j, k directions.
    """
    return load_tiff(filename)
        


def load_tiff_multi(filename, im_num=None, verbose=True):
    """Load 3D TIFF from file

    Args:
        filename (string): Path of image file
        im_num (int): Im_Num of stack in lif file.

    Returns:
        array (numpy.ndarray): Array of image data
        spacing (tuple): Voxel spacing (in microns) in the i, j, k directions.
    """
    return load_tiff(filename)



def load_tiff_vec(filename, im_num=None, verbose=True):
    """Load 3D TIFF from file

    Args:
        filename (string): Path of image file
        im_num (int): Im_Num of stack in lif file.

    Returns:
        array (numpy.ndarray): Array of image data
        spacing (tuple): Voxel spacing (in microns) in the i, j, k directions.
    """
    return load_tiff(filename)



def load_tiff_f(filename, im_num=None, verbose=True):
    """Load 3D TIFF from file

    Args:
        filename (string): Path of image file
        im_num (int): Im_Num of stack in lif file.

    Returns:
        array (numpy.ndarray): Array of image data
        spacing (tuple): Voxel spacing (in microns) in the i, j, k directions.
    """
    return load_tiff(filename)


if __name__=="__main__":
    image3d, spacing = load_tiff_stack(sys.argv[1])
    write_tiff(sys.argv[2], image3d, spacing)

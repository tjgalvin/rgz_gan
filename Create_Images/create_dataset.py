#!/usr/bin/env python
"""Script to read in and convert FITS images to numpy binaries, and create
a mask containing the location of the RGZ IR host
"""
import pandas as pd
import numpy as np 
import astropy.units as u
import pickle
from tqdm import tqdm
from astropy.io import fits
from astropy.coordinates import SkyCoord

np.random.seed(42) # Order is needed in the Universe


FIRST_PIX = 1.8*u.arcsecond # Pixel size of FIRST survey. Square pixels
FIRST_FWHM = 5*u.arcsecond / FIRST_PIX
FIRST_SIG = FIRST_FWHM / 2.355

def Gaussian(yx: tuple, cen: tuple, sigma: float):
    """Generate a Gaussian with an offset relative to the center of
    the image. 
    
    Arguments:
        yx {tuple} -- Output of meshgrid of the y and x coordinates
        cen {tuple} -- Position of the IR relative to the center of the image
        sigma {float} -- Size of the 2D Gaussian
    """
    y, x = yx
    y0, x0 = cen

    # Do the meshing
    A = 1 / (2*sigma**2)
    eq =  np.exp(-A*((x-x0)**2 + (y-y0)**2)) #Gaussian
    
    return eq


def make_ir_host(center_pos: SkyCoord, ir_pos: SkyCoord, shape: tuple):
    """Create a mask highlighting the IR host position
    
    Arguments:
        center_pos {SkyCoord} -- Center position of the current image
        ir_pos {SkyCoord} -- IR host position
        shape {tuple} -- Image dimension size to make mask to
    """
    offsets = center_pos.spherical_offsets_to(ir_pos)

    cen_xy = np.array(shape) // 2

    # ra increases right to left. Arrays do not.
    dx = -(offsets[0].to(u.arcsecond) / FIRST_PIX)
    dy = offsets[1].to(u.arcsecond) / FIRST_PIX
    
    yx = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    mask = Gaussian(yx, (cen_xy[0] + dy, cen_xy[1] + dx), FIRST_SIG)

    return mask


def get_fits(f: str):
    """Read in the FITS data for FIRST

    f {str} -- The filename to read in
    """
    with fits.open(f) as hdu1:
        img = hdu1[0].data.copy()

    return img


def main(files: pd.DataFrame, out_path: str, *args, 
        first_path: str='Images/first',
        wise_path: str='Images/wise_reprojected', 
        **kwargs):
    """Run the preprocessing on the set of files
    
    Arguments:
        files {pd.DataFrame} -- Dataframe with catalogue of objects/filenames
        out_path {str} -- Name of the file to create
        first_path {str} -- Relative path containing the first fits images
        wise_path {str} -- Relative path containing the wise reprojected fits images

    Raises:
        Exception -- Catch all used for removing files that fail preprocessing
    """

    shape = get_fits(f"{first_path}/{files['filename'].values[0]}").shape
    height, width = shape

    print(f'Derived height, width: {height}, {width}')

    failed  = [] 
    to_dump = {'first':[], 'wise':[], 'host':[], 'row':[]}

    for (index, row) in tqdm(df.iterrows()):
        f = row['filename']
        try:
            # If no reliable IR position, move a long
            if row['consensus.ir_ra'] == -99:
                continue

            img_first = get_fits(f"{first_path}/{f}")
            img_wise  = get_fits(f"{wise_path}/{f}")

            radio_pos   = SkyCoord(ra=row['radio.ra']*u.deg, dec=row['radio.dec']*u.deg, frame='icrs')
            ir_host_pos = SkyCoord(ra=row['consensus.ir_ra']*u.deg, dec=row['consensus.ir_dec']*u.deg, frame='icrs')        

            mask = make_ir_host(radio_pos, ir_host_pos, (height, width))

            to_dump['first'].append(img_first.astype('f'))
            to_dump['wise'].append(img_wise.astype('f'))
            to_dump['host'].append(mask.astype('f'))
            to_dump['row'].append(row)
        
        except ValueError as ve:
            print(f"{f}: {ve}")
            failed.append(row)
        except FileNotFoundError as fe:
            print(f"{f}: {fe}")
            failed.append(row)

        except Exception as e:
            raise e

    # Pickle up the objects
    with open(out_path, 'wb') as of:
        pickle.dump(out_path, to_dump)

    return 


if __name__ == '__main__':

    df = pd.read_csv('FIRST_Cata_Images.csv')

    main(df, 'rgz_images.pkl')
    
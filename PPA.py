# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 22:40:05 2014

@author: Themos Tsikas, Jack Richmond
"""

import sys
import os
import shutil
import time
import configparser
import tkinter.dialog
import numpy
import scipy
import tkinter
import tkinter.messagebox
import tkinter.filedialog
import optparse
import requests
import subprocess
import threading

import seestar_api


from nova_module import *

from PIL import ImageFont, ImageDraw, Image

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.io import fits
from astropy import wcs

from urllib.request import urlopen


def help_f():
    '''
    Our help window
    '''
    tkinter.messagebox.showinfo("Help", "Still to come...")


def about_f():
    '''
    our about window
    '''
    tkinter.messagebox.showinfo('About',
                          'PhotoPolarAlign v2.0.0 \n' +
                          'Copyright © 2024 Themos Tsikas, ' +
                          'Jack Richmond, ' +
                          'JC Lesaint, ' +
                          'Antony Linares')

def scale_frm_wcs(fn):
    hdu = fits.open(fn)
    head = hdu[0].header
    return scale_frm_header(head)

def parity_frm_header(head):
    '''
    look in the plate-solution header for the parity information
    '''
    try:
        # nova's wcs files have the parity in the comments
        comments = head['COMMENT']
        size = (len(comments))
        for i in range(0, size):
            if comments[i][0:6] == 'parity':
                tkns = comments[i].split(' ')
                return int(tkns[1])
        # if not found, return 1
        return 1
    except KeyError:
        return 1

    
def scale_frm_header(head):
    '''
    look in the plate-solution header for the scale information
    '''
    try:
        cdelt1 = abs(head["CDELT1"])
        # cdelt1 is in degree/px, convert it in arcsec/px
        return float(cdelt1)*3600.0
    except KeyError:
        # nova's wcs files have the scale in the comments
        comments = head['COMMENT']
        size = (len(comments))
        scale = -1
        for i in range(0, size):
            if comments[i][0:5] == 'scale':
                tkns = comments[i].split(' ')
                return float(tkns[1])


def dec_frm_header(head):
    '''
    look in header for width and height of image
   '''
    # nova's and AstroArt's wcs files have CRVAL2
    dec = head['CRVAL2']
    return dec


def wid_hei_frm_header(head):
    '''
    look in header for width and height of image
   '''
    try:
        # nova's wcs files have IMAGEW / IMAGEH
        width = head['IMAGEW']
        height = head['IMAGEH']
        return width, height
    except KeyError:
        try:
            # AstroArt's fits files have NAXIS1 / NAXIS2
            width = head['NAXIS1']
            height = head['NAXIS2']
            return width, height
        except KeyError:
            return 0, 0

def decdeg2dms(dd):
    mnt,sec = divmod(dd*3600,60)
    deg,mnt = divmod(mnt,60)
    return deg,mnt,sec

def cross(crd, img, colour):
    '''
    Annotate with a cross for the RA axis
    '''
    draw = ImageDraw.Draw(img)
    coords = crd[0]
    ax1 = coords[0]
    ay1 = coords[1]
    draw.line((ax1 - 30, ay1 - 30) + (ax1 + 30, ay1 + 30),
              fill=colour, width=3)
    draw.line((ax1 + 30, ay1 - 30) + (ax1 - 30, ay1 + 30),
              fill=colour, width=3)


def circle(centre, img, colour, name):
    '''
    Annotate with a circle
    '''
    font = ImageFont.load('symb24.pil')
    draw = ImageDraw.Draw(img)
    cen = centre[0]
    ax1 = cen[0]
    ay1 = cen[1]
    draw.ellipse((ax1 - 20, ay1 - 20, ax1 + 20, ay1 + 20),
                 fill=None, outline=colour)
    draw.text((ax1 + 30, ay1), name, fill=colour, font=font)


def cpcircle(centre, img, scl):
    '''
    Annotate with target circles
    '''
    font = ImageFont.load('helvR24.pil')
    draw = ImageDraw.Draw(img)
    cen = centre[0]
    ax1 = cen[0]
    ay1 = cen[1]
    number = [5, 10, 20, 40]
    for i in number:
        rad = (i*60)/scl
        draw.ellipse((ax1 - rad, ay1 - rad, ax1 + rad, ay1 + rad),
                     fill=None, outline='Green')
        draw.text((ax1 + (rad*26)/36, ay1 + (rad*26/36)), str(i),
                  font=font)
    draw.line((ax1 - 30, ay1) + (ax1 - 4, ay1), fill='Green', width=2)
    draw.line((ax1 +4, ay1) + (ax1 + 30, ay1), fill='Green', width=2)
    draw.line((ax1, ay1 - 30) + (ax1, ay1 - 4),fill='Green', width=2)
    draw.line((ax1, ay1 + 4) + (ax1, ay1 + 30),fill='Green', width=2)


class PhotoPolarAlign(tkinter.Frame):
    '''
    Our application as a class
    '''

    def dwarf_bar(self, txt, txt_process = "", txt_info = ""):
        '''
        Update the Dwarf bar
        '''
        self.dwarf_status_msg = txt
        self.dwarf_status_msg_process = txt_process
        self.dwarf_status_msg_info = txt_info
        msg_display = self.dwarf_status_msg
        if self.dwarf_status_msg_process:
            msg_display += " | " + self.dwarf_status_msg_process
        if self.dwarf_status_msg_info:
            msg_display += " : " + self.dwarf_status_msg_info
        self.wdwco.config(text=msg_display)
        self.wdwco.update()

    def stat_bar(self, txt):
        '''
        Update the Status bar
        '''
        self.stat_msg = txt
        self.wstat.config(text=self.stat_msg)
        self.wstat.update()
    
    def write_config_file(self):
        '''
        Update the user preferences file
        '''
        # the API key
        if not self.config.has_section('nova'):
            self.config.add_section('nova')
        self.config.set('nova', 'apikey', str(self.apikey.get()))
        # the image directory
        if not self.config.has_section('file'):
            self.config.add_section('file')
        self.config.set('file', 'imgdir', str(self.imgdir))
        # the geometry
        if not self.config.has_section('appearance'):
            self.config.add_section('appearance')
        self.config.set('appearance', 'geometry',
                        str(self.myparent.winfo_geometry()))
        # the operating options
        if not self.config.has_section('operations'):
            self.config.add_section('operations')
        self.config.set('operations','restrict scale',
                        str(self.restrict_scale.get()))
        # the local solve options
        if not self.config.has_section('local'):
            self.config.add_section('local')
        # self.config.set('local','shell',
        #                 str(self.local_shell.get()))
        # self.config.set('local','downscale',
        #                 str(self.local_downscale.get()))
        # self.config.set('local','configfile',
        #                 str(self.local_configfile.get()))
        # self.config.set('local','scale_units',
        #                 str(self.local_scale_units.get()))
        # self.config.set('local','scale_low',
        #                 str(self.local_scale_low.get()))
        # self.config.set('local','scale_hi',
        #                 str(self.local_scale_hi.get()))
        # self.config.set('local','xtra',
        #                 str(self.local_xtra.get()))
        #
        self.config.set("local","astap_path",
                        str(self.astap_path.get()))
        with open(self.cfgfn, 'w') as cfgfile:
            self.config.write(cfgfile)
        cfgfile.close()

    def settings_destroy(self):
        '''
        User asked to close the Settings
        '''
        self.write_config_file()
        self.wvar4.configure(text=('%.3s...........' % self.apikey.get()))
        self.settings_win.destroy()
                                
        
    def settings_open(self):
        '''
        Our Settings window
        '''
        # create child window
        win = tkinter.Toplevel()
        self.settings_win = win
        win.geometry('480x600')
        win.title('Settings')
        # get the API key information
        frm = tkinter.LabelFrame(win, borderwidth=2, relief='ridge', text='nova.astrometry.net')
        frm.pack(side='top', ipadx=20, padx=20, fill='x')
        nxt = tkinter.Label(frm, text='API Key')
        nxt.grid(row=0, column=0, pady=4, sticky='w')
        nxt =tkinter. Entry(frm, textvariable=self.apikey)
        nxt.grid(row=0, column=1, pady=4)
        nxt = tkinter.Label(frm, text='Restrict scale')
        nxt.grid(row=1, column=0, pady=4, sticky='w')
        nxt = tkinter.Checkbutton(frm, var=self.restrict_scale)
        nxt.grid(row=1, column=1, pady=4)

        frm = tkinter.LabelFrame(win, borderwidth=2, relief='ridge', text='Local solver Configuration')
        frm.pack(side='top', ipadx=20, padx=20, fill='x')
        
        nxt = tkinter.Label(frm, text='ASTAP Path')
        nxt.grid(row=0, column=0, pady=4, sticky='w')
        nxt = tkinter.Entry(frm, textvariable=self.astap_path)
        nxt.grid(row=0, column=1, pady=4)

        
        tkinter.Button(win, text='OK', command=self.settings_destroy).pack(pady=4)

    def copy_file_in_current_directory(self, source_filename, destination_filename):
        """
        Copies a file from source_filename to destination_filename in the current directory.
        :param source_filename: The name of the source file to copy.
        :param destination_filename: The name where the file should be copied to.
        """
        try:
            current_directory = os.getcwd()
            source_path = os.path.join(current_directory, source_filename)
            destination_path = os.path.join(current_directory, destination_filename)
            shutil.copy(source_path, destination_path)
            print(f"File copied from {source_filename} to {destination_filename}")
            return True
        except FileNotFoundError:
            print(f"Source file {source_filename} not found in the current directory.")
        except PermissionError:
            print(f"Permission denied. Unable to copy to {destination_filename}.")
        except Exception as e:
            print(f"Error copying file: {e}")

        return False

    def get_file_modification_time(self, file_path):
        return os.path.getmtime(file_path)

    def quit_method(self):
        '''
        User wants to quit
        '''
        # Signal the threads to stop
        print ("Quit program")

        self.write_config_file()
        self.myparent.destroy()

    def happy_with(self, wcs, img):
        '''
        check that .wcs (wcs) is compatible with .jpg (img)
        '''
        if os.path.exists(wcs):
            # DBG print wcs, 'exists'
            # check timestamps
            # DBG print os.stat(wcs).st_atime, os.stat(wcs).st_mtime, os.stat(wcs).st_ctime, 'wcs'
            # DBG print os.stat(img).st_atime, os.stat(img).st_mtime, os.stat(img).st_ctime, 'img'
            if os.stat(wcs).st_mtime> os.stat(img).st_mtime:
                return True
        return False

    def get_file(self, hint):
        '''
        Check if Dwarf is connected to take a photo
        '''

        img = False
        if self.dwarf_status:
            self.seestar_take_photo()
            time.sleep(1)
            self.dwarf_status_msg_process = "Downloading photo"
            self.dwarf_status_msg_info = "In progress"
            self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
            img = self.seestar_download_photo()

            if img:
                print ("Image Saved: ", img) 
                self.dwarf_status_msg_process = "Download photo"
                self.dwarf_status_msg_info = "Success"
                self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
            else:
                self.dwarf_status_msg_process = "Download photo"
                self.dwarf_status_msg_info = "Error"
                self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
                return False 
        else:
            '''
            User wants to select an image file
            '''
            options = {}
            options['filetypes'] = [('FIT files', '.fit .fits .FIT .FITS'),
                                    ('all files', '.*')]
            options['initialdir'] = self.imgdir
            titles = {}
            titles['v'] = 'The vertical image of the Celestial Pole region'
            titles['h'] = 'The horizontal image of the Celestial Pole region'
            titles['i'] = 'The horizontal image after Alt/Az adjustment'
            options['title'] = titles[hint]
            img = tkinter.filedialog.askopenfilename(**options)

        if img:
            print("img: ", img)
            wcs = os.path.splitext(img)[0] + '.wcs'
            if self.happy_with(wcs, img):
                self.update_solved_labels(hint, 'active')
            else:
                self.update_solved_labels(hint, 'disabled')
            self.imgdir = os.path.dirname(img)
            if hint == 'v':
                self.vimg_fn = img
                self.vwcs_fn = wcs
                self.havev = True
                self.wvar1.configure(text=os.path.basename(img))
                self.wvfn.configure(bg='green', activebackground='green')
            elif hint == 'h':
                self.himg_fn = img
                self.hwcs_fn = wcs
                self.haveh = True
                self.wvar2.configure(text=os.path.basename(img))
                self.whfn.configure(bg='green', activebackground='green')
            elif hint == 'i':
                self.iimg_fn = img
                self.iwcs_fn = wcs
                self.havei = True
                self.wvar3.configure(text=os.path.basename(img))
                self.wifn.configure(bg='green', activebackground='green')

    def update_scale(self, hint):
        try: 
            if hint == 'v':
                self.scale = scale_frm_wcs(self.vwcs_fn)
            elif hint == 'h':
                self.scale = scale_frm_wcs(self.hwcs_fn)
            elif hint == 'i':
                self.scale = scale_frm_wcs(self.iwcs_fn)
            self.havescale = True
            self.wvar5.configure(text=('%.2f' % self.scale))
        except:
            self.havescale = False
            self.wvar5.configure(text='--.--')
            return

    def solve(self, hint, solver):
        '''
        Solve an image
        '''
        if hint == 'h' or hint == 'v':
            if self.vimg_fn == self.himg_fn:
                self.stat_bar( ('Image filenames coincide - Check the Image ' +
                                'filenames'))
                return
        if hint == 'h':
            aimg = self.himg_fn
            awcs = self.hwcs_fn
        if hint == 'v':
            aimg = self.vimg_fn
            awcs = self.vwcs_fn
        if hint == 'i':
            aimg = self.iimg_fn
            awcs = self.iwcs_fn
        try:
            open(aimg)
        except IOError:
            self.stat_bar( ("couldn't open the image - Check the Image " +
                            'filename' + aimg))
            return
        self.stat_bar( 'Solving image...')
        if solver=='nova':
            self.img2wcs(self.apikey.get(), aimg, awcs, hint)
        if solver=='local':
            self.limg2wcs(aimg, awcs, hint)
        self.update_scale(hint)
        
    def update_display(self, cpcrd, the_scale):
        '''
        update Computed displayed quantities
        '''
        axis = self.axis
        x1a = axis[0]
        y1a = axis[1]
        x2a = cpcrd[0][0]
        y2a = cpcrd[0][1]
        self.scale = the_scale
        self.havescale = True
        self.wvar5.configure(text=('%.2f' % the_scale))
        self.wvar6.configure(text=str(int(x1a))+','+str(int(y1a)))
        self.wvar7.configure(text=(str(int(x2a)) +',' + str(int(y2a))))
        err = the_scale*numpy.sqrt((x1a-x2a)**2 + (y1a-y2a)**2)/60.0
        self.wvar8.configure(text=('%.2f' % err))
        if x2a > x1a:
            inst = 'Right '
        else:
            inst = 'Left '
        ddeg = abs(x2a - x1a)*the_scale/3600.0
        inst = inst + ('%02d:%02d:%02d' % decdeg2dms(ddeg))
        self.wvar9.configure(text=inst)
        if y2a > y1a:
            inst = inst + ' Down '
        else:
            inst = inst + ' Up '
        ddeg = abs(y2a - y1a)*the_scale/3600.0
        inst = inst + ('%02d:%02d:%02d' % decdeg2dms(ddeg))
        self.wvar9.configure(text=inst)

    def annotate_imp(self):
        '''
        Annotate the improvement image
        '''
        if self.iimg_fn == self.himg_fn:
            self.stat_bar( ('Image filenames coincide - Check the Image ' +
                            'filenames'))
            return
        try:
            imi = Image.open(self.iimg_fn)
            imh = Image.open(self.himg_fn)
            # Load the FITS hdulist using astropy.io.fits
            hdulisti = fits.open(self.iwcs_fn)
            hdulisth = fits.open(self.hwcs_fn)
        except IOError:
            return
        axis = self.axis
        try:
            axis[0]
        except:
            self.stat_bar("don't know where Polar Axis is - Find Polar Axis")
            return
        self.stat_bar('Annotating...')
        headi = hdulisti[0].header
        headh = hdulisth[0].header
        wcsi = wcs.WCS(headi)
        now = Time.now()
        if self.hemi == 'N':
            cp = SkyCoord(ra=0, dec=90, frame='fk5', unit='deg', equinox=now)
        else:
            cp = SkyCoord(ra=0, dec=-90, frame='fk5', unit='deg', equinox=now)
        cpj2000 = cp.transform_to(FK5(equinox='J2000'))
        cpskycrd = numpy.array([[cpj2000.ra.deg, cpj2000.dec.deg]],
                               numpy.float_)
        cpcrdi = wcsi.wcs_world2pix(cpskycrd, 1)
        scalei = scale_frm_header(headi)
        widthi = imi.width
        heighti = imi.height
        widthh = imh.width
        heighth = imh.height

        #widthi, heighti = wid_hei_frm_header(headi)
        #if wid_hei_frm_header(headi) != wid_hei_frm_header(headh) :
        if widthi != widthh or heighti != heighth:
            self.stat_bar('Incompatible image dimensions...')
            return
        if parity_frm_header(headi) == 0 :
            self.stat_bar('Wrong parity...')
            return
        self.update_display(cpcrdi, scalei)
        cpcircle(cpcrdi, imi, scalei)
        cross([axis], imi, 'Red')
        if self.hemi == 'N':
            poli = wcsi.wcs_world2pix(self.polaris, 1)
            lami = wcsi.wcs_world2pix(self.lam, 1)
            circle(poli, imi, 'White', 'a')
            circle(lami, imi, 'Orange', 'l')
            left = int(min(cpcrdi[0][0], poli[0][0], lami[0][0], axis[0]))
            right = int(max(cpcrdi[0][0], poli[0][0], lami[0][0], axis[0]))
            bottom = int(min(cpcrdi[0][1], poli[0][1], lami[0][1], axis[1]))
            top = int(max(cpcrdi[0][1], poli[0][1], lami[0][1], axis[1]))
        else:
            ori = wcsi.wcs_world2pix(self.chi, 1)
            whi = wcsi.wcs_world2pix(self.sigma, 1)
            rei = wcsi.wcs_world2pix(self.red, 1)
            circle(whi, imi, 'White', 's')
            circle(ori, imi, 'Orange', 'c')
            circle(rei, imi, 'Red', '!')
            left = int(min(cpcrdi[0][0], ori[0][0], whi[0][0], axis[0]))
            right = int(max(cpcrdi[0][0], ori[0][0], whi[0][0], axis[0]))
            bottom = int(min(cpcrdi[0][1], ori[0][1], whi[0][1], axis[1]))
            top = int(max(cpcrdi[0][1], ori[0][1], whi[0][1], axis[1]))
        margin = int(2500/scalei)
        xl = max(1, left - margin)
        xr = min(widthi, right + margin)
        yt = min(heighti, top + margin)
        yb = max(1, bottom - margin)
        croppedi = imi.crop((xl, yb, xr, yt))
        croppedi.load()
        crop_fn = os.path.splitext(self.iimg_fn)[0] + '_cropi.ppm'
        croppedi.save(crop_fn, 'PPM')
        self.create_imgwin(crop_fn, self.iimg_fn)
        self.stat_bar( 'Idle')

    def annotate(self):
        '''
        Find RA axis and Annotate the pair of horiz/vertical images
        '''
        if self.vimg_fn == self.himg_fn:
            self.stat_bar( ('Image filenames coincide - Check the Image ' +
                            'filenames'))
            return
        try:
            imh = Image.open(self.himg_fn)
            imv = Image.open(self.vimg_fn)
            # Load the FITS hdulist using astropy.io.fits
            hdulistv = fits.open(self.vwcs_fn)
            hdulisth = fits.open(self.hwcs_fn)
        except IOError:
            return
        self.stat_bar( 'Finding RA axis...')
        # Parse the WCS keywords in the primary HDU
        headv = hdulistv[0].header
        headh = hdulisth[0].header
        wcsv = wcs.WCS(headv)
        wcsh = wcs.WCS(headh)
        decv = dec_frm_header(headv)
        dech = dec_frm_header(headh)
        if decv > 65 and dech > 65:
            self.hemi = 'N'
        elif decv < -65 and dech < -65:
            self.hemi = 'S'
        else:
            self.stat_bar( 'Nowhere near (>25 deg) the Poles!')
            return
        now = Time.now()
        if self.hemi == 'N':
            cp = SkyCoord(ra=0, dec=90, frame='fk5', unit='deg', equinox=now)
        else:
            cp = SkyCoord(ra=0, dec=-90, frame='fk5', unit='deg', equinox=now)

        # CP now, in J2000 coordinates, precess
        cpj2000 = cp.transform_to(FK5(equinox='J2000'))
        # sky coordinates
        cpskycrd = numpy.array([[cpj2000.ra.deg, cpj2000.dec.deg]],
                               numpy.float_)
        # pixel coordinates
        cpcrdh = wcsh.wcs_world2pix(cpskycrd, 1)
        if self.hemi == 'N':
            print('Northern Celestial Pole', dech)
        else:
            print('Southern Celestial Pole', dech)
        scaleh = scale_frm_header(headh)
        widthh = imh.width
        heighth = imh.height
        widthv = imv.width
        heightv = imv.height
        #widthh, heighth = wid_hei_frm_header(headh)
        #if wid_hei_frm_header(headh) != wid_hei_frm_header(headv):
        if widthh != widthv or heighth != heightv:
            self.stat_bar( 'Incompatible image dimensions...')
            return
        if parity_frm_header(headh) == 0 or parity_frm_header(headv) == 0 :
            self.stat_bar( 'Wrong parity...')
            return
        
        def displacement(coords):
            '''
            the movement of a sky object in the two images
            '''
            pixcrd1 = numpy.array([coords], numpy.float_)
            skycrd = wcsv.wcs_pix2world(pixcrd1, 1)
            pixcrd2 = wcsh.wcs_world2pix(skycrd, 1)
            return pixcrd2 - pixcrd1
        axis = scipy.optimize.broyden1(displacement, [widthh/2, heighth/2])
        self.axis = axis
        self.update_display(cpcrdh, scaleh)
        #
        self.stat_bar( 'Annotating...')
        cpcircle(cpcrdh, imh, scaleh)
        cross([axis], imh, 'Red')
        # add reference stars
        if self.hemi == 'N':
            polh = wcsh.wcs_world2pix(self.polaris, 1)
            lamh = wcsh.wcs_world2pix(self.lam, 1)
            circle(polh, imh, 'White', 'a')
            circle(lamh, imh, 'Orange', 'l')
            left = int(min(cpcrdh[0][0], polh[0][0], lamh[0][0], axis[0]))
            right = int(max(cpcrdh[0][0], polh[0][0], lamh[0][0], axis[0]))
            bottom = int(min(cpcrdh[0][1], polh[0][1], lamh[0][1], axis[1]))
            top = int(max(cpcrdh[0][1], polh[0][1], lamh[0][1], axis[1]))
        else:
            orh = wcsh.wcs_world2pix(self.chi, 1)
            whh = wcsh.wcs_world2pix(self.sigma, 1)
            reh = wcsh.wcs_world2pix(self.red, 1)
            circle(whh, imh, 'White', 's')
            circle(orh, imh, 'Orange', 'c')
            circle(reh, imh, 'Red', '!')
            left = int(min(cpcrdh[0][0], orh[0][0], whh[0][0], axis[0]))
            right = int(max(cpcrdh[0][0], orh[0][0], whh[0][0], axis[0]))
            bottom = int(min(cpcrdh[0][1], orh[0][1], whh[0][1], axis[1]))
            top = int(max(cpcrdh[0][1], orh[0][1], whh[0][1], axis[1]))
        margin = int(2500/scaleh)
        xl = max(1, left - margin)
        xr = min(widthh, right + margin)
        yt = min(heighth, top + margin)
        yb = max(1, bottom - margin)
        croppedh = imh.crop((xl, yb, xr, yt))
        croppedh.load()
        crop_fn = os.path.splitext(self.himg_fn)[0] + '_croph.ppm'
        croppedh.save(crop_fn, 'PPM')
        self.create_imgwin(crop_fn, self.himg_fn)
        self.stat_bar( 'Idle')

    def create_imgwin(self, img_fn, title):
        '''
        creates a window to display an image
        '''
        # create child window
        img = tkinter.PhotoImage(file=img_fn)
        win = tkinter.Toplevel()
        wwid = min(800, img.width())
        whei = min(800, img.height())
        win.geometry(('%dx%d' % (wwid+28, whei+28)))
        win.title(os.path.basename(title))
        frame = tkinter.Frame(win, bd=0)
        frame.pack()
        xscrollbar = tkinter.Scrollbar(frame, orient='horizontal')
        xscrollbar.pack(side='bottom', fill='x')
        yscrollbar = tkinter.Scrollbar(frame, orient='vertical')
        yscrollbar.pack(side='right', fill='y')
        canvas = tkinter.Canvas(frame, bd=0, width=wwid, height=whei,
                        scrollregion=(0, 0, img.width(), img.height()),
                        xscrollcommand=xscrollbar.set,
                        yscrollcommand=yscrollbar.set)
        canvas.pack(side='top', fill='both', expand=1)
        canvas.create_image(0, 0, image=img, anchor='nw')
        xscrollbar.config(command=canvas.xview)
        yscrollbar.config(command=canvas.yview)
        frame.pack()
        # next statement is important! creates reference to img
        canvas.img = img

    def update_solved_labels(self, hint, sta):
        '''
        updates displayed Solved labels
        '''
        if hint == 'v':
            widget = self.wvok
        elif hint == 'h':
            widget = self.whok
        elif hint == 'i':
            widget = self.wiok
        # oldstate = widget.config()['state'][4]
        if (sta == 'active'):
            widget.configure(state='active', bg='green',
                             activebackground='green',
                             highlightbackground='green')
        elif (sta == 'disabled'):
            widget.configure(state='disabled', bg='red',
                             activebackground='red',
                             highlightbackground='red')
        widget.update()

    def seestar_download_photo(self):
        print("seestar_download_photo")
        url = seestar_api.get_last_image_url().partition(".jpg")[0] + ".fit"
        img_data = requests.get(url).content
        img_name = url.split("/")[-1]
        if not os.path.exists("Images"):
            os.mkdir("Images")
        with open("Images/"+img_name, "wb") as handler:
            handler.write(img_data)
        # delete images from seestar memory
        seestar_api.delete_image("Seestar_polar_align-sub")
        seestar_api.delete_image("Seestar_polar_align")
        return "Images/"+img_name

    def seestar_take_photo(self):
        print("seestar_take_photo")
        self.dwarf_status_msg_process = "Starting Seestar Take Photo..."
        self.dwarf_status_msg_info = ""
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
        setting = seestar_api.get_setting()
        exposure_stack = setting["Value"]["result"]["exp_ms"]["stack_l"]
        horizontal_calibration_enable = setting["Value"]["result"]["auto_3ppa_calib"]
        seestar_api.set_exposure(5000)
        seestar_api.set_horizontal_calibration(False)
        seestar_api.set_lp_filter(False)
        seestar_api.set_target_name("Seestar_polar_align")
        time.sleep(1)
        
        seestar_api.start_stack()
        self.dwarf_status_msg_process = "Take Photo..."
        self.dwarf_status_msg_info = "In Progress"
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
        time.sleep(2)
        if seestar_api.get_wheel_position() == 0:
            self.dwarf_status_msg_process = "Seestar is taking dark frames..."
            self.dwarf_status_msg_info = "Waiting"
            self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
            # Seestar is doing dark frames
            # wait
            while seestar_api.get_wheel_position() == 0:
                time.sleep(2)
            
            self.dwarf_status_msg_process = "Take Photo..."
            self.dwarf_status_msg_info = "In Progress"
            self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
            time.sleep(3)
        time.sleep(8)
        seestar_api.stop_stack()
        time.sleep(1)
        seestar_api.set_exposure(exposure_stack)
        seestar_api.set_horizontal_calibration(horizontal_calibration_enable)
        self.dwarf_status_msg_process = "Take Photo..."
        self.dwarf_status_msg_info = "Success"
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)

    def seestar_start_up(self):
        print("seestar_start_up")
        self.dwarf_status_msg_process = "Starting Seestar Start-up Sequence..."
        self.dwarf_status_msg_info = ""
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
        seestar_api.start_up_sequence()
        self.dwarf_status_msg_process = "Seestar Start-up Sequence"
        self.dwarf_status_msg_info = "Success"
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)

    def seestar_autofocus(self):
        print("seestar_autofocus")
        self.dwarf_status_msg_process = "Starting Seestar Autofocus..."
        self.dwarf_status_msg_info = ""
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
        seestar_api.auto_focus()
        self.dwarf_status_msg_process = "Seestar Autofocus"
        self.dwarf_status_msg_info = "Success"
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)

    def seestar_stop_autofocus(self):
        print("seestar_stop_autofocus")
        self.dwarf_status_msg_process = "Stopping Seestar Autofocus..."
        self.dwarf_status_msg_info = ""
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
        seestar_api.stop_auto_focus()
        self.dwarf_status_msg_process = "Seestar Stop Autofocus"
        self.dwarf_status_msg_info = "Success"
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)

    def seestar_move_ra_90(self, right):
        print("seestar_move_ra_90")
        self.dwarf_status_msg_process = "Moving RA Axis by 90 degrees..."
        self.dwarf_status_msg_info = ""
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)
        seestar_api.move_ra_90(right)
        self.dwarf_status_msg_process = "RA Axis moved"
        self.dwarf_status_msg_info = "Success"
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)

    def seestar_connect(self):
        try:
            if seestar_api.is_seestar_connected():
                self.dwarf_status_msg = "Connected"
                self.dwarf_status = True
            else:
                self.dwarf_status_msg = "Not connected"
                self.dwarf_status = False
        except requests.exceptions.ConnectionError:
            self.dwarf_status_msg = "Not connected"
            self.dwarf_status = False
        self.dwarf_bar(self.dwarf_status_msg, self.dwarf_status_msg_process, self.dwarf_status_msg_info)

    def create_widgets(self, master=None):
        '''
        creates the main window components
        '''
        self.myparent = master
        self.myparent.title('Photo Polar Alignment')
        #
        self.menubar = tkinter.Menu(master)
        self.filemenu = tkinter.Menu(self.menubar, tearoff=0)
        self.seestarmenu = tkinter.Menu(self.menubar, tearoff=0)
        self.helpmenu = tkinter.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label='File', menu=self.filemenu)
        self.menubar.add_cascade(label='Seestar', menu=self.seestarmenu)
        self.menubar.add_cascade(label='Help', menu=self.helpmenu)
        self.filemenu.add_command(label='Settings...',
                                  command=self.settings_open)
        self.filemenu.add_command(label='Exit', command=self.quit_method)
        self.seestarmenu.add_command(label='Connect...',
                                  command=lambda:threading.Thread(target=self.seestar_connect).start())
        self.seestarmenu.add_command(label='Move RA 90° Right...',
                                  command=lambda:threading.Thread(target=lambda:self.seestar_move_ra_90(True)).start())
        self.seestarmenu.add_checkbutton(label='Move RA 90° Left...',
                                         command=lambda:threading.Thread(target=lambda:self.seestar_move_ra_90(False)).start())
        self.seestarmenu.add_command(label='Start-up Sequence...',
                                  command=self.seestar_start_up)
        self.autofocusmenu = tkinter.Menu(self.seestarmenu, tearoff=0)
        self.seestarmenu.add_cascade(label='Autofocus', menu=self.autofocusmenu)

        self.autofocusmenu.add_command(label='Autofocus...',
                                  command=self.seestar_autofocus)
        self.autofocusmenu.add_command(label='Stop Autofocus...',
                                  command=self.seestar_stop_autofocus)

        self.helpmenu.add_command(label='Help', command=help_f)
        self.helpmenu.add_command(label='About...', command=about_f)
        self.myparent.config(menu=self.menubar)
        # #################################################################
        self.wfrop = tkinter.LabelFrame(master, text='Operations')
        self.wfrop.pack(side='top', fill='x')
        #
        nxt = tkinter.Button(self.wfrop, image=self.vicon, command=lambda : threading.Thread(target=lambda: self.get_file('v')).start())
        nxt.grid(row=0, column=0, sticky='ew', padx=10, pady=4, rowspan=3)
        self.wvfn = nxt
        nxt = tkinter.Button(self.wfrop, text='Nova', command=lambda : self.solve('v','nova'))
        nxt.grid(row=0, column=1, sticky='ew', padx=10, pady=4)
        self.wvsol = nxt
        nxt = tkinter.Button(self.wfrop, text='Local', command=lambda : self.solve('v','local'))
        nxt.grid(row=1, column=1, sticky='ew', padx=10, pady=4)
        self.wlvsol = nxt
        nxt = tkinter.Label(self.wfrop, text='Solved', state='disabled')
        nxt.grid(row=2, column=1, sticky='ew', padx=10, pady=4)
        self.wvok = nxt
        #
        nxt = tkinter.Button(self.wfrop, image=self.hicon, command=lambda : threading.Thread(target=lambda: self.get_file('h')).start())
        nxt.grid(row=3, column=0, sticky='ew', padx=10, pady=4, rowspan=3)
        self.whfn = nxt
        nxt = tkinter.Button(self.wfrop, text='Nova', command=lambda : self.solve('h','nova'))
        nxt.grid(row=3, column=1, sticky='ew', padx=10, pady=4)
        self.whsol = nxt
        nxt = tkinter.Button(self.wfrop, text='Local', command=lambda : self.solve('h','local'))
        nxt.grid(row=4, column=1, sticky='ew', padx=10, pady=4)
        self.wlhsol = nxt
        nxt = tkinter.Label(self.wfrop, text='Solved', state='disabled')
        nxt.grid(row=5, column=1, sticky='ew', padx=10, pady=4)
        self.whok = nxt
        #
        nxt = tkinter.Button(self.wfrop, text='Find Polar Axis',
                     command=self.annotate)
        nxt.grid(row=6, column=0, sticky='ew', padx=10, pady=4, columnspan=2)
        self.wann = nxt
        #
        nxt = tkinter.Button(self.wfrop, image=self.iicon, command=lambda : threading.Thread(target=lambda: self.get_file('i')).start())
        nxt.grid(row=3, column=3, sticky='ew', padx=10, pady=4, rowspan=3)
        self.wifn = nxt
        nxt = tkinter.Button(self.wfrop, text='Nova', command=lambda : self.solve('i','nova'))
        nxt.grid(row=3, column=4, sticky='ew', padx=10, pady=4)
        self.wisol = nxt
        nxt = tkinter.Button(self.wfrop, text='Local', command=lambda : self.solve('i','local'))
        nxt.grid(row=4, column=4, sticky='ew', padx=10, pady=4)
        self.wlisol = nxt
        nxt = tkinter.Label(self.wfrop, text='Solved', state='disabled')
        nxt.grid(row=5, column=4, sticky='ew', padx=10, pady=4)
        self.wiok = nxt
        #
        nxt = tkinter.Button(self.wfrop, text='Show Improvement',
                     command=self.annotate_imp)
        nxt.grid(row=6, column=3, sticky='ew', padx=10, pady=4, columnspan=2)
        self.wanni = nxt
        # #################################################################

        nxt = tkinter.LabelFrame(master, borderwidth=2, relief='ridge',
                         text='Info')
        nxt.pack(side='top', fill='x')
        self.wfrvar = nxt
        nxt = tkinter.Label(self.wfrvar, text = 'Given')
        nxt.grid(row=0, column=1, columnspan=2, sticky='w')
        nxt = tkinter.Label(self.wfrvar, anchor='w', text='Vertical:')
        nxt.grid(row=1, column=0, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text='---------')
        nxt.grid(row=1, column=1, sticky='e')
        self.wvar1 = nxt
        nxt = tkinter.Label(self.wfrvar, text='Horizontal:')
        nxt.grid(row=2, column=0, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text='---------')
        nxt.grid(row=2, column=1, sticky='e')
        self.wvar2 = nxt
        nxt = tkinter.Label(self.wfrvar, text='Improved:')
        nxt.grid(row=3, column=0, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text='---------')
        nxt.grid(row=3, column=1, sticky='e')
        self.wvar3 = nxt
        nxt = tkinter.Label(self.wfrvar, text='API key:')
        nxt.grid(row=4, column=0, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text=('%.3s...........' % self.apikey.get()))
        nxt.grid(row=4, column=1, sticky='e')
        self.wvar4 = nxt

        nxt = tkinter.Label(self.wfrvar, text = 'Computed')
        nxt.grid(row=0, column=3, columnspan=2, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text='Scale (arcsec/pixel):')
        nxt.grid(row=1, column=2, sticky='w')
        if self.havescale:
            nxt = tkinter.Label(self.wfrvar, text=self.scale)
        else:
            nxt = tkinter.Label(self.wfrvar, text='--.--')
        nxt.grid(row=1, column=3, sticky='e')
        self.wvar5 = nxt
        nxt = tkinter.Label(self.wfrvar, text='RA axis position:')
        nxt.grid(row=2, column=2, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text='---,---')
        nxt.grid(row=2, column=3, sticky='e')
        self.wvar6 = nxt
        nxt = tkinter.Label(self.wfrvar, text='CP position:')
        nxt.grid(row=3, column=2, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text='---,---')
        nxt.grid(row=3, column=3, sticky='e')
        self.wvar7 = nxt
        nxt = tkinter.Label(self.wfrvar, text='Error (arcmin):')
        nxt.grid(row=4, column=2, sticky='w')
        nxt = tkinter.Label(self.wfrvar, text='--.--')
        nxt.grid(row=4, column=3, sticky='e')
        self.wvar8 = nxt
        # #################################################################
        nxt = tkinter.LabelFrame(master, borderwidth=2, relief='ridge',
                         text='Move (dd:mm:ss)')
        nxt.pack(side='top', fill='x')
        self.wfrmo = nxt
        nxt = tkinter.Label(self.wfrmo, anchor='center', font='-weight bold -size 14')
        nxt.pack(anchor='center')
        self.wvar9 = nxt
        # #################################################################
        nxt = tkinter.LabelFrame(master, borderwidth=2, relief='ridge',
                         text='Seestar Status')
        nxt.pack(side='top', fill='x')
        self.wfdwco = nxt
        nxt = tkinter.Label(self.wfdwco, anchor='w', text=self.dwarf_status_msg)
        nxt.pack(anchor='w')
        self.wdwco = nxt

        # #################################################################
        nxt = tkinter.LabelFrame(master, borderwidth=2, relief='ridge', text='Status')
        nxt.pack(side='bottom', fill='x')
        self.wfrst = nxt
        nxt = tkinter.Label(self.wfrst, anchor='w', text=self.stat_msg)
        nxt.pack(anchor='w')
        self.wstat = nxt

    def __init__(self, master=None):
        
        self.master = master
        self.master.protocol("WM_DELETE_WINDOW", self.quit_method)

        # a F8Ib 2.0 mag star, Alpha Ursa Minoris
        self.polaris = numpy.array([[037.954561, 89.264109]], numpy.float_)
        #
        # a M1III 6.4 mag star, Lambda Ursa Minoris
        self.lam = numpy.array([[259.235229, 89.037706]], numpy.float_)
        #
        # a F0III 5.4 mag star, Sigma Octans
        self.sigma = numpy.array([[317.195164, -88.956499]], numpy.float_)
        #
        # a K3IIICN 5.3 mag star, Chi Octans
        self.chi = numpy.array([[283.696388, -87.605843]], numpy.float_)
        #
        # a M1III 7.2 mag star, HD90104
        self.red = numpy.array([[130.522862, -89.460536]], numpy.float_)
        #
        # the pixel coords of the RA axis, if solution exists
        self.axis = None
        self.havea = False
        # the Settings window
        self.settings_win = None
        # the User preferences file
        self.cfgfn = 'PPA.ini'

        self.astap_path = tkinter.StringVar()

        # self.local_shell = tkinter.StringVar()
        # self.local_downscale = tkinter.IntVar()
        # self.local_configfile = tkinter.StringVar()
        # self.local_scale_units = tkinter.StringVar()
        # self.local_scale_low = tkinter.DoubleVar()
        # self.local_scale_hi = tkinter.DoubleVar()
        # self.local_xtra = tkinter.StringVar()
        

        # Read the User preferences
        self.config = configparser.ConfigParser()
        self.config.read(self.cfgfn)
        # ...the key
        try:
            k_ini = self.config.get('nova', 'apikey')
        except :
            k_ini = None
        print ("k_ini", k_ini)
        self.apikey = tkinter.StringVar(value=k_ini)
        # ...the Image directory
        try: 
            self.imgdir = self.config.get('file', 'imgdir')
        except :        
            self.imgdir = None
        # ...geometry
        try:
            self.usergeo = self.config.get('appearance', 'geometry')
        except :
            self.usergeo = None
        master.geometry(self.usergeo)
        # do we want to help solves by restricting the scale once we have an estimate
        self.restrict_scale = tkinter.IntVar(value=0)
        try:
            self.restrict_scale.set(self.config.getint('operations','restrict scale'))
        except:
            self.restrict_scale.set(value=0)
            
        # the filenames of images
        self.vimg_fn = ''
        self.havev = False
        self.himg_fn = ''
        self.haveh = False
        self.iimg_fn = ''
        self.havei = False
        # the filenames of the .wcs solutions
        self.vwcs_fn = ''
        self.hwcs_fn = ''
        self.iwcs_fn = ''
        # the button icons
        self.vicon = tkinter.PhotoImage(file='v2_2.ppm')
        self.hicon = tkinter.PhotoImage(file='h2_2.ppm')
        self.iicon = tkinter.PhotoImage(file='i2_2.ppm')
        # the solved image scale
        self.havescale = False
        self.scale = None
        # the discovered hemisphere
        self.hemi = None
        # initialise attributes set elsewhere
        self.menubar = None
        self.helpmenu = None
        self.filemenu = None
        self.wfrop = None
        self.wvfn = None
        self.wvsol = None
        self.wlvsol = None
        self.wvok = None

        self.whfn = None
        self.whsol = None
        self.wlhsol = None
        self.whok = None

        self.wifn = None
        self.wisol = None
        self.wlisol = None
        self.wiok = None

        self.wann = None
        self.wanni = None

        self.wfr2 = None
        self.wfrvar = None
        self.wvar1 = None
        self.wvar2 = None
        self.wvar3 = None
        self.wvar4 = None
        self.wfrcomp = None
        self.wvar5 = None
        self.wvar6 = None
        self.wvar7 = None
        self.wvar8 = None

        self.wfrmo = None
        self.wvar9 = None

        self.wfdwco = None
        self.wdwco = None

        self.wfrst = None
        self.wstat = None

        self.myparent = None

        self.dwarf_status = False
        self.dwarf_status_msg = "Not connected"
        self.dwarf_status_msg_process = ""
        self.dwarf_status_msg_info = ""

        self.stat_msg = 'Idle'
        tkinter.Frame.__init__(self, master)
        self.create_widgets(master)
        # check local solver
        self.wlvsol.configure(state='disabled')
        self.wlhsol.configure(state='disabled')
        self.wlisol.configure(state='disabled')
        try:
            # self.local_shell.set(self.config.get('local','shell'))
            # self.local_downscale.set(self.config.getint('local','downscale'))
            # self.local_configfile.set(self.config.get('local','configfile'))
            # self.local_scale_units.set(self.config.get('local','scale_units'))
            # self.local_scale_low.set(self.config.getfloat('local','scale_low'))
            # self.local_scale_hi.set(self.config.getfloat('local','scale_hi'))
            # self.local_xtra.set(self.config.get('local','xtra'))

            self.astap_path.set(self.config.get("local","astap_path"))

            #exit_status = os.system(f"{str(self.astap_path.get())} -f > /dev/null")
            #if exit_status != 0:
            if not os.path.exists(self.astap_path.get()):
                print("Can't use local astap solver, check astap path in settings")
            else:
                self.wlvsol.configure(state='active')
                self.wlhsol.configure(state='active')
                self.wlisol.configure(state='active')
        except:
            # self.local_shell.set('')
            # self.local_downscale.set(1)
            # self.local_configfile.set('')
            # self.local_scale_units.set('')
            # self.local_scale_low.set(0)
            # self.local_scale_hi.set(0)
            # self.local_xtra.set('')

            self.astap_path.set("")
        if (not self.apikey.get() or self.apikey.get()=='') and (not self.astap_path or self.astap_path.get()==""):
            self.settings_open()
        self.pack()
        #

    def limg2wcs(self, filename, wcsfn, hint):
        t_start = time.time()
        #if (('OSTYPE' in os.environ and os.environ['OSTYPE']=='linux') or
        #    (os.uname()[0]=='Linux') or
        #    ('OSTYPE' in os.environ and os.environ['OSTYPE']=='darwin') or
        #    ('OS'     in os.environ and os.environ['OS']    =='Windows_NT')):

        print('___________________________________________________________')
        # cmd = 'solve-field -b ' + self.local_configfile.get()
        # if self.havescale and self.restrict_scale.get()==1:
        #     up_lim = self.scale*1.05
        #     lo_lim = self.scale*0.95
        #     cmd = cmd + (' -u app -L %.2f -H %.2f ' % (lo_lim,  up_lim))
        # else:
        #     cmd = cmd + ' -u ' + self.local_scale_units.get()
        #     cmd = cmd + (' -L %.2f' % self.local_scale_low.get())
        #     cmd = cmd + (' -H %.2f' % self.local_scale_hi.get())
        # if self.local_downscale.get() != 1:    
        #     cmd = cmd + (' -z %d' % self.local_downscale.get())
        # cmd = cmd + ' ' + self.local_xtra.get()
        # cmd = cmd + ' -O '
        # cmd = cmd + ' \\"%s\\"'
        # template = ((self.local_shell.get() % cmd))
        # # print template
        # cmd = (template % filename)
        subprocess.run([self.astap_path.get(),"-f",filename,"-r","30","-wcs","-speed","slow","-ra","0","-spd","180"])
        print('___________________________________________________________')
        if self.happy_with(wcsfn, filename):
            self.update_solved_labels(hint, 'active')
        else:
            # try south hemisphere
            subprocess.run([self.astap_path.get(),"-f",filename,"-r","30","-wcs","-speed","slow","-ra","0","-spd","0"])
            if self.happy_with(wcsfn, filename):
                self.update_solved_labels(hint, 'active')
            else:
                self.update_solved_labels(hint, 'disabled')
        self.update_scale(hint)
        self.stat_bar('Idle')
        print('local solve time ' + str(time.time()-t_start))
        print('___________________________________________________________')            

    
    def img2wcs(self, ankey, filename, wcsfn, hint):
        '''
        Plate solves one image
        '''
        t_start = time.time()
        parser = optparse.OptionParser()
        parser.add_option('--server', dest='server',
                        default=NovaClient.default_url,
                        help='Set server base URL (eg, %default)')
        parser.add_option('--apikey', '-k', dest='apikey',
                        help='API key for Astrometry.net web service; if not' +
                        'given will check AN_API_KEY environment variable')
        parser.add_option('--upload', '-u', dest='upload', help='Upload a file')
        parser.add_option('--wait', '-w', dest='wait', action='store_true',
                        help='After submitting, monitor job status')
        parser.add_option('--wcs', dest='wcs',
                        help='Download resulting wcs.fits file, saving to ' +
                        'given filename; implies --wait if --urlupload or' +
                        '--upload')
        parser.add_option('--kmz', dest='kmz',
                        help='Download resulting kmz file, saving to given ' +
                        'filename; implies --wait if --urlupload or --upload')
        parser.add_option('--urlupload', '-U', dest='upload_url',
                        help='Upload a file at specified url')
        parser.add_option('--scale-units', dest='scale_units',
                        choices=('arcsecperpix', 'arcminwidth', 'degwidth',
                                'focalmm'),
                        help='Units for scale estimate')
        parser.add_option('--scale-lower', dest='scale_lower', type=float,
                        help='Scale lower-bound')
        parser.add_option('--scale-upper', dest='scale_upper', type=float,
                        help='Scale upper-bound')
        parser.add_option('--scale-est', dest='scale_est', type=float,
                        help='Scale estimate')
        parser.add_option('--scale-err', dest='scale_err', type=float,
                        help='Scale estimate error (in PERCENT), eg "10" if' +
                        'you estimate can be off by 10%')
        parser.add_option('--ra', dest='center_ra', type=float, help='RA center')
        parser.add_option('--dec', dest='center_dec', type=float,
                        help='Dec center')
        parser.add_option('--radius', dest='radius', type=float,
                        help='Search radius around RA,Dec center')
        parser.add_option('--downsample', dest='downsample_factor', type=int,
                        help='Downsample image by this factor')
        parser.add_option('--parity', dest='parity', choices=('0', '1'),
                        help='Parity (flip) of image')
        parser.add_option('--tweak-order', dest='tweak_order', type=int,
                        help='SIP distortion order (default: 2)')
        parser.add_option('--crpix-center', dest='crpix_center',
                        action='store_true', default=None,
                        help='Set reference point to center of image?')
        parser.add_option('--sdss', dest='sdss_wcs', nargs=2,
                        help='Plot SDSS image for the given WCS file; write ' +
                        'plot to given PNG filename')
        parser.add_option('--galex', dest='galex_wcs', nargs=2,
                        help='Plot GALEX image for the given WCS file; write' +
                        'plot to given PNG filename')
        parser.add_option('--substatus', '-s', dest='sub_id',
                        help='Get status of a submission')
        parser.add_option('--jobstatus', '-j', dest='job_id',
                        help='Get status of a job')
        parser.add_option('--jobs', '-J', dest='myjobs', action='store_true',
                        help='Get all my jobs')
        parser.add_option('--jobsbyexacttag', '-T', dest='jobs_by_exact_tag',
                        help='Get a list of jobs associated with a given' +
                        'tag--exact match')
        parser.add_option('--jobsbytag', '-t', dest='jobs_by_tag',
                        help='Get a list of jobs associated with a given tag')
        parser.add_option('--private', '-p', dest='public', action='store_const',
                        const='n', default='y',
                        help='Hide this submission from other users')
        parser.add_option('--allow_mod_sa', '-m', dest='allow_mod',
                        action='store_const', const='sa', default='d',
                        help='Select license to allow derivative works of ' +
                        'submission, but only if shared under same conditions ' +
                        'of original license')
        parser.add_option('--no_mod', '-M', dest='allow_mod', action='store_const',
                        const='n', default='d',
                        help='Select license to disallow derivative works of ' +
                        'submission')
        parser.add_option('--no_commercial', '-c', dest='allow_commercial',
                        action='store_const', const='n', default='d',
                        help='Select license to disallow commercial use of' +
                        ' submission')
        # load opt with defaults, as above
        opt, args = parser.parse_args([''.split()])
        # add given arguments
        opt.wcs = wcsfn
        opt.apikey = ankey
        opt.upload = filename
        if self.havescale and self.restrict_scale.get() == 1:
            opt.scale_units = 'arcsecperpix'
            opt.scale_est = ('%.2f' % self.scale)
            opt.scale_err = 5
        # DEBUG print opt
        print('with estimated scale', opt.scale_est)
        args = {}
        args['apiurl'] = opt.server
        clnt = NovaClient(**args)
        try:
            clnt.login(opt.apikey)
        except RequestError as URLError:
            self.stat_bar(("Couldn't log on to nova.astrometry.net " +
                            '- Check the API key'))
            return
        if opt.upload or opt.upload_url:
            if opt.wcs or opt.kmz:
                opt.wait = True
            kwargs = dict()
            if opt.scale_lower and opt.scale_upper:
                kwargs.update(scale_lower=opt.scale_lower,
                            scale_upper=opt.scale_upper,
                            scale_type='ul')
            elif opt.scale_est and opt.scale_err:
                kwargs.update(scale_est=opt.scale_est,
                            scale_err=opt.scale_err,
                            scale_type='ev')
            elif opt.scale_lower or opt.scale_upper:
                kwargs.update(scale_type='ul')
                if opt.scale_lower:
                    kwargs.update(scale_lower=opt.scale_lower)
                if opt.scale_upper:
                    kwargs.update(scale_upper=opt.scale_upper)

            for key in ['scale_units', 'center_ra', 'center_dec', 'radius',
                        'downsample_factor', 'tweak_order', 'crpix_center', ]:
                if getattr(opt, key) is not None:
                    kwargs[key] = getattr(opt, key)
            if opt.parity is not None:
                kwargs.update(parity=int(opt.parity))
            if opt.upload:
                upres = clnt.upload(opt.upload, **kwargs)
            stat = upres['status']
            if stat != 'success':
                print('Upload failed: status', stat)
                print(upres)
                sys.exit(-1)
            opt.sub_id = upres['subid']
        if opt.wait:
            if opt.job_id is None:
                if opt.sub_id is None:
                    print("Can't --wait without a submission id or job id!")
                    sys.exit(-1)
                while True:
                    stat = clnt.sub_status(opt.sub_id, justdict=True)
                    # print 'Got status:', stat
                    jobs = stat.get('jobs', [])
                    if len(jobs):
                        for j in jobs:
                            if j is not None:
                                break
                        if j is not None:
                            print('Selecting job id', j)
                            opt.job_id = j
                            break
                    time.sleep(5)
            success = False
            while True:
                stat = clnt.job_status(opt.job_id, justdict=True)
                # print 'Got job status:', stat
                # TODO : stat may be None! should recover
                if stat.get('status', '') in ['success']:
                    success = (stat['status'] == 'success')
                    break
                time.sleep(5)
            if success:
                clnt.job_status(opt.job_id)
                retrieveurls = []
                if opt.wcs:
                    # We don't need the API for this, just construct URL
                    url = opt.server.replace('/api/', '/wcs_file/%i' % opt.job_id)
                    retrieveurls.append((url, opt.wcs))
                for url, fne in retrieveurls:
                    print('Retrieving file from', url)
                    fle = urlopen(url)
                    txt = fle.read()
                    wfl = open(fne, 'wb')
                    wfl.write(txt)
                    wfl.close()
                    print('Wrote to', fne)
                    self.update_solved_labels(hint, 'active')
                    self.stat_bar('Idle')
                    print('nova solve time ' + str(time.time()-t_start))
                    print('___________________________________________________________')            
            opt.job_id = None
            opt.sub_id = None
        if opt.sub_id:
            print(clnt.sub_status(opt.sub_id))
        if opt.job_id:
            print(clnt.job_status(opt.job_id))
        if opt.jobs_by_tag:
            tag = opt.jobs_by_tag
            print(clnt.jobs_by_tag(tag, None))
        if opt.jobs_by_exact_tag:
            tag = opt.jobs_by_exact_tag
            print(clnt.jobs_by_tag(tag, 'yes'))
        if opt.myjobs:
            jobs = clnt.myjobs()
            print(jobs)

ROOT = tkinter.Tk()
ROOT.geometry('600x500+300+300')
APP = PhotoPolarAlign(master=ROOT)
ROOT.mainloop()

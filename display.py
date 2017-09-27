import argparse
import sys
import os
import numpy as np
#### import the simple module from the paraview
from paraview.simple import *


description = "ParaView python script to generate 3D event displays from 3D spacepoint data created with pixy_roimux. Run with ParaView pvpython."
parser = argparse.ArgumentParser(description=description)
parser.add_argument("hits_filename", help="CSV file containing the 3D hits.")
parser.add_argument("-p", "--pca_filename", help="CSV file containing the PCA data. If not specified, no PCA line is drawn.")
parser.add_argument("-o", "--plot_filename", help="Save plot to file. If not specified, the interactive display is started.")
parser.add_argument("-c", "--colour_column", help="Specify the colour column. Q for charge and A for ambiguities (default: %(default)s).", default="Q")
args = parser.parse_args()

if args.plot_filename:
    if args.plot_filename[-6:] == ".webgl":
        create_plot = "view"
    elif args.plot_filename[-4:] == ".png":
        create_plot = "screenshot"
    else:
        print("ERROR: Unrecognised output file type: " + args.plot_filename + "!")
        sys.exit(1)
else:
    create_plot = ""

sampleTime = .21 # us
driftSpeed = .21 # cm/us
pixelPitch = .248 # cm
readoutCenterX = 18.5 * pixelPitch
readoutCenterY = 18.5 * pixelPitch
tpcLength = 58.6 # cm
tpcRadius = 5.05 # cm


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
hits_csv = CSVReader(FileName=[args.hits_filename])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=hits_csv)
tableToPoints1.XColumn = 'X'
tableToPoints1.YColumn = 'Y'
tableToPoints1.ZColumn = 'Z'

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1609, 834]

# set active view
SetActiveView(renderView1)

# create a new 'Glyph'
# sphere
glyph1 = Glyph(Input=tableToPoints1,
    GlyphType='Sphere')
glyph1.Scalars = ['POINTS', args.colour_column]
glyph1.Vectors = ['POINTS', 'None']
glyph1.ScaleMode = 'off'
glyph1.ScaleFactor = .3
glyph1.GlyphMode = 'All Points'
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphType.ThetaResolution = 8
glyph1.GlyphType.PhiResolution = 8

# create a new 'Glyph'
# box
#boxWidth = pixelPitch
#boxHeight = sampleTime * driftSpeed

#glyph1 = Glyph(Input=tableToPoints1,
#    GlyphType='Box')
#glyph1.Scalars = ['POINTS', args.colour_column]
#glyph1.Vectors = ['POINTS', 'None']
#glyph1.ScaleMode = 'off'
#glyph1.ScaleFactor = 1.0
#glyph1.GlyphMode = 'All Points'
#glyph1.GlyphTransform = 'Transform2'
#glyph1.GlyphType.XLength = boxWidth
#glyph1.GlyphType.YLength = boxWidth
#glyph1.GlyphType.ZLength = boxHeight
#glyph1.GlyphType.Center = [(boxWidth / 2.), (boxWidth / 2.), (boxHeight / 2.)]

# get color transfer function/color map for args.colour_column
qLUT = GetColorTransferFunction(args.colour_column)

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.ColorArrayName = ['POINTS', args.colour_column]
glyph1Display.LookupTable = qLUT
glyph1Display.OSPRayScaleArray = args.colour_column
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 59.83425195217133
glyph1Display.SelectScaleArray = args.colour_column
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GaussianRadius = 29.917125976085664
glyph1Display.SetScaleArray = ['POINTS', args.colour_column]
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', args.colour_column]
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.Orientation = [180.0, 0.0, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for args.colour_column
qPWF = GetOpacityTransferFunction(args.colour_column)

# set active source
SetActiveSource(hits_csv)

# create a new 'Cylinder'
# coordinates of the cylinder center
cylX = 0.
cylY = 0.
cylZ = 0.

cylinder1 = Cylinder()

# Properties modified on cylinder1
cylinder1.Resolution = 60
cylinder1.Height = tpcLength
cylinder1.Radius = tpcRadius
# Note the permutation due to the rotation of the cylinder (default is along Y axis).
cylinder1.Center = [cylX, cylZ, - cylY]

# show data in view
cylinder1Display = Show(cylinder1, renderView1)
# trace defaults for the display properties.
cylinder1Display.ColorArrayName = [None, '']
cylinder1Display.OSPRayScaleArray = 'Normals'
cylinder1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cylinder1Display.SelectOrientationVectors = 'None'
cylinder1Display.ScaleFactor = 60.0
cylinder1Display.SelectScaleArray = 'None'
cylinder1Display.GlyphType = 'Arrow'
cylinder1Display.GaussianRadius = 30.0
cylinder1Display.SetScaleArray = [None, '']
cylinder1Display.ScaleTransferFunction = 'PiecewiseFunction'
cylinder1Display.OpacityArray = [None, '']
cylinder1Display.OpacityTransferFunction = 'PiecewiseFunction'
cylinder1Display.Orientation = [90.0, 0.0, 0.0]
cylinder1Display.Opacity = 0.05
cylinder1Display.DiffuseColor = [0.0, 1.0, 0.0]

if args.pca_filename:
    pca = np.loadtxt(args.pca_filename, delimiter=",")
    l1 = (0 - tpcLength / 2. - pca[0][2]) / pca[1][2]
    l2 = (0 + tpcLength / 2. - pca[0][2]) / pca[1][2]
    line_start = pca[0] + l1 * pca[1]
    line_end = pca[0] + l2 * pca[1]
    
    # create a new 'Line'
    line1 = Line()
    
    # Properties modified on line1
    line1.Point1 = line_start
    line1.Point2 = line_end
    
    # show data in view
    line1Display = Show(line1, renderView1)
    # trace defaults for the display properties.
    line1Display.ColorArrayName = [None, '']
    line1Display.OSPRayScaleArray = 'Texture Coordinates'
    line1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    line1Display.SelectOrientationVectors = 'None'
    line1Display.ScaleFactor = 0.1
    line1Display.SelectScaleArray = 'None'
    line1Display.GlyphType = 'Arrow'
    line1Display.GaussianRadius = 0.05
    line1Display.SetScaleArray = [None, '']
    line1Display.ScaleTransferFunction = 'PiecewiseFunction'
    line1Display.OpacityArray = [None, '']
    line1Display.OpacityTransferFunction = 'PiecewiseFunction'
    line1Display.Orientation = [180.0, 0.0, 0.0]

# set active source
#SetActiveSource(glyph1)
renderView1.Background = [0.67, 0.67, 0.67]
Render()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-115.95241218116812, 6.2190223816782236, 0.]
renderView1.CameraFocalPoint = [4.660341811180115, 6.2190223816782236, 0.]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 312.16877815484816
renderView1.ViewSize = [800, 1000]

# reset view to fit data
#renderView1.ResetCamera()

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

if create_plot == "view":
    ExportView(args.plot_filename,
            view=renderView1)
elif create_plot == "screenshot":
    SaveScreenshot(args.plot_filename,
            magnification=2,
            quality=100,
            view=renderView1)
else:
    Interact(view=renderView1)

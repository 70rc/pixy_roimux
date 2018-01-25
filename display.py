import argparse
import json
import sys
import numpy as np
import vtk
#### import the simple module from the paraview
import paraview.simple as pv


description = "ParaView python script to generate 3D event displays from 3D spacepoint data created with pixy_roimux. Run with ParaView pvpython."
parser = argparse.ArgumentParser(description=description)
parser.add_argument("runParamsFile")
parser.add_argument("hitsFile", help="CSV file containing the 3D hits.")
parser.add_argument("-p", "--pcaFile", help="CSV file containing the PCA data. If not specified, no PCA line is drawn.")
parser.add_argument("-o", "--plotFile", help="Save plot to file. If not specified, the interactive display is started.")
parser.add_argument("-c", "--colourColumn", help="Specify the colour column. Q for charge and A for ambiguities (default: %(default)s).", default="Q")
args = parser.parse_args()

if args.plotFile:
    if args.plotFile[-6:] == ".webgl":
        createPlot = "view"
    elif args.plotFile[-4:] == ".png":
        createPlot = "screenshot"
    else:
        print("ERROR: Unrecognised output file type: " + args.plotFile + "!")
        sys.exit(1)
else:
    createPlot = None

with open(args.runParamsFile, "r") as runParamsFile:
    runParams = json.load(runParamsFile)

tpcLength = runParams["driftLength"]
tpcRadius = runParams["tpcRadius"]
pixelPitch = runParams["pixelPitch"]

chargeData = np.transpose(np.loadtxt(args.hitsFile, delimiter=',', skiprows=1))[3]
chargeScale = 1. / (np.mean(chargeData) + np.std(chargeData))
scaleFactor = 1. * pixelPitch * chargeScale

#### disable automatic camera reset on 'Show'
pv._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
hitsCsv = pv.CSVReader(FileName=[args.hitsFile])

# create a new 'Table To Points'
tableToPoints1 = pv.TableToPoints(Input=hitsCsv)
tableToPoints1.XColumn = 'X'
tableToPoints1.YColumn = 'Y'
tableToPoints1.ZColumn = 'Z'

# find view
renderView1 = pv.FindViewOrCreate('RenderView1', viewtype='RenderView')

# set active view
pv.SetActiveView(renderView1)

# create a new 'Glyph'
# sphere
glyph1 = pv.Glyph(Input=tableToPoints1, GlyphType='Sphere')
glyph1.Scalars = ['POINTS', 'Q']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = scaleFactor
glyph1.GlyphMode = 'All Points'
#glyph1.GlyphType.ThetaResolution = 8
#glyph1.GlyphType.PhiResolution = 8

# get color transfer function/color map for args.colourColumn
cLUT = pv.GetColorTransferFunction(args.colourColumn)
nc = vtk.vtkNamedColors()
rgb = [0. for i in range(3)]
if args.colourColumn == 'Q':
    cLUT.ApplyPreset("Plasma (matplotlib)")
else:
    nc.GetColorRGB("Green", rgb)
    rgbPoints = [0.] + rgb  # Accepted
    nc.GetColorRGB("Magenta", rgb)
    rgbPoints += [1.] + rgb  # Rejected Hit
    nc.GetColorRGB("Blue", rgb)
    rgbPoints += [2.] + rgb  # Rejected Ambiguity
    cLUT.RGBPoints = rgbPoints

# show data in view
glyph1Display = pv.Show(glyph1, renderView1)
glyph1Display.ColorArrayName = ['POINTS', args.colourColumn]
glyph1Display.LookupTable = cLUT
if args.pcaFile:
    glyph1Display.Opacity = .1

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# set active source
pv.SetActiveSource(None)

cylinder1 = pv.Cylinder()

# Properties modified on cylinder1
cylinder1.Resolution = 60
cylinder1.Height = tpcLength
cylinder1.Radius = tpcRadius

# show data in view
cylinder1Display = pv.Show(cylinder1, renderView1)
cylinder1Display.Orientation = [90.0, 0.0, 0.0]
cylinder1Display.Opacity = 0.05
nc.GetColorRGB("Yellow", rgb)
cylinder1Display.DiffuseColor = rgb

if args.pcaFile:
    pca = np.loadtxt(args.pcaFile, delimiter=",")
    avePos = pca[0]
    direction = pca[1]
    direction /= np.linalg.norm(direction)
    relOffset = avePos[2] / direction[2]
    cylPos = avePos - relOffset * direction

    cylinder2 = pv.Cylinder()
    cylinder2.Resolution = 60
    cylinder2.Height = tpcLength
    cylinder2.Radius = pixelPitch / 5.
    cylinder2Display = pv.Show(cylinder2, renderView1)
    if args.colourColumn == 'Q':
        nc.GetColorRGB("Lime", rgb)
    else:
        nc.GetColorRGB("Red", rgb)
    cylinder2Display.DiffuseColor = rgb
    cylinder2Display.Orientation = [np.rad2deg(np.arcsin(direction[2])),
                                    0.,
                                    (np.rad2deg(np.arctan2(direction[1], direction[0])) - 90.)]
    cylinder2Display.Position = cylPos

renderView1.Background = [0.67, 0.67, 0.67]
renderView1.Update()

viewAngle = 10.
renderView1.CameraViewAngle = viewAngle
renderView1.CameraPosition = [(-1.1 * tpcLength / (2. * np.tan(np.deg2rad(viewAngle / 2.)))), 0., 0.]
renderView1.CameraFocalPoint = [0., 0., 0.]
renderView1.CameraViewUp = [0.0, 0.0, -1.0]
renderView1.ViewSize = [800, 1000]

if createPlot == "view":
    pv.ExportView(args.plotFile,
            view=renderView1)
elif createPlot == "screenshot":
    pv.SaveScreenshot(args.plotFile,
            magnification=2,
            quality=100,
            view=renderView1)
else:
    pv.Interact(view=renderView1)

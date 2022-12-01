# script-version: 2.0
# Catalyst state generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1326, 1117]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [-0.006819963455200195, -0.006789207458496094, 0.9329147338867188]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-30.323604248818963, -30.995370661316226, 34.63367653619927]
renderView1.CameraFocalPoint = [-0.006819963455200108, -0.0067892074584954614, 0.9329147338867205]
renderView1.CameraViewUp = [0.14609634768515117, 0.6592351860978479, 0.7376075017267621]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 14.211831152800906
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.Shadows = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1326, 1117)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
meshpvd = PVDReader(registrationName='mesh', FileName='/Users/vmateevitsi/Desktop/NekRS/post/polaris/pb146/2022.09.13/mesh.pvd')
meshpvd.PointArrays = ['temperature', 'pressure', 'velocity_x', 'velocity_y', 'velocity_z']

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=meshpvd)
calculator1.ResultArrayName = 'velocity'
calculator1.Function = '(velocity_x*iHat)+(velocity_y*jHat)+(velocity_z*kHat)'

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=calculator1)
contour1.ContourBy = ['POINTS', 'temperature']
contour1.Isosurfaces = [40.0, 50.0, 60.0]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=calculator1,
    SeedType='Point Cloud')
streamTracer1.Vectors = ['POINTS', 'velocity']
streamTracer1.MaximumStreamlineLength = 19.660363741904426

# init the 'Point Cloud' selected for 'SeedType'
streamTracer1.SeedType.Center = [-0.00018615544906097625, -0.00015386535008010682, 0.9222175736278677]
streamTracer1.SeedType.NumberOfPoints = 1000
streamTracer1.SeedType.Radius = 6.0

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=streamTracer1,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'velocity']
glyph1.ScaleArray = ['POINTS', 'velocity']
glyph1.ScaleFactor = 0.5
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'Every Nth Point'
glyph1.Stride = 100

# create a new 'Clean to Grid'
cleantoGrid1 = CleantoGrid(registrationName='CleantoGrid1', Input=meshpvd)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from streamTracer1
streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'pressure'
pressureLUT = GetColorTransferFunction('pressure')
pressureLUT.RGBPoints = [-5.3056056129471285, 0.231373, 0.298039, 0.752941, -1.6962987647850474, 0.865003, 0.865003, 0.865003, 1.9130080833770338, 0.705882, 0.0156863, 0.14902]
pressureLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'
streamTracer1Display.ColorArrayName = ['POINTS', 'pressure']
streamTracer1Display.LookupTable = pressureLUT
streamTracer1Display.SelectTCoordArray = 'None'
streamTracer1Display.SelectNormalArray = 'None'
streamTracer1Display.SelectTangentArray = 'None'
streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer1Display.SelectOrientationVectors = 'Normals'
streamTracer1Display.ScaleFactor = 1.9658089637756349
streamTracer1Display.SelectScaleArray = 'AngularVelocity'
streamTracer1Display.GlyphType = 'Arrow'
streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracer1Display.GaussianRadius = 0.09829044818878174
streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
streamTracer1Display.OSPRayScaleFunction.Points = [0.8387917876243591, 0.0, 0.5, 0.0, 1.2149438656535572, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer1Display.ScaleTransferFunction.Points = [-38.91939354593831, 0.0, 0.5, 0.0, 77.80703505194674, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer1Display.OpacityTransferFunction.Points = [-38.91939354593831, 0.0, 0.5, 0.0, 77.80703505194674, 1.0, 0.5, 0.0]

# show data from glyph1
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('velocity')
velocityLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.7180433714489944, 0.865003, 0.865003, 0.865003, 3.436086742897989, 0.705882, 0.0156863, 0.14902]
velocityLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'velocity']
glyph1Display.LookupTable = velocityLUT
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.OSPRayScaleArray = 'AngularVelocity'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'Normals'
glyph1Display.ScaleFactor = 35.117750549316405
glyph1Display.SelectScaleArray = 'AngularVelocity'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'AngularVelocity'
glyph1Display.GaussianRadius = 1.7558875274658203
glyph1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'AngularVelocity']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [0.8387917876243591, 0.0, 0.5, 0.0, 1.2149438656535572, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-86.67754077607128, 0.0, 0.5, 0.0, 100.891465474227, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-86.67754077607128, 0.0, 0.5, 0.0, 100.891465474227, 1.0, 0.5, 0.0]

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.AnnotationsInitialized = 1
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '1', '2', '3', '0', '1', '2', '3']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['FIELD', 'vtkBlockColors']
contour1Display.LookupTable = vtkBlockColorsLUT
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = 'temperature'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'velocity'
contour1Display.ScaleFactor = 1.1990900596851712
contour1Display.SelectScaleArray = 'temperature'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'temperature'
contour1Display.GaussianRadius = 0.05995450298425856
contour1Display.SetScaleArray = ['POINTS', 'temperature']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'temperature']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour1Display.OSPRayScaleFunction.Points = [0.8387917876243591, 0.0, 0.5, 0.0, 1.2149438656535572, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-1.6962987647850474, 0.0, 0.5, 0.0, -1.6960545778274536, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-1.6962987647850474, 0.0, 0.5, 0.0, -1.6960545778274536, 1.0, 0.5, 0.0]

# show data from cleantoGrid1
cleantoGrid1Display = Show(cleantoGrid1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cleantoGrid1Display.Representation = 'Surface'
cleantoGrid1Display.ColorArrayName = ['POINTS', '']
cleantoGrid1Display.Opacity = 0.1
cleantoGrid1Display.RenderPointsAsSpheres = 1
cleantoGrid1Display.SelectTCoordArray = 'None'
cleantoGrid1Display.SelectNormalArray = 'None'
cleantoGrid1Display.SelectTangentArray = 'None'
cleantoGrid1Display.OSPRayScaleArray = 'pressure'
cleantoGrid1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cleantoGrid1Display.SelectOrientationVectors = 'None'
cleantoGrid1Display.ScaleFactor = 1.9660363741904427
cleantoGrid1Display.SelectScaleArray = 'None'
cleantoGrid1Display.GlyphType = 'Arrow'
cleantoGrid1Display.GlyphTableIndexArray = 'None'
cleantoGrid1Display.GaussianRadius = 0.09830181870952213
cleantoGrid1Display.SetScaleArray = ['POINTS', 'pressure']
cleantoGrid1Display.ScaleTransferFunction = 'PiecewiseFunction'
cleantoGrid1Display.OpacityArray = ['POINTS', 'pressure']
cleantoGrid1Display.OpacityTransferFunction = 'PiecewiseFunction'
cleantoGrid1Display.DataAxesGrid = 'GridAxesRepresentation'
cleantoGrid1Display.PolarAxes = 'PolarAxesRepresentation'
cleantoGrid1Display.ScalarOpacityUnitDistance = 0.10115246366203111
cleantoGrid1Display.OpacityArrayName = ['POINTS', 'pressure']
cleantoGrid1Display.ExtractedBlockIndex = 2

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cleantoGrid1Display.OSPRayScaleFunction.Points = [0.8387917876243591, 0.0, 0.5, 0.0, 1.2149438656535572, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cleantoGrid1Display.ScaleTransferFunction.Points = [-5.3056056129471285, 0.0, 0.5, 0.0, 1.9130080833770338, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cleantoGrid1Display.OpacityTransferFunction.Points = [-5.3056056129471285, 0.0, 0.5, 0.0, 1.9130080833770338, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for vtkBlockColorsLUT in view renderView1
vtkBlockColorsLUTColorBar = GetScalarBar(vtkBlockColorsLUT, renderView1)
vtkBlockColorsLUTColorBar.WindowLocation = 'UpperLeftCorner'
vtkBlockColorsLUTColorBar.Title = 'vtkBlockColors'
vtkBlockColorsLUTColorBar.ComponentTitle = ''

# set color bar visibility
vtkBlockColorsLUTColorBar.Visibility = 1

# get color legend/bar for pressureLUT in view renderView1
pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)
pressureLUTColorBar.Title = 'pressure'
pressureLUTColorBar.ComponentTitle = ''

# set color bar visibility
pressureLUTColorBar.Visibility = 1

# get color legend/bar for velocityLUT in view renderView1
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
velocityLUTColorBar.WindowLocation = 'UpperRightCorner'
velocityLUTColorBar.Title = 'velocity'
velocityLUTColorBar.ComponentTitle = 'Magnitude'

# set color bar visibility
velocityLUTColorBar.Visibility = 1

# show color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')
velocityPWF.Points = [0.0, 0.0, 0.5, 0.0, 3.436086742897989, 1.0, 0.5, 0.0]
velocityPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'pressure'
pressurePWF = GetOpacityTransferFunction('pressure')
pressurePWF.Points = [-5.3056056129471285, 0.0, 0.5, 0.0, 1.9130080833770338, 1.0, 0.5, 0.0]
pressurePWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')
vtkBlockColorsPWF.Points = [0.8387917876243591, 0.0, 0.5, 0.0, 1.2149438656535572, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
jPG1 = CreateExtractor('JPG', renderView1, registrationName='JPG1')
# trace defaults for the extractor.
# init the 'JPG' selected for 'Writer'
jPG1.Writer.FileName = 'RenderViewBlock1_%.6ts%cm.jpg'
jPG1.Writer.ImageResolution = [1920, 1080]
jPG1.Writer.Format = 'JPEG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(contour1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)

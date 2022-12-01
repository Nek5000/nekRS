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
renderView1.ViewSize = [1193, 741]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.0, 0.0, 10.000000000000004]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [17.639524985899328, 2.940458190596826, 44.35819099638215]
renderView1.CameraFocalPoint = [0.0, 0.0, 10.000000000000004]
renderView1.CameraViewUp = [-0.02922701018963403, 0.9970956092164056, -0.07032871360077621]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 10.024968827881715
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1193, 741)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
meshpvd = PVDReader(registrationName='mesh')
meshpvd.PointArrays = ['pressure', 'velocity_x', 'velocity_y', 'velocity_z']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from meshpvd
meshpvdDisplay = Show(meshpvd, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'velocity_x'
velocity_xLUT = GetColorTransferFunction('velocity_x')
velocity_xLUT.RGBPoints = [4.5434101973085134e-11, 0.231373, 0.298039, 0.752941, 0.024999517000358726, 0.865003, 0.865003, 0.865003, 0.04999903395528335, 0.705882, 0.0156863, 0.14902]
velocity_xLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'velocity_x'
velocity_xPWF = GetOpacityTransferFunction('velocity_x')
velocity_xPWF.Points = [4.5434101973085134e-11, 0.0, 0.5, 0.0, 0.04999903395528335, 1.0, 0.5, 0.0]
velocity_xPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
meshpvdDisplay.Representation = 'Surface'
meshpvdDisplay.ColorArrayName = ['POINTS', 'velocity_x']
meshpvdDisplay.LookupTable = velocity_xLUT
meshpvdDisplay.RenderPointsAsSpheres = 1
meshpvdDisplay.SelectTCoordArray = 'None'
meshpvdDisplay.SelectNormalArray = 'None'
meshpvdDisplay.SelectTangentArray = 'None'
meshpvdDisplay.OSPRayScaleArray = 'pressure'
meshpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
meshpvdDisplay.SelectOrientationVectors = 'None'
meshpvdDisplay.ScaleFactor = 2.000000000000001
meshpvdDisplay.SelectScaleArray = 'None'
meshpvdDisplay.GlyphType = 'Arrow'
meshpvdDisplay.GlyphTableIndexArray = 'None'
meshpvdDisplay.GaussianRadius = 0.10000000000000003
meshpvdDisplay.SetScaleArray = ['POINTS', 'pressure']
meshpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
meshpvdDisplay.OpacityArray = ['POINTS', 'pressure']
meshpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
meshpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
meshpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
meshpvdDisplay.ScalarOpacityFunction = velocity_xPWF
meshpvdDisplay.ScalarOpacityUnitDistance = 0.15363493625908156
meshpvdDisplay.OpacityArrayName = ['POINTS', 'pressure']

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
meshpvdDisplay.OSPRayScaleFunction.Points = [0.8387917876243591, 0.0, 0.5, 0.0, 1.2149438656535572, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
meshpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
meshpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocity_xLUT in view renderView1
velocity_xLUTColorBar = GetScalarBar(velocity_xLUT, renderView1)
velocity_xLUTColorBar.Title = 'velocity_x'
velocity_xLUTColorBar.ComponentTitle = ''

# set color bar visibility
velocity_xLUTColorBar.Visibility = 1

# show color legend
meshpvdDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
jPG1 = CreateExtractor('JPG', renderView1, registrationName='JPG1')
# trace defaults for the extractor.
# init the 'JPG' selected for 'Writer'
jPG1.Writer.FileName = 'RenderView1_%.6ts%cm.jpg'
jPG1.Writer.ImageResolution = [1193, 741]
jPG1.Writer.Format = 'JPEG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(jPG1)
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

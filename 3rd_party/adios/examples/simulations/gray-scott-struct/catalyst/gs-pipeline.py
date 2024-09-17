from paraview.simple import *
from paraview import print_info

# catalyst options
from paraview.catalyst import Options

options = Options()

# print start marker
print_info("begin '%s'", __name__)

# directory under which to save all extracts
# generated using Extractors defined in the pipeline
# (optional, but recommended)
options.ExtractsOutputDirectory = "."
SaveExtractsUsingCatalystOptions(options)

view = CreateRenderView()
# when using Fides, registrationName is always 'fides'
producer = TrivialProducer(registrationName="fides")
display = Show(producer)
view.ResetCamera()
ColorBy(display, ("POINTS", "U"))
display.RescaleTransferFunctionToDataRange(True, False)
display.SetScalarBarVisibility(view, True)
transFunc = GetColorTransferFunction("U")
transFunc.RescaleOnVisibilityChange = 1

display.SetRepresentationType("Surface")

clip = Clip(registrationName="clip1", Input=producer)
clip.ClipType = "Plane"
clip.Scalars = ["POINTS", "U"]
clip.ClipType.Origin = [1.5, 1.5, 1.5]
clip.ClipType.Normal = [0.0, 1.0, 0.0]
clipDisplay = Show(clip, view, "UnstructuredGridRepresentation")

Hide(producer, view)
view.ResetCamera()

camera = GetActiveCamera()
camera.Azimuth(45)
camera.Elevation(45)

# the extractor will save the view on each time step
extractor = CreateExtractor("PNG", view, registrationName="PNG1")
extractor.Writer.FileName = "output-{timestep}.png"
extractor.Writer.ImageResolution = [800, 800]


def catalyst_execute(info):
    print_info("in '%s::catalyst_execute'", __name__)
    global producer
    producer.UpdatePipeline()
    print("updating pipeline and saving image")

    # print("-----------------------------------")
    # print("executing (cycle={}, time={})".format(info.cycle, info.time))
    # print("bounds:", producer.GetDataInformation().GetBounds())
    # print("U-range:", producer.PointData['U'].GetRange(0))
    # print("V-range:", producer.PointData['V'].GetRange(0))


# print end marker
print_info("end '%s'", __name__)

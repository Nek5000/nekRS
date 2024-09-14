import argparse
from paraview.simple import *

from paraview import print_info

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst

options = catalyst.Options()
options.GlobalTrigger = "TimeStep"
options.EnableCatalystLive = 1
options.CatalystLiveTrigger = "TimeStep"
options.ExtractsOutputDirectory = "/tmp"
# options.ExtractsOutputDirectory = '.'

DISABLE_EXTRACTOR = False

# setup and returns the view
def SetupRenderView():
    # Create a new 'Render View'
    view = CreateView("RenderView")

    camera = GetActiveCamera()
    camera.Azimuth(45)
    camera.Elevation(45)

    SetActiveView(None)

    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name="Layout #1")
    layout1.AssignView(0, view)
    layout1.SetSize(824, 656)

    # restore active view
    SetActiveView(view)
    return view


# A different proxy type is used depending on whether you're using
# catalyst, or doing post hoc visualization
# this returns the catalyst producer
def SetupCatalystProducer():
    producer = TrivialProducer(registrationName="fides")
    return producer


# this returns the fides proxy which should be used for post hoc vis
# i.e., you're reading .bp files
def SetupFidesReader(json, bp, sst):
    if json is None:
        # in this case the bp file must contain the Fides attributes
        fides = FidesReader(StreamSteps=1, FileName=bp)
        return fides

    fides = FidesJSONReader(StreamSteps=1, FileName=json)
    if sst:
        fides.DataSourceEngines = ["source", "SST"]
    # 'source' is the name of the ADIOS data source in the JSON data model
    fides.DataSourcePath = ["source", bp]
    # required to update the fides reader
    fides.UpdatePipelineInformation()
    return fides


# takes in a producer and view and sets up the visualization pipeline
def SetupVisPipeline(producer, view):
    view.ResetCamera()

    # get color transfer function/color map for 'U'
    uLUT = GetColorTransferFunction("V")
    uLUT.AutomaticRescaleRangeMode = "Clamp and update every timestep"
    uLUT.RescaleOnVisibilityChange = 1

    contour1 = Contour(registrationName="Contour1", Input=producer)
    contour1.ContourBy = ["POINTS", "V"]
    contour1.Isosurfaces = [0.1, 0.3, 0.5, 0.7]
    contour1.PointMergeMethod = "Uniform Binning"

    contour1Display = Show(contour1, view, "GeometryRepresentation")

    contour1Display.Representation = "Surface"
    contour1Display.ColorArrayName = ["POINTS", "V"]

    Hide(producer, view)

    # trace defaults for the display properties.
    contour1Display.Representation = "Surface"
    contour1Display.ColorArrayName = ["POINTS", "V"]
    contour1Display.SetScaleArray = ["POINTS", "V"]
    contour1Display.ScaleTransferFunction = "PiecewiseFunction"
    contour1Display.LookupTable = uLUT

    uLUTColorBar = GetScalarBar(uLUT, view)
    uLUTColorBar.Title = "V"
    uLUTColorBar.ComponentTitle = ""

    # show color legend
    contour1Display.SetScalarBarVisibility(view, True)
    return contour1, contour1Display


# sets up an extractor for writing out PNG images
def SetupExtractor(view):
    # create extractor
    pNG1 = CreateExtractor("PNG", view, registrationName="PNG1")
    # trace defaults for the extractor.
    pNG1.Trigger = "TimeStep"

    # init the 'PNG' selected for 'Writer'
    pNG1.Writer.FileName = "output_{timestep:06d}.png"
    pNG1.Writer.ImageResolution = [800, 800]
    pNG1.Writer.Format = "PNG"


# Catalyst uses this to update the pipeline and print out some info on each time step
# you don't need to call this directly in your script; ParaView Catalyst will call it for you
def catalyst_execute(info):
    print_info("in '%s::catalyst_execute'", __name__)
    global pipeline, display
    pipeline.UpdatePipeline()
    display.RescaleTransferFunctionToDataRange()

    print_info("executing (cycle={}, time={})".format(info.cycle, info.time))
    print_info("U-range: {}".format(producer.PointData["U"].GetRange(0)))
    print_info("V-range: {}".format(producer.PointData["V"].GetRange(0)))


def ParseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-j",
        "--json_filename",
        help="path to Fides JSON file",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-b", "--bp_filename", help="path to bp file", type=str, required=True
    )
    parser.add_argument("--staging", help="use SST engine", action="store_true")
    args = parser.parse_args()
    return args


def StreamingVis(args):
    # adios/fides step status
    NotReady = 1
    EndOfStream = 2

    # setup the reader, view, pipeline
    fides = SetupFidesReader(args.json_filename, args.bp_filename, args.staging)
    view = SetupRenderView()

    step = 0
    while True:
        status = NotReady
        while status == NotReady:
            # must call PrepareNextStep to get Fides ready to read the
            # next step
            fides.PrepareNextStep()
            fides.UpdatePipelineInformation()
            status = fides.NextStepStatus
        if status == EndOfStream:
            # done reading the file
            return
        if step == 0:
            # set up the pipeline on the first step
            pipeline, display = SetupVisPipeline(fides, view)

        # need to update the pipeline and then save the output
        pipeline.UpdatePipeline()
        display.RescaleTransferFunctionToDataRange()
        output = f"output-{step:05d}.png"
        SaveScreenshot(output, view, ImageResolution=[800, 800])
        step += 1


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    print("in __main__()")
    args = ParseArgs()
    StreamingVis(args)
else:
    # in this case we're running from Catalyst
    view = SetupRenderView()
    producer = SetupCatalystProducer()
    pipeline, display = SetupVisPipeline(producer, view)

    # normally not needed, but a bug fix is in progress to fix an issue
    # when using ParaView Live with extractors
    if not DISABLE_EXTRACTOR:
        SetupExtractor(view)

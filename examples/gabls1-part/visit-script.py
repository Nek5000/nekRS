#! visit -cli -nowin -s

import sys

## config ##
start_state = 0
stop_state = 100

## Plots particles ##
# Need to use the post-processed version of the particles
OpenDatabase("part*.3d database", 0)
AddPlot("Scatter", "X", 1, 1)
ScatterAtts = ScatterAttributes()
ScatterAtts.var1 = "X"
ScatterAtts.var1Role = ScatterAtts.Coordinate0  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var2 = "Y"
ScatterAtts.var2Role = ScatterAtts.Coordinate1  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var3 = "Z"
ScatterAtts.var3Role = ScatterAtts.Coordinate2  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var4 = "Color"
ScatterAtts.var4Role = ScatterAtts.Color  # Coordinate0, Coordinate1, Coordinate2, Color, None
#ScatterAtts.var4MinFlag = 1
#ScatterAtts.var4Min = particle_color_min
#ScatterAtts.var4MaxFlag = 1
#ScatterAtts.var4Max = particle_color_max
ScatterAtts.pointType = ScatterAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.pointSize = 0.022
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorByColorTable  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.colorTableName = "Paired3"
ScatterAtts.invertColorTable = 1
SetPlotOptions(ScatterAtts)

## Plot fluid data ##
OpenDatabase("gabls-part.nek5000", 0)
AddPlot("Pseudocolor", "temperature")
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = 0.0015
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 0.005
PseudocolorAtts.colorTableName = "gray"
PseudocolorAtts.lightingFlag = 0
SetPlotOptions(PseudocolorAtts)
AddOperator("ThreeSlice", 0)
ThreeSliceAtts = ThreeSliceAttributes()
ThreeSliceAtts.x = 3.99
ThreeSliceAtts.y = 0.10
ThreeSliceAtts.z = 0
SetOperatorOptions(ThreeSliceAtts, 0, 0)

## Correlate fluid and particles ##
CreateDatabaseCorrelation("Correlation01",("gabls-part.nek5000", "part*.3d database"), 0)

w = GetWindowInformation()
print(w.timeSliders)
print(w.activeTimeSlider)
sys.stdout.flush()


## Remove Annotations ##
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes3D.visible = 0
AnnotationAtts.axes3D.triadFlag = 0
AnnotationAtts.legendInfoFlag = 0
AnnotationAtts.axes3D.bboxFlag = 0
AnnotationAtts.axes3D.bboxLocation = (0, 4, 0, 4, 0, 4)
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.databaseInfoFlag = 0
SetAnnotationAttributes(AnnotationAtts)

## Render ##
DrawPlots()

View3DAtts = GetView3D()
View3DAtts.viewNormal = (-0.866135954355074, 0.132289153675643, 0.4819834939012 )
View3DAtts.viewUp = (0.0949584329674016, 0.990326910828711, -0.101170676071859)
View3DAtts.focus = (2, 2, 2)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 3.4641
View3DAtts.nearPlane = -6.9282
View3DAtts.farPlane = 6.9282
View3DAtts.imagePan = (0.00863090175204485, 0.0675639974313456)
View3DAtts.imageZoom = 1.0
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (2, 2, 2)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)

SaveAtts = SaveWindowAttributes()
SaveAtts.format = SaveAtts.PNG
SaveAtts.fileName = "gabls-part"
SaveAtts.width=2048
SaveAtts.height=2048
SetSaveWindowAttributes(SaveAtts)

DrawPlots()
for state in range(start_state, stop_state):
  SetTimeSliderState(state)
  SaveWindow()





The fluid and particles can be rendered using VisIT and the included script.
First, the `visnek` from Nek5000 to generate `gabls-part.nek5000`.  Next, run
`visit-script.py` with VisIT.  This is done through the command line by running
`visit -cli -nowin -s visit-script.py` from this directory.  To use the script
from the graphical interface, choose the "Controls>Launch CLI" menu option, then
use the prompt to evaluate `Source("/path/to/NekRS/examples/gabls-part/visit-script.py")`.

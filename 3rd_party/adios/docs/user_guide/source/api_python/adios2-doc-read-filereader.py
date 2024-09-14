import numpy as np
from adios2 import FileReader

with FileReader("cfd.bp") as s:
    # inspect variables
    vars = s.available_variables()
    for name, info in vars.items():
        print("variable_name: " + name, end=" ")
        for key, value in info.items():
            print("\t" + key + ": " + value, end=" ")
        print()
    print()

    nproc = s.read("nproc")
    print(f"nproc is {nproc} of type {type(nproc)} with ndim {nproc.ndim}")

    # read variables return a numpy array with corresponding selection
    steps = int(vars["physical_time"]["AvailableStepsCount"])
    physical_time = s.read("physical_time", step_selection=[0, steps])
    print(
        f"physical_time is {physical_time} of type {type(physical_time)} with "
        f"ndim {physical_time.ndim} shape = {physical_time.shape}"
    )

    steps = int(vars["temperature"]["AvailableStepsCount"])
    temperature = s.read("temperature", step_selection=[0, steps])
    temp_unit = s.read_attribute("temperature/unit")
    print(f"temperature array size is {temperature.size} of shape {temperature.shape}")
    print(f"temperature unit is {temp_unit} of type {type(temp_unit)}")

    steps = int(vars["pressure"]["AvailableStepsCount"])
    pressure = s.read("pressure", step_selection=[0, steps])
    press_unit = s.read_attribute("pressure/unit")

    print()

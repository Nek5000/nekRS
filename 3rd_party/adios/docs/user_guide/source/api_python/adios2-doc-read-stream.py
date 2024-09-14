import numpy as np
from adios2 import Stream

with Stream("cfd.bp", "r") as s:
    # steps comes from the stream
    for _ in s.steps():
        # track current step
        print(f"Current step is {s.current_step()}")

        # inspect variables in current step
        for name, info in s.available_variables().items():
            print("variable_name: " + name, end=" ")
            for key, value in info.items():
                print("\t" + key + ": " + value, end=" ")
            print()

        if s.current_step() == 0:
            nproc = s.read("nproc")
            print(f"nproc is {nproc} of type {type(nproc)} with ndim {nproc.ndim}")

        # read variables return a numpy array with corresponding selection
        physical_time = s.read("physical_time")
        print(f"physical_time is {physical_time} of type {type(physical_time)}")
        temperature = s.read("temperature")
        temp_unit = s.read_attribute("temperature/unit")
        print(f"temperature array size is {temperature.size} of shape {temperature.shape}")
        print(f"temperature unit is {temp_unit} of type {type(temp_unit)}")
        pressure = s.read("pressure")
        press_unit = s.read_attribute("unit", "pressure")
        print(f"pressure unit is {press_unit} of type {type(press_unit)}")
        print()

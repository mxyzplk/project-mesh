# -*- coding: utf-8 -*-
from config import Config
from surface import Surface

config_file = "config.txt"

conf = Config(config_file)

imesh = Surface(conf.i_mesh_type)
omesh = Surface(conf.o_mesh_type)

# Input mesh reading
imesh.read_grids(conf.data_dir, conf.i_grids_file)
imesh.read_elements(conf.data_dir, conf.i_elements_file)
imesh.calc_area()
imesh.allocate_press(conf.np)

# Output mesh reading
omesh.read_grids(conf.data_dir, conf.o_grids_file)
omesh.read_elements(conf.data_dir, conf.o_elements_file)
omesh.mapgrids(imesh.centers, imesh.nelements, conf.plane)
omesh.calc_area()


# Input pressures reading
for i in range(conf.np):
    imesh.read_press(conf.data_dir, conf.data['files'][i], conf.i_press_type)

    if (conf.i_press_type == 1):    # pressure on elements
        omesh.projectmesh(imesh.centers, imesh.area, imesh.press, conf.plane)
    elif (conf.i_press_type == 0):   # pressure on grids
        imesh.pressure_on_elements_centers()
        omesh.projectmesh(imesh.centers, imesh.area, imesh.pressures_on_centers, conf.plane)
        
    fout = "forces_" + str(i + 1) + ".txt"
    omesh.write_projected_mesh(fout)
    
    fout = "forces_" + str(i + 1) + "_input_int.txt"
    imesh.intmesh(fout)
    
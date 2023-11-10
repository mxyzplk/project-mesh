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
imesh.allocate_press(conf.np)

# Output mesh reading
omesh.read_grids(conf.data_dir, conf.o_grids_file)
omesh.read_elements(conf.data_dir, conf.o_elements_file)

# Input pressures reading
for i in range(conf.np):
    imesh.read_press(conf.data_dir, conf.data['files'][i], conf.i_press_type)
    omesh.projectmesh(imesh.elements, imesh.area, imesh.press, conf.plane)
    fout = "forces_" + str(i)
    omesh.write_projected_mesh(fout)
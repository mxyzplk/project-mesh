# -*- coding: utf-8 -*-
import yaml
import os

class Config:
    def __init__(self, filename):

        # directories
        self.main_dir = os.path.dirname(os.path.abspath(__file__))
        self.resources_dir = os.path.join(self.main_dir, '../resources')
        self.results = os.path.join(self.main_dir, '../results')
        self.configpath = os.path.join(self.resources_dir, filename)        
        
        with open(self.configpath, 'r') as file:
            self.data = yaml.load(file, Loader=yaml.FullLoader)
        
        # Number of pressure files
        self.np = len(self.data['files'])
                        
        self.name = self.data['component']['name']
        self.i_grids_file = self.data['component']['input_grids']
        self.i_elements_file = self.data['component']['input_elements']
        self.o_grids_file = self.data['component']['output_grids']
        self.o_elements_file = self.data['component']['output_elements']        
        self.i_mesh_type = self.data['component']['input_mesh_type']
        self.i_press_type = self.data['component']['input_press_type']
        self.o_mesh_type = self.data['component']['output_mesh_type']
        self.o_press_type = self.data['component']['output_press_type']        
        self.data_dir = self.data['component']['data_dir']
        self.plane = self.data['component']['plane']
            
        
        

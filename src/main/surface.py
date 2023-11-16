# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy.spatial.qhull import Delaunay


class Surface:
    def __init__(self, mesh_type):
        
        # data
        self.ngrids = 0
        self.nelements = 0
        self.mesh_type = int(mesh_type)
        
        # directories
        self.main_dir = os.path.dirname(os.path.abspath(__file__))
        self.resources_dir = os.path.join(self.main_dir, '../resources')
        self.results_dir = os.path.join(self.main_dir, '../results')         

        
    def read_grids(self, filedir, filename):
        
        filepath = os.path.join(self.resources_dir, filedir, filename)
        
        with open(filepath, 'r') as file:
            line = file.readline()
            temp = line.split()
            self.ngrids = int(temp[0])                  # Number of grids
            self.grids = np.empty((self.ngrids, 3))
            self.ids = np.empty(self.ngrids)
            
            for i in range(self.ngrids):
                line = file.readline()
                temp = line.split()
                self.ids[i] = int(temp[0])              # Reading grid id   
                self.grids[i, 0] = float(temp[1])       # Reading grid x coordinate  
                self.grids[i, 1] = float(temp[2])       # Reading grid y coordinate 
                self.grids[i, 2] = float(temp[3])       # Reading grid z coordinate 
                
                
    def read_elements(self, filedir, filename):
        
        filepath = os.path.join(self.resources_dir, filedir, filename)
        
        with open(filepath, 'r') as file:
            line = file.readline()
            temp = line.split()
            self.nelements = int(temp[0])                  # Number of elements
            self.elements = np.empty((self.nelements, self.mesh_type))
            self.centers = np.empty((self.nelements, 3))
            self.area = np.empty((self.nelements, 4))

            for i in range(self.nelements):
                line = file.readline()
                temp = line.split()
                for j in range(self.mesh_type):
                    self.elements[i, j] = int(temp[j])       # Reading elements
                    
        self.calc_panel_centers()
                    
                    
    def read_press(self, filedir, filename, ptype):
        
        filepath = os.path.join(self.resources_dir, filedir, filename)
        
        with open(filepath, 'r') as file:
            if (int(ptype) == 0):  # press on grids
                self.press = np.empty(self.ngrids)
                for i in range(self.ngrids):
                        line = file.readline()
                        temp = line.split()
                        self.press[i] = float(temp[1])
                        
            elif (int(ptype) == 1):  # press on element centers
                self.press = np.empty(self.nelements)
                for i in range(self.nelements):
                        line = file.readline()
                        temp = line.split()
                        self.press[i] = float(temp[1])
     

    def is_point_inside_panel(self, panel_grids, center, plane):
        g1, g2 = 0, 1  # Default values for plane xy

        if int(plane) == 1:  # Adjust values for plane xz
            g2 = 2
            

        grids = panel_grids[:, [g1, g2]]

        point_to_check = center[[g1, g2]]
        
        triangulation = Delaunay(grids)
        simplex_index = triangulation.find_simplex(point_to_check, tol=1e-10)
        
        is_inside = simplex_index != -1
       
        return is_inside


    def mapgrids(self, points, npi, plane):
        
        self.mapg = np.zeros((npi))
        check_points = np.zeros(npi, dtype=bool)
        ch = 0
        self.mapg[:] = -1
        
        for i in range(self.nelements):  # For each cell in target mesh
            panel_grids = self.grids[self.elements[i, :].astype(int) - 1, :]   # grid array        

            for j in range(npi): # Evaluate each point with pressure from the input mesh
                            
                if not check_points[j]:
                    check = self.is_point_inside_panel(panel_grids, points[j, :], plane)
                    
                if (check == True and check_points[j] == False):
                    check_points[j] = True
                    self.mapg[j] = i     
                    ch += 1
                    print(str(ch) + "/" + str(npi))
                    
        filepath = os.path.join(self.results_dir, "map-log.txt")
        
        with open(filepath, 'w') as file:
            file.write("Missing grids\n")
            for j in range(npi):
                if (self.mapg[j] == -1):
                    file.write("{:8d} {:12.6f} {:12.6f} {:12.6f}\n".format(j, points[j, 0], points[j, 1], points[j, 2]))
                      

    def projectmesh(self, points, areas, pressures, plane):
        
        # points: points with pressure in original mesh
        # pressures: pressures from original mesh
        # areas: panel areas from original mesh
        
        self.pcenter = np.zeros(self.nelements)
        self.cforce = np.zeros((self.nelements, 3))
        
        npi = len(pressures)
        
        for i in range(self.nelements):  # For each cell in target mesh
            
            forces = np.zeros(3)
            
            for j in range(npi): # Evaluate each point with pressure from the input mesh
                if (self.mapg[j] == i):
                    self.pcenter[i] += pressures[j]
                    forces[:] += pressures[j] * areas[j, 0:3]
            
            self.cforce[i, :] = forces[:]
            

 

    def write_projected_mesh(self, fout):
        
        filepath = os.path.join(self.results_dir, fout)
        
        fout2 = fout.replace(".txt", "_output_int.txt")
        
        with open(filepath, 'w') as file:
            for i in range(self.nelements):
                file.write("{:12.6f} {:12.6f} {:12.6f} {:.6e} {:.6e} {:.6e} {:.7e} {:.7e} {:.7e} {:.7e}\n".format(self.centers[i, 0], \
                                                                                                                       self.centers[i, 0], \
                                                                                                                       self.centers[i, 0], \
                                                                                                                       self.cforce[i, 0], \
                                                                                                                       self.cforce[i, 1], \
                                                                                                                       self.cforce[i, 2], \
                                                                                                                       self.area[i, 0], \
                                                                                                                       self.area[i, 1], \
                                                                                                                       self.area[i, 2], \
                                                                                                                       self.area[i, 3]))                 
        
        self.intmesh2(fout2)
    
  
     
    def pressure_on_elements_centers(self):
        
        self.pressures_on_centers = np.empty(self.nelements)
        
        for i in range(self.nelements):
            if (self.mesh_type == 3):
                self.pressures_on_centers[i] = (self.grids[int(self.press[i, 0]) - 1] + self.grids[int(self.press[i, 1]) - 1] + self.grids[int(self.press[i, 2]) - 1]) / 3
            elif (self.mesh_type == 4):
                self.pressures_on_centers[i] = 0.25 * (self.grids[int(self.press[i, 0]) - 1] + self.grids[int(self.press[i, 1]) - 1] + self.grids[int(self.press[i, 2]) - 1] + self.grids[int(self.press[i, 3]) - 1])
        
        
    def allocate_press(self, npt):
        self.press = np.zeros(npt)


    def calc_panel_centers(self):
        
        for i in range(self.nelements):
            if (self.mesh_type == 3):
                self.centers[i, 0] = (self.grids[int(self.elements[i, 0]) - 1, 0] + self.grids[int(self.elements[i, 1]) - 1, 0] + self.grids[int(self.elements[i, 2]) - 1, 0]) / 3
                self.centers[i, 1] = (self.grids[int(self.elements[i, 0]) - 1, 1] + self.grids[int(self.elements[i, 1]) - 1, 1] + self.grids[int(self.elements[i, 2]) - 1, 1]) / 3
                self.centers[i, 2] = (self.grids[int(self.elements[i, 0]) - 1, 2] + self.grids[int(self.elements[i, 1]) - 1, 2] + self.grids[int(self.elements[i, 2]) - 1, 2]) / 3
            elif (self.mesh_type == 4):
                self.centers[i, 0] = 0.25 * (self.grids[int(self.elements[i, 0]) - 1, 0] + self.grids[int(self.elements[i, 1]) - 1, 0] + self.grids[int(self.elements[i, 2]) - 1, 0] + self.grids[int(self.elements[i, 3]) - 1, 0])
                self.centers[i, 1] = 0.25 * (self.grids[int(self.elements[i, 0]) - 1, 1] + self.grids[int(self.elements[i, 1]) - 1, 1] + self.grids[int(self.elements[i, 2]) - 1, 1] + self.grids[int(self.elements[i, 3]) - 1, 1])
                self.centers[i, 2] = 0.25 * (self.grids[int(self.elements[i, 0]) - 1, 2] + self.grids[int(self.elements[i, 1]) - 1, 2] + self.grids[int(self.elements[i, 2]) - 1, 2] + self.grids[int(self.elements[i, 3]) - 1, 2])
    
        
    
    def calc_area(self):
         
        for i in range(self.nelements):
            if (self.mesh_type == 3):
                points = tuple(self.grids[[int(self.elements[i, 0]) - 1, int(self.elements[i, 1]) - 1, int(self.elements[i, 2]) - 1], 0:3])
                self.area[i, 0:3] = self.calculate_triangle_area_3d(points)
            elif (self.mesh_type == 4):
                points = tuple(self.grids[[int(self.elements[i, 0]) - 1, int(self.elements[i, 1]) - 1, int(self.elements[i, 2]) - 1, int(self.elements[i, 3]) - 1], 0:3])
                self.area[i, 0:3] = self.calculate_quadrilateral_area_3d(points)
            self.area[i, 3] = (self.area[i, 0] ** 2 + self.area[i, 1] ** 2 + self.area[i, 2] ** 2) ** 0.5
                
                
                
                
                
    def calculate_triangle_area_3d(self, points):
        A, B, C = points
        AB = np.array(B) - np.array(A)
        AC = np.array(C) - np.array(A)  

        area_vector = 0.5 * np.cross(AB, AC)
        
        return area_vector
        


    def calculate_quadrilateral_area_3d(self, points):
        A, B, C, D = points
        AC = np.array(C) - np.array(A)
        BD = np.array(D) - np.array(B)  

        area_vector = 0.5 * np.cross(AC, BD)       
        
        return area_vector
    
    
    
    def write_press_grids(self, fout):
        
        filepath = os.path.join(self.results_dir, fout)
        filepath_vtk = filepath.replace(".txt", ".vtk")
        
        with open(filepath, 'w') as file: 
            for i in range(self.ngrids):
                file.write("{:12.6f}\n".format(self.press[i]))
                
                
        with open(filepath_vtk, 'w') as file:
            
            file.write("# vtk DataFile Version 3.0\n")
            file.write("vtk output\n")
            file.write("ASCII\n")
            file.write("DATASET UNSTRUCTURED_GRID\n\n")
            file.write("POINTS {:8d} float\n".format(self.ngrids))
            
            for i in range(int(self.ngrids)):
                file.write("{:12.6f} {:12.6f} {:12.6f}\n".format(self.grids[i, 0], self.grids[i, 1], self.grids[i, 2]))
                
            file.write("\n")
            file.write("CELLS {:8d} {:8d}\n".format(self.nelements, self.nelements * (int(self.mesh_type + 1))))
            
            for i in range(self.nelements):
                if (int(self.mesh_type) == 3):
                    file.write("{:8d} {:8d} {:8d} {:8d}\n".format(int(self.mesh_type), self.elements[i, 0], self.elements[i, 1], self.elements[i, 2] ))
                elif (int(self.mesh_type) == 4):
                    file.write("{:8d} {:8d} {:8d} {:8d} {:8d}\n".format(int(self.mesh_type), int(self.elements[i, 0])-1, int(self.elements[i, 1])-1, int(self.elements[i, 2])-1, int(self.elements[i, 3])-1))
                    
            file.write("\n")
            file.write("CELL_TYPES {:8d}\n".format(self.nelements))
            
            for i in range(self.nelements):
                file.write("9\n")

            file.write("\n")
            
            file.write("POINT_DATA {:8d}\n".format(self.ngrids))
            file.write("SCALARS dcp float\n")
            file.write("LOOKUP_TABLE default\n")
            file.write("\n")
            for i in range(self.ngrids):
                file.write("{:12.6f} \n".format(self.press[i]))
            
    def write_press_elements(self, fout):
        
        filepath = os.path.join(self.results_dir, fout)
        filepath_vtk = filepath.replace(".txt", ".vtk")
        
        with open(filepath, 'w') as file: 
            for i in range(self.nelements):
                file.write("{:12.6f}\n".format(self.press[i])) 

        with open(filepath_vtk, 'w') as file:
            
            file.write("# vtk DataFile Version 3.0\n")
            file.write("vtk output\n")
            file.write("ASCII\n")
            file.write("DATASET UNSTRUCTURED_GRID\n\n")
            file.write("POINTS {:8d} float\n".format(self.ngrids))
            
            for i in range(int(self.ngrids)):
                file.write("{:12.6f} {:12.6f} {:12.6f}\n".format(self.grids[i, 0], self.grids[i, 1], self.grids[i, 2]))
                
            file.write("\n")
            file.write("CELLS {:8d} {:8d}\n".format(self.nelements, self.nelements * (int(self.mesh_type + 1))))
            
            for i in range(self.nelements):
                if (int(self.mesh_type) == 3):
                    file.write("{:8d} {:8d} {:8d} {:8d}\n".format(int(self.mesh_type), int(self.elements[i, 0]) - 1, int(self.elements[i, 1]) - 1 , int(self.elements[i, 2]) - 1))
                elif (int(self.mesh_type) == 4):
                    file.write("{:8d} {:8d} {:8d} {:8d} {:8d}\n".format(int(self.mesh_type), int(self.elements[i, 0]) - 1, int(self.elements[i, 1]) - 1, int(self.elements[i, 2]) - 1, int(self.elements[i, 3]) - 1))
                    
            file.write("\n")
            file.write("CELL_TYPES {:8d}\n".format(self.nelements))
            
            for i in range(self.nelements):
                file.write("9\n")

            file.write("\n")
            
            file.write("CELL_DATA {:8d}\n".format(self.nelements))
            file.write("SCALARS dcp float\n")
            file.write("LOOKUP_TABLE default\n")
            file.write("\n")
            for i in range(self.nelements):
                file.write("{:12.6f} \n".format(self.press[i]))
                
                

    def intmesh(self, fout):
        
        filepath = os.path.join(self.results_dir, fout)
        self.forces = np.zeros((self.nelements, 3))
        self.moments = np.zeros((self.nelements, 3))
        self.distance = np.zeros((self.nelements, 3))
        self.refpoint = np.zeros(3)
        self.refcs = np.zeros((4, 3))
        self.refcsvec = np.zeros((3, 3))
        self.rotated_forces = np.zeros((self.nelements, 3))
        self.rotated_moments = np.zeros((self.nelements, 3))
        self.integrated_forces = np.zeros(3)
        self.integrated_moments = np.zeros(3)
        
        cspath = os.path.join(self.resources_dir, "cs.txt")

        with open(cspath, 'r') as f1:
            #
            line = f1.readline()
            #
            line = f1.readline()
            temp = line.split()
            self.refpoint[0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()            
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[0, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[1, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[2, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[3, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            #
            line = f1.readline()
            temp = line.split()            
            self.aref = float(temp[0])
            #
            line = f1.readline()
            temp = line.split()            
            self.cref = float(temp[0])
            #
            line = f1.readline()
            temp = line.split()            
            self.bref = float(temp[0])            
                                                
            
        for i in range(3):
            self.refcsvec[0, i] = self.refcs[1, i] - self.refpoint[i]    # Coordinate System Vectors
            self.refcsvec[1, i] = self.refcs[2, i] - self.refpoint[i]
            self.refcsvec[2, i] = self.refcs[3, i] - self.refpoint[i]


        for i in range(self.nelements):
            self.forces[i, 0] = self.press[i] * self.area[i, 0] # Fx
            self.forces[i, 1] = self.press[i] * self.area[i, 1] # Fy
            self.forces[i, 2] = self.press[i] * self.area[i, 2] # Fz
            for j in range(3):
                self.distance[i, j] = self.centers[i, j] - self.refpoint[j]
            
            #a = np.cross(self.forces[i, 0:3], self.distance)
            #print(a)
            self.moments[i, 0:3] = np.cross(self.forces[i, 0:3], self.distance[i, 0:3])  # Mx, My, Mz
            
            #b = np.dot(self.refcsvec, self.forces)
            #print(b)
            self.rotated_forces[i, 0:3] = np.dot(self.refcsvec, self.forces[i, 0:3])
            self.rotated_moments[i, 0:3] = np.dot(self.refcsvec, self.moments[i, 0:3])
        

        self.integrated_forces[0] = np.sum(self.rotated_forces[:, 0]) / self.aref
        self.integrated_forces[1] = np.sum(self.rotated_forces[:, 1]) / self.aref
        self.integrated_forces[2] = np.sum(self.rotated_forces[:, 2]) / self.aref
        
        self.integrated_moments[0] = np.sum(self.rotated_forces[:, 0]) / (self.aref * self.bref)
        self.integrated_moments[1] = np.sum(self.rotated_forces[:, 1]) / (self.aref * self.cref)
        self.integrated_moments[2] = np.sum(self.rotated_forces[:, 2]) / (self.aref * self.bref)
        
        
        with open(filepath, 'w') as f2:    
            f2.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(self.integrated_forces[0], \
                                                                                      self.integrated_forces[1], \
                                                                                      self.integrated_forces[2], \
                                                                                      self.integrated_moments[0], \
                                                                                      self.integrated_moments[1], \
                                                                                      self.integrated_moments[2]))
            
            
            
    def intmesh2(self, fout):
        
        filepath = os.path.join(self.results_dir, fout)
        self.forces = np.zeros((self.nelements, 3))
        self.moments = np.zeros((self.nelements, 3))
        self.distance = np.zeros((self.nelements, 3))
        self.refpoint = np.zeros(3)
        self.refcs = np.zeros((4, 3))
        self.refcsvec = np.zeros((3, 3))
        self.rotated_forces = np.zeros((self.nelements, 3))
        self.rotated_moments = np.zeros((self.nelements, 3))
        self.integrated_forces = np.zeros(3)
        self.integrated_moments = np.zeros(3)
        
        cspath = os.path.join(self.resources_dir, "cs.txt")

        with open(cspath, 'r') as f1:
            #
            line = f1.readline()
            #
            line = f1.readline()
            temp = line.split()
            self.refpoint[0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()            
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[0, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[1, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[2, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            temp = line.split()
            self.refcs[3, 0:3] = list(map(float, temp[0:3]))
            #
            line = f1.readline()
            #
            line = f1.readline()
            temp = line.split()            
            self.aref = float(temp[0])
            #
            line = f1.readline()
            temp = line.split()            
            self.cref = float(temp[0])
            #
            line = f1.readline()
            temp = line.split()            
            self.bref = float(temp[0])            
                                                
            
        for i in range(3):
            self.refcsvec[0, i] = self.refcs[1, i] - self.refpoint[i]    # Coordinate System Vectors
            self.refcsvec[1, i] = self.refcs[2, i] - self.refpoint[i]
            self.refcsvec[2, i] = self.refcs[3, i] - self.refpoint[i]


        for i in range(self.nelements):
            self.forces[i, 0] = self.cforce[i, 0] # Fx
            self.forces[i, 1] = self.cforce[i, 1] # Fy
            self.forces[i, 2] = self.cforce[i, 2] # Fz
            for j in range(3):
                self.distance[i, j] = self.centers[i, j] - self.refpoint[j]
            
            #a = np.cross(self.forces[i, 0:3], self.distance)
            #print(a)
            self.moments[i, 0:3] = np.cross(self.forces[i, 0:3], self.distance[i, 0:3])  # Mx, My, Mz
            
            #b = np.dot(self.refcsvec, self.forces)
            #print(b)
            self.rotated_forces[i, 0:3] = np.dot(self.refcsvec, self.forces[i, 0:3])
            self.rotated_moments[i, 0:3] = np.dot(self.refcsvec, self.moments[i, 0:3])
        

        self.integrated_forces[0] = np.sum(self.rotated_forces[:, 0]) / self.aref
        self.integrated_forces[1] = np.sum(self.rotated_forces[:, 1]) / self.aref
        self.integrated_forces[2] = np.sum(self.rotated_forces[:, 2]) / self.aref
        
        self.integrated_moments[0] = np.sum(self.rotated_forces[:, 0]) / (self.aref * self.bref)
        self.integrated_moments[1] = np.sum(self.rotated_forces[:, 1]) / (self.aref * self.cref)
        self.integrated_moments[2] = np.sum(self.rotated_forces[:, 2]) / (self.aref * self.bref)
        
        
        with open(filepath, 'w') as f2:    
            f2.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(self.integrated_forces[0], \
                                                                                      self.integrated_forces[1], \
                                                                                      self.integrated_forces[2], \
                                                                                      self.integrated_moments[0], \
                                                                                      self.integrated_moments[1], \
                                                                                      self.integrated_moments[2]))
            
            
            
        
import math
import numpy as np

class FEM2D():
    def __init__(self, file):
        self.input_file = file

    def parse_input(self):
        pass

    def fill_global_dof(self):
        nodes = None
        n_nodes = None
        dof_count = None

        dof_array = np.zeros((n_nodes,3))

        for node_num, node in enumerate(nodes):
            for constraint_num, constraint in enumerate(node['constraints']):
                if constraint == 0:
                    dof_count += 1
                    dof_array[node_num, constraint_num] = dof_count

        return dof_count, dof_array
    
    def get_transformation_matrix(self, elem):
        theta = elem['rotation']

        transformation_matrix = np.zeros((6,6))

        transformation_matrix[1, 1] = math.cos(theta)
        transformation_matrix[1, 2] = math.sin(theta)
        transformation_matrix[2, 1] = -1*math.sin(theta)
        transformation_matrix[2, 2] = math.cos(theta)
        transformation_matrix[3, 3] = 1

        for i_row in range(3):
            for i_col in range(3):
                transformation_matrix[i_row+3, i_col+3] = transformation_matrix[i_row, i_col]

    def get_elem_stiffness_matrix(self, elem):
        A = elem['area']
        I = elem['inertia']
        L = elem['length']
        E = elem['modulus']

        transformation_matrix = self.get_transformation_matrix(elem)
        stiffness_matrix = np.zeros((6,6))

        stiffness_matrix[0] = np.array([A*E/L, 0, 0, -A*E/L, 0, 0])
        stiffness_matrix[1] = np.array([0, 12*E*I/L**3, 6*E*I/L**2, 0, -12*E*I/L**3, 6*E*I/L**2])
        stiffness_matrix[2] = np.array([0, 6*E*I/L**2, 4*E*I/L, 0, -6*E*I/L**2, 2*E*I/L])
        stiffness_matrix[3] = np.array([-A*E/L, 0, 0, A*E/L, 0, 0])
        stiffness_matrix[4] = np.array([0, -12*E*I/L**3, -6*E*I/L**2, 0, +12*E*I/L**3, -6*E*I/L**2])
        stiffness_matrix[5] = np.array([0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, 4*E*I/L])

        stiffness_matrix = np.matmul(np.matmul(np.transpose(transformation_matrix), stiffness_matrix), transformation_matrix)

        return stiffness_matrix
    
    def assemble_global_stiffness(self):
        n_dof, dof_array = self.fill_global_dof()
        nodes = None
        elems = None

        global_stiffness_matrix = np.zeros((n_dof,n_dof))
        for elem in elems:
            local_stiffness_matrix = self.get_elem_stiffness_matrix(elem)
            local_dof1 = -1
            for i_node1 in range(2):
                node_num1 = elem['nodes'][i_node1]
                for i_dof1 in range(3):
                    local_dof2 = local_dof1
                    local_dof1  += 1
                    global_dof1 = dof_array[node_num1, i_dof1]
                    for i_node2 in range(i_node1, 2):
                        node_num2 = elem['nodes'][i_node2]
                        b_dof2 = 1
                        if i_node2 == i_node1: b_dof2 = i_dof1
                        for i_dof2 in range(b_dof2, 3):
                            local_dof2 += 1
                            global_dof2 = dof_array[node_num2, i_dof2]
                            if global_dof1 != 0 and global_dof2 != 0:
                                global_stiffness_matrix[global_dof1, global_dof2] += local_stiffness_matrix[local_dof1, local_dof2]

        for i_row in range(n_dof):
            for i_col in range(i_row+1, n_dof):
                global_stiffness_matrix[i_col, i_row] = global_stiffness_matrix[i_row, i_col]

        return global_stiffness_matrix
    
    def assemble_load_matrix(self):
        n_dof, dof_array = self.fill_global_dof()
        nodes = None
        elems = None

        node_loads = np.zeros((n_dof, 1))
        fixed_end_loads = np.zeros((n_dof, 1))
        prescribed_displacement_loads = np.zeros((n_dof, 1))

        # Assemble the nodal load matrix
        for node_num, node in enumerate(nodes):
            for dof_num in range(len(node['constraints'])):
                if dof_array[node_num, dof_num] != 0:
                    node_loads[dof_array[node_num, dof_num]] += node['point_loads'][dof_num]

        # Assemble the fixed-end load matrix
        for elem in elems:
            local_dof = -1
            for node_num, node in enumerate(elem['nodes']):
                current_node = node['node_number']
                for dof_num in range(len(node['contraints'])):
                    local_dof += 1
                    if dof_array[current_node, dof_num] != 0:
                        fixed_end_loads[current_node, dof_num] += elem['end_loads'][local_dof]

        
        # Assemble the prescribed displacement load matrix
        for elem_num, elem in enumerate(elems):
            delta = np.zeros((6,1))
            for dof_num, displacement in enumerate(elem['nodes']['displacements']):
                delta[dof_num, 0] = displacement
            local_stiffness_matrix = self.get_elem_stiffness_matrix(elem)
            elem_displacement_loads = np.transpose(np.matmul(local_stiffness_matrix, delta))

            local_dof = -1
            for node_num, node in enumerate(elem['nodes']):
                current_node = node['node_number']
                for dof_num in range(len(node['constraints'])):
                    local_dof += 1
                    if dof_array[current_node, dof_num] != 0:
                        prescribed_displacement_loads[current_node, dof_num] += elem_displacement_loads[elem_num, local_dof]

        return node_loads, fixed_end_loads, prescribed_displacement_loads


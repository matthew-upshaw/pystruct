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

        dof_num = np.zeros((n_nodes,2))

        for node_num, node in enumerate(nodes):
            for constraint_num, constraint in enumerate(node['constraints']):
                if constraint == 0:
                    dof_count += 1
                    dof_num[node_num, constraint_num] = dof_count

        return dof_count, dof_num

    def get_elem_stiffness_matrix(self, elem):
        A = elem['area']
        I = elem['inertia']
        L = elem['length']
        E = elem['modulus']

        stiffness_matrix = np.zeros((6,6))

        stiffness_matrix[0] = np.array([A*E/L, 0, 0, -A*E/L, 0, 0])
        stiffness_matrix[1] = np.array([0, 12*E*I/L**3, 6*E*I/L**2, 0, -12*E*I/L**3, 6*E*I/L**2])
        stiffness_matrix[2] = np.array([0, 6*E*I/L**2, 4*E*I/L, 0, -6*E*I/L**2, 2*E*I/L])
        stiffness_matrix[3] = np.array([-A*E/L, 0, 0, A*E/L, 0, 0])
        stiffness_matrix[4] = np.array([0, -12*E*I/L**3, -6*E*I/L**2, 0, +12*E*I/L**3, -6*E*I/L**2])
        stiffness_matrix[5] = np.array([0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, 4*E*I/L])

        return stiffness_matrix
    
    def assemble_global_stiffness(self):
        n_dof = None
        nodes = None
        elems = None

        global_stiffness_matrix = np.zeros((n_dof,n_dof))
        for elem in elems:
            K = self.get_elem_stiffness_matrix(elem)
            for node_num, node in enumerate(elem['nodes']):
                pass

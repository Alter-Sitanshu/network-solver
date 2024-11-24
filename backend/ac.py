import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# Define frequency of AC analysis
FREQUENCY = 100  # in Hz
OMEGA = 2 * np.pi * FREQUENCY  # Angular frequency in rad/s

class Components:
    def __init__(self, u, v, resistance=0, current_source=0, voltage_source=0, capacitance=0, inductance=0):
        self.edge = (u, v)
        self.resistance = resistance
        self.current_source = current_source
        self.voltage_source = voltage_source
        self.capacitance = capacitance
        self.inductance = inductance

    def impedance(self):
        """ Calculate the AC impedance of the component at the given frequency. """
        if self.resistance:
            return self.resistance  # Resistor impedance (real)
        elif self.capacitance:
            return 1 / (1j * OMEGA * self.capacitance)  # Capacitor impedance
        elif self.inductance:
            return 1j * OMEGA * self.inductance  # Inductor impedance
        else:
            return 0  # Ideal source or short-circuit

def main():
    # Define AC components with impedance calculations
    r1 = Components(1, 2, resistance=3000)
    c1 = Components(2, 3, capacitance=3e-5)  # Capacitance in Farads
    l1 = Components(3, 4, inductance=3e-2)   # Inductance in Henries
    v1 = Components(1, 4, voltage_source=-50)  # Voltage source amplitude

    components_list = [r1, c1, l1, v1]

    # Perform AC analysis
    Bf, Cf, Current_source, Voltage_source, Impedance, graph = getTieSetCutSet(components_list)
    Impedance_current, Branch_current = get_impedence_current(Current_source, Voltage_source, Impedance, Bf)
    Branch_voltages, Twig_voltages = get_voltages(Cf, Current_source, Voltage_source, Impedance)

    time_domain_voltages = {edge: phasor_to_time(edge, x) for edge,x in zip(graph.edges, Branch_voltages)}
    time_domain_current = {edge: phasor_to_time(edge, x) for edge,x in zip(graph.edges, Impedance_current)}
    
    time = np.linspace(0, 1 / FREQUENCY, 1000)
    fig, ax = plt.subplots(2, 3, figsize = (12,6))
    ax[0][0].plot(time, time_domain_voltages[(1,2)])
    ax[0][0].grid(True)
    ax[0][1].plot(time, time_domain_voltages[(2,3)])
    ax[0][1].grid(True)
    ax[0][2].plot(time, time_domain_voltages[(3,4)])
    ax[0][2].grid(True)
    ax[1][0].plot(time, time_domain_current[(1,2)])
    ax[1][0].grid(True)
    ax[1][1].plot(time, time_domain_current[(2,3)])
    ax[1][1].grid(True)
    ax[1][2].plot(time, time_domain_current[(3,4)])
    ax[1][2].grid(True)

    plt.show()

def getTieSetCutSet(components_list):
    node_list = [comp.edge for comp in components_list]
    Impedance = []
    Current_source = []
    Voltage_source = []

    graph = nx.DiGraph()
    graph.add_edges_from(node_list)

    for edge in graph.edges:
        for comp in components_list:
            if comp.edge == edge:
                Impedance.append(comp.impedance())
                Current_source.append(comp.current_source)
                Voltage_source.append(comp.voltage_source)

    Impedance = np.array(Impedance, dtype=complex)
    Current_source = np.array(Current_source, dtype=complex)
    Voltage_source = np.array(Voltage_source, dtype=complex)

    # Define tree and cotree for the graph
    non_directed_graph = graph.to_undirected()
    tree_graph = nx.minimum_spanning_tree(non_directed_graph)

    directed_tree_graph = nx.DiGraph()
    directed_tree_graph.add_edges_from(tree_graph.edges)
    cotree_edges = graph.edges - directed_tree_graph.edges()

    cotree_graph = nx.DiGraph()
    cotree_graph.add_edges_from(cotree_edges)

    Bf = np.zeros((len(cotree_edges), len(graph.edges)), dtype=int)
    Cf = np.zeros((len(tree_graph.edges()), len(graph.edges)), dtype=int)
    edge_index = {edge: i for i, edge in enumerate(graph.edges)}

    # Fill Bf and Cf matrices
    for i, edge in enumerate(cotree_graph.edges):
        tree_copy = directed_tree_graph.copy()
        tree_copy.add_edge(*edge)
        cycle = list(nx.find_cycle(tree_copy, orientation="ignore"))
        for start, end, _ in cycle:
            sign = 1 if (start, end) == edge else -1
            Bf[i, edge_index[(start, end)]] = sign

    for i, edge in enumerate(directed_tree_graph.edges()):
        T_minus_edge = directed_tree_graph.copy()
        T_minus_edge.remove_edge(*edge)
        components = list(nx.connected_components(T_minus_edge.to_undirected()))
        component_map = {node: idx for idx, component in enumerate(components) for node in component}
        for start, end in graph.edges:
            sign = 1 if component_map[start] < component_map[end] else -1
            Cf[i, edge_index[(start, end)]] = sign

    return Bf, Cf, Current_source, Voltage_source, Impedance, graph

def get_impedence_current(Current_source, Voltage_source, Impedance, Bf):
    Is = Current_source.T
    Vs = Voltage_source.T
    R = np.diag(Impedance)

    B = Bf
    Z = np.linalg.inv(B @ R @ B.T)
    Loop_current = -Z @ B @ R @ Is - Z @ B @ Vs
    Branch_current = B.T @ Loop_current
    Impedence_current = Branch_current + Is

    return Impedence_current, Branch_current

def get_voltages(Cf, Current_source, Voltage_source, Impedance):
    Y = np.diag(Impedance)
    I = Current_source.T
    V = Voltage_source.T

    for i in range(len(Impedance)):
        if Y[i, i] != 0:
            Y[i, i] = 1.0 / Y[i, i]
        else:
            Y[i, i] = 1e12  # Large impedance for open circuit

    C = Cf
    P = np.linalg.inv(C @ Y @ C.T)
    Twig_voltages = P @ C @ I + P @ C @ Y @ V
    Branch_voltages = C.T @ Twig_voltages

    return Branch_voltages, Twig_voltages

def phasor_to_time(edge, Branch_values):
    I_magnitude = np.abs(Branch_values)
    if abs(Branch_values.imag) < 1e-10:
        Branch_values = I_magnitude + 0j # making the imaginary part zero
    I_phase = np.angle(Branch_values)
    time = np.linspace(0, 1 / FREQUENCY, 1000)
    I_time = I_magnitude * np.sin(OMEGA * time + I_phase)
    print(edge, " : ", I_phase)

    return I_time
if __name__ == "__main__":
    main()

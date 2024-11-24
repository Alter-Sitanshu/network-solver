# %%
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# %%
class Components:
    def __init__(self,u, v, resistance=0 , current_source=0, voltage_source=0, capacitance = 0, inductance = 0):
        self.edge = (u,v)
        self.resistance = resistance 
        self.current_source = current_source
        self.voltage_source = voltage_source
        self.capacitance = capacitance
        self.inductance = inductance
    def pt(self):
        return f"current_source={self.current_source},voltage_source={self.voltage_source})"

def main(circuit,simulationType,simulationTimeStep,stamps):
    for comp in circuit:
        print("Component:",comp,"\n")
    print("simulationType:",simulationType)
    print("simulationTimeStep:",simulationTimeStep,"\n")

    # Get components list
    components_list = get_components_list(circuit)

    # Q = get_transient(components_list,h=simulationTimeStep,stamps=stamps,simulationType=simulationType)
    # result = solve(components_list,h=simulationTimeStep,stamps=stamps)
    # result = convert_res_to_json(result)
    # return result

    if simulationType=="dc":
        dc_operating_point_response,graph = get_dc_operating_point_response(components_list, simulationType)
        result = convert_res_to_json(dc_operating_point_response)
        return result
    elif simulationType=="transient":
         prev_responses,edge_list = get_transient_response(components_list,h=simulationTimeStep,stamps=stamps,simulationType=simulationType)
         result = convert_res_to_json(prev_responses)
         return result

def convert_res_to_json(response):

    json_compatible_result = {
    str(key): [[float(value) for value in sublist] for sublist in values]
    for key, values in response.items()
    }
    
    return json_compatible_result

def get_components_list(circuit):
    components_list = []

    for component in circuit:
        if component["category"]=="resistor":
            comp = Components(u = component["node1"], v = component["node2"], resistance=float(component["value"]))
            components_list.append(comp)
        elif component["category"]=="capacitor":
            comp = Components(u = component["node1"], v = component["node2"], capacitance=float(component["value"]))
            components_list.append(comp)
        elif component["category"]=="inductor":
            comp = Components(u = component["node1"], v = component["node2"], inductance=float(component["value"]))
            components_list.append(comp)
        elif component["category"]=="voltageSource":
            comp = Components(u = component["node1"], v = component["node2"], voltage_source=float(component["value"]))
            components_list.append(comp)
        elif component["category"]=="currentSource":
            comp = Components(u = component["node1"], v = component["node2"], current_source=float(component["value"]))
            components_list.append(comp)

    return components_list

def get_dc_operating_point_response(components_list,simulationType):
    # If simulation type is dc
    if simulationType=="dc":
        for component in components_list:
            if component.resistance==0 and component.capacitance!=0:
                component.resistance=float(1000000000)
            elif component.resistance==0 and component.capacitance==0 and component.inductance!=0:
                component.resistance = float(0)
    else:
        # To get the very first response
        for component in components_list:
            if component.resistance==0 and component.capacitance!=0:
                component.resistance=float(0)
            elif component.resistance==0 and component.capacitance==0 and component.inductance!=0:
                component.resistance = float(1000000000)

    Bf,Cf,Current_source,Voltage_source,Resistance,graph = getTieSetCutSet(components_list)
    edge_list = graph.edges
    dc_operating_point_response: dict[tuple, list] = {edge: [] for edge in edge_list}

    Impedence_current,Branch_current = get_impedence_current(Current_source,Voltage_source,Resistance,Bf)
    Branch_voltages,Twig_voltages= get_voltages(Cf,Current_source,Voltage_source,Resistance)

    for index, edge in enumerate(edge_list):
        dc_operating_point_response[edge].append([Branch_voltages[index], Branch_current[index]])

    return dc_operating_point_response,graph

def get_transient_response(components_list,h,stamps=15,simulationType="transient"):
    # store previous responses -> BranchCurrent and Branch Voltages

    # To get the very first response
    first_response,graph = get_dc_operating_point_response(components_list,simulationType)
    edge_list = graph.edges
    prev_responses: dict[tuple, list] = {edge: [] for edge in edge_list}

    prev_responses = first_response

#  CHANGING THE CAPACITANCE AND INDUCTORS INTO RESISTANCES
    capacitance: list[Components] = []
    inductors: list[Components] = []
    for component in components_list:
        if component.resistance == 0 and component.capacitance!=0:
            component.resistance= float(h)/float(component.capacitance)
            capacitance.append(component)

        if component.inductance!=0:
            component.resistance = float(component.inductance)/float(h)
            inductors.append(component)

    for i in range(stamps):
        for c in capacitance:
            Vn = float(prev_responses[c.edge][-1][0])
            model_current_source = float(Vn)/float(c.resistance)
            c.current_source = model_current_source
        for l in inductors:
            In = float(prev_responses[l.edge][-1][1])
            model_volt = float(In) * float(l.resistance)
            l.voltage_source = -model_volt

        for comp in components_list:
            print(comp.pt())
            
        
        cir_response,graph = get_dc_operating_point_response(components_list,simulationType)
        edge_list = graph.edges

        for index, edge in enumerate(edge_list):
            prev_responses[edge].append(cir_response[edge][0])

    return prev_responses,edge_list

def solve(components_list,h=0.1,stamps=20):
    max_node_number = 0
    voltage_source_idx = 0

    for component in components_list:
        if component.edge[0]>max_node_number:
            max_node_number = component.edge[0]
        elif component.edge[1]>max_node_number:
            max_node_number = component.edge[1]

    edge_list = []
    for component in components_list:
        edge_list.append(component.edge)

    prev_responses: dict[tuple, list] = {edge: [] for edge in edge_list}

    # first response
    for component in components_list:
            if component.resistance==0 and component.capacitance!=0:
                component.resistance = float(0.0000000001)
            elif component.resistance==0 and component.capacitance==0 and component.inductance!=0:
                component.resistance=float(100000000)

    first_response,S = get_transient(components_list=components_list)
    print(first_response)

    # Save the first response
    for component in components_list:
        edge = component.edge
        res = first_response.T[0][edge[1]-1]-first_response.T[0][edge[0]-1]
        res_I = res/(component.resistance) if (component.resistance != 0) else 0
        prev_responses[edge].append((res,res_I))

    capacitance: list[Components] = []
    inductors: list[Components] = []
    for component in components_list:
        if component.capacitance!=0:
            component.resistance= float(h)/float(component.capacitance)
            capacitance.append(component)

        if component.inductance!=0:
            component.resistance = float(component.inductance)/float(h)
            inductors.append(component)

    for i in range(stamps):
        voltage_source_idx = 0
        for c in capacitance:
            Vn = float(prev_responses[c.edge][-1][0])
            model_current_source = float(Vn)/float(c.resistance)
            c.current_source = model_current_source

        for l in inductors:
            In = float(prev_responses[l.edge][-1][1])
            l.current_source = -In

        response,S = get_transient(components_list=components_list)
        # print(response)
        for component in components_list:
            edge = component.edge
            res = response.T[0][edge[1]-1]-response.T[0][edge[0]-1]
            if component.resistance != 0:
                res_I = res/(component.resistance)
            elif component.resistance == 0 and component.voltage_source != 0:
                res_I = response.T[0][max_node_number + voltage_source_idx-1]
                voltage_source_idx += 1
            elif component.current_source != 0:
                res_I = component.current_source

            prev_responses[edge].append((res,res_I))

    # for edge in edge_list:
    #     print(edge,": ",prev_responses[edge],"\n")

    return prev_responses

def get_transient(components_list):

    max_node_number = 0
    voltage_source_count = 0
    gnd_node = 1

    for component in components_list:
        if component.voltage_source!=0:
            voltage_source_count = voltage_source_count+1

        if component.edge[0]>max_node_number:
            max_node_number = component.edge[0]
        elif component.edge[1]>max_node_number:
            max_node_number = component.edge[1]
        
    G = np.zeros((max_node_number+voltage_source_count,max_node_number+voltage_source_count))
    S = np.zeros((1,max_node_number+voltage_source_count))

    # For gnd_node
    G[0][0] = 1

    voltage_source_idx = 0

    for component in components_list:
        u = component.edge[0]
        v = component.edge[1]

        if component.resistance!=0:
            g = float(1/component.resistance)

            if u== gnd_node:
                G[v-1][v-1] += g
            else:
                G[u-1][u-1] += g
                G[u-1][v-1] += -g
                G[v-1][v-1] += g
                G[v-1][u-1] += -g

        if component.voltage_source!=0:
            if u==gnd_node:
                G[max_node_number+voltage_source_idx][v-1]= 1
                G[v-1][max_node_number+voltage_source_idx]= 1
            else:
                G[max_node_number+voltage_source_idx][u-1]= 1
                G[u-1][max_node_number+voltage_source_idx]= 1
                G[max_node_number+voltage_source_idx][v-1]= -1
                G[v-1][max_node_number+voltage_source_idx]= -1

            S[0][max_node_number+voltage_source_idx] = component.voltage_source

            voltage_source_idx = voltage_source_idx+1

        if component.current_source!=0:
            if u==gnd_node:
                S[0][v-1] = component.current_source
            else:
                S[0][u-1] = component.current_source
                S[0][v-1] = -component.current_source
    
    # print(G)
    # print(S)
    Q = np.dot(np.linalg.inv(G),S.T)
    # print(Q)

    return Q,S

def getTieSetCutSet(components_list):

    node_list = [comp.edge for comp in components_list]
    Resistance = []
    Current_source = []
    Voltage_source = []

    graph = nx.DiGraph()
    graph.add_edges_from(node_list)
    pos = nx.planar_layout(graph)

    # %%

    for edge in graph.edges:
        for comp in components_list:
            if comp.edge == edge:
                Resistance.append(comp.resistance)
                Current_source.append(comp.current_source)
                Voltage_source.append(comp.voltage_source)
                
    Resistance = np.array(Resistance)
    Current_source = np.array(Current_source)
    Voltage_source = np.array(Voltage_source)

    non_directed_graph = graph.to_undirected()
    tree_graph = nx.minimum_spanning_tree(non_directed_graph)

    directed_tree_graph = nx.DiGraph()
    directed_tree_graph.add_edges_from(tree_graph.edges)
    cotree_edges = graph.edges() - directed_tree_graph.edges()

    cotree_graph = nx.DiGraph()
    cotree_graph.add_edges_from(cotree_edges)

    Bf = np.zeros((len(cotree_edges), len(graph.edges())), dtype=int)
    Cf = np.zeros((len(tree_graph.edges()), len(graph.edges())), dtype=int)
    edge_index = {edge: i for i, edge in enumerate(graph.edges())}


    for i,edge in enumerate(cotree_graph.edges):
        tree_copy = directed_tree_graph.copy()
        tree_copy.add_edge(*edge)
        cycle = list(nx.find_cycle(tree_copy, orientation="ignore"))

        co_tree_orient = [i[2] for i in cycle if (i[0],i[1]) == edge][0]
        for (start, end, orientation) in cycle:
            sign = 1 if orientation == co_tree_orient else -1
            Bf[i, edge_index[(start, end)]] = sign

    for i, edge in enumerate(directed_tree_graph.edges()):
        T_minus_edge = directed_tree_graph.copy()
        
        T_minus_edge.remove_edge(*edge)
        T_undirected = T_minus_edge.to_undirected().copy()

        components = list(nx.connected_components(T_undirected))
        component_map = {node: idx for idx, component in enumerate(components) for node in component}
        for start, end in graph.edges():
            if component_map[start] > component_map[end]:
                sign = -1
                if start < end:
                    Cf[i, edge_index[(start, end)]] = sign
            elif component_map[start] < component_map[end]:
                sign = 1
                if start < end:
                    Cf[i, edge_index[(start, end)]] = sign

    return Bf,Cf,Current_source,Voltage_source,Resistance,graph

def get_impedence_current(Current_source,Voltage_source,Resistance,Bf):
    Is = Current_source.T
    Vs = Voltage_source.T
    R = np.diag(Resistance)

    B = Bf
    B=B.astype(float);
    R=R.astype(float);
    Is=Is.astype(float);
    Vs=Vs.astype(float);
    Z=np.linalg.inv(np.dot(B,np.dot(R,B.T)))
    Loop_current = -np.dot(Z,np.dot(B,np.dot(R,Is)))-np.dot(Z,np.dot(B,Vs))
    Branch_current=np.dot(B.T,Loop_current)
    Impedence_current = Branch_current + Is

    return Impedence_current,Branch_current

def get_voltages(Cf,Current_source,Voltage_source,Resistance):
    Y= np.diag(Resistance)
    Y= Y.astype(float)
    I=Current_source.T
    I= I.astype(float)
    V= Voltage_source.T
    V=V.astype(float)

    for i in range(0,len(Resistance)):
        if Y[i][i]!=0:
            Y[i][i]=(1.0/Y[i][i])
        else:
            Y[i][i]= float(1000000000)
    C=Cf
    C=C.astype(float)

    P=np.linalg.inv(np.dot(C,np.dot(Y,C.T)))
    Twig_voltages = np.dot(P,np.dot(C,I))+ np.dot(P,np.dot(C,np.dot(Y,V)))
    Branch_voltages = np.dot(C.T,Twig_voltages)

    return Branch_voltages,Twig_voltages

if __name__ == "__main__":
    main()




"""The present code provides a framework to run pseudoflow for ultimate pit limit definition.
It has been developed by S. Avalos from the Geometallurygical Group at Queen's University as part of a PhD program.
The code is not free of bugs but running end-to-end. 
Any comments and further improvements are well recevied at: 17saa6@queensu.ca
July 16th, 2020.
Geomet Group - Queen's University - Canada"""


import numpy as np
import networkx as NetX
import pseudoflow as pf
import time

def Pseudoflow_UPL(BM, nx, ny, nz, VarIn, VarOut):
    print("Beginning Pseudoflow")
    start_UPL = time.time() 
    source = 0
    sink = np.int(nx*ny*nz + 1)
    
    # Graph creation
    Graph = NetX.DiGraph()
    
    # External arcs creation by external function. Source - Nodes, Nodes - Sink
    Graph = CreateExternalArcs(BM, nx, ny, nz, Graph=Graph, Var=VarIn)
    
    # Internal arcs creation by external function. Block precedence (1x5 or 1x9)
    for ind_z in range(nz - 1):
        pos_z = nz - ind_z - 2
        for pos_y in range(ind_z + 1, ny - ind_z - 1):
            for pos_x in range(ind_z + 1, nx - ind_z - 1):
                # Precedence of 5 blocks
                Graph = CreateInternalArcs1x5(pos_x, pos_y, pos_z, nx, ny, Graph=Graph)
                # Precedence of 9 blocks
                #Graph = CreateInternalArcs1x9(pos_x, pos_y, pos_z, nx, ny, Graph=Graph)
    
    # Solving the minimum cut problem via pf.hpf solver
    RangeLambda = [0]
    breakpoints, cuts, info = pf.hpf(Graph, source, sink, const_cap="const", mult_cap="mult", lambdaRange=RangeLambda, roundNegativeCapacity=False)
    
    #Going over the cuts.items finding the nodes inside the resulting UPL.
    B = {x:y for x, y in cuts.items() if y == [1] and x!=0}
    InsideList = list(B.keys())
    
    # Set all blocks as zero
    BM[:,VarOut] = 0 

    for indUPL in range(len(InsideList)): 
        # Set blocks inside UPL as one
        BM[np.int(InsideList[indUPL] -1),VarOut] = 1
    
    print("--> Pseudoflow time: --%s seconds " % (np.around((time.time() - start_UPL), decimals=2)))  

    return BM

def CreateExternalArcs(BM, nx, ny, nz, Graph, Var):
    Sink = np.int(nx*ny*nz + 1)

    for t_z in range(nz):
        pos_z = nz - t_z - 1
        for t_y in range(t_z, ny-t_z):
            for t_x in range(t_z,nx-t_z):
                p_i = 1 + t_x + nx*t_y + ny*nx*pos_z 
                Capacity = np.absolute(np.around(BM[p_i-1,Var], decimals=2))
                if BM[p_i-1,Var] < 0: #Negative local Economic Value
                    Graph.add_edge(p_i, Sink, const=Capacity, mult=1)
                else:
                    Graph.add_edge(0, p_i, const=Capacity, mult=1)
    return Graph

def CreateInternalArcs1x9(pos_x, pos_y, pos_z, nx, ny, Graph):
    p_0 =  1 + pos_x + nx*pos_y + ny*nx*pos_z
    
    p_1 =  1 + (pos_x-1) + nx*(pos_y-1) + ny*nx*(pos_z+1)
    p_2 =  1 + (pos_x) + nx*(pos_y-1) + ny*nx*(pos_z+1)
    p_3 =  1 + (pos_x+1) + nx*(pos_y-1) + ny*nx*(pos_z+1)
    p_4 =  1 + (pos_x-1) + nx*(pos_y) + ny*nx*(pos_z+1)
    p_5 =  1 + (pos_x) + nx*(pos_y) + ny*nx*(pos_z+1)
    p_6 =  1 + (pos_x+1) + nx*(pos_y) + ny*nx*(pos_z+1)
    p_7 =  1 + (pos_x-1) + nx*(pos_y+1) + ny*nx*(pos_z+1)
    p_8 =  1 + (pos_x) + nx*(pos_y+1) + ny*nx*(pos_z+1)
    p_9 =  1 + (pos_x+1) + nx*(pos_y+1) + ny*nx*(pos_z+1)
    
    Graph.add_edge(p_0, p_1, const=99e99, mult=1)
    Graph.add_edge(p_0, p_2, const=99e99, mult=1)
    Graph.add_edge(p_0, p_3, const=99e99, mult=1)
    Graph.add_edge(p_0, p_4, const=99e99, mult=1)
    Graph.add_edge(p_0, p_5, const=99e99, mult=1)
    Graph.add_edge(p_0, p_6, const=99e99, mult=1)
    Graph.add_edge(p_0, p_7, const=99e99, mult=1)
    Graph.add_edge(p_0, p_8, const=99e99, mult=1)
    Graph.add_edge(p_0, p_9, const=99e99, mult=1)
    
    return Graph

def CreateInternalArcs1x5(pos_x, pos_y, pos_z, nx, ny, Graph):
    p_0 =  1 + pos_x + nx*pos_y + ny*nx*pos_z
    p_2 =  1 + (pos_x) + nx*(pos_y-1) + ny*nx*(pos_z+1)
    p_4 =  1 + (pos_x-1) + nx*(pos_y) + ny*nx*(pos_z+1)
    p_5 =  1 + (pos_x) + nx*(pos_y) + ny*nx*(pos_z+1)
    p_6 =  1 + (pos_x+1) + nx*(pos_y) + ny*nx*(pos_z+1)
    p_8 =  1 + (pos_x) + nx*(pos_y+1) + ny*nx*(pos_z+1)

    Graph.add_edge(p_0, p_2, const=99e99, mult=1)
    Graph.add_edge(p_0, p_4, const=99e99, mult=1)
    Graph.add_edge(p_0, p_5, const=99e99, mult=1)
    Graph.add_edge(p_0, p_6, const=99e99, mult=1)
    Graph.add_edge(p_0, p_8, const=99e99, mult=1)

    return Graph

def main():    
    print("Start")
    start_time = time.time() 
    nx, xmn, xsiz = 44, 24300, 16
    ny, ymn, ysiz = 62, 24800, 16
    nz, zmn, zsiz = 26, 3600, 16

    BlockModel = np.loadtxt("BM_matrix.txt") # Import Block Model
    BlockModel = Pseudoflow_UPL(BM=BlockModel, nx=nx, ny=ny, nz=nz, VarIn=4, VarOut=5)
    
    '''Save Block Model'''
    np.savetxt(fname="BM_matrix.txt", X=BlockModel, fmt='%.3f', delimiter='\t')	
          
    return print("--%s seconds of the whole process-" % (np.around((time.time() - start_time), decimals=2)))  


if __name__ == "__main__":
    main()
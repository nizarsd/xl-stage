######################### Graph Generator v1.1 ##################
#																#
# 				   		 University of York					    #
# 				          Graceful Project					    #
# 						Created by Nizar Dahir				    #
# 						    Baghdad, 2020					    #
#																#
#################################################################
# This script does the same job of graph_gen.m script 
# In addition to generating the TG text files
# TG files include TG edgees formated at "Src. Dest. Bandwidth"
# Now implemented in Python
# Nizar Dahir 9/4/2020

# This code expects:-
#1- A dirctory named "graphs" at the same dirctory at which the code is running
#2- Graphviz installed an added to  PATH    

import numpy as np 
import random as rn 
from graphviz import Source
import os

# change this to the graphviz path in your computer (not necessary of gaphviz is already in PATH)
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'


def onesidenorm(mu,sig):
    y = 0	
    while (y < mu):
        y =  np.random.normal(mu,sig);

    return y
    
    
def limitednorm(mu,sig,mn,mx):	
    y = 0
    while (y <= mn or y >= mx):
        y =  np.random.normal(mu,sig);

    return y 
     
# Set the seed to repeat results
nodes = 26

# Probability of connection
prob=0.1

# Number of graphs		
N = 3


label_fsize = 20

# min_actual_ranks = 3
# minrank = 3
# maxrank = 3

# Egde properties
MAX_BW = 500	# max bandwdth
MIN_BW = 0.1    # min. bw
MEAN_BW = 50    # mean of bw
VAR_BW = 50     # variance of bw

# Rank properiteis 
min_actual_ranks = 7 #int(nodes**0.5) 
minrank = min_actual_ranks # int(min_actual_ranks*1.2) 
maxrank = min_actual_ranks # int(min_actual_ranks*1.5)

rn.seed(10)
np.random.seed(10)

for g in range(0,N):
    valid = False;
    while (not valid):
        conn = np.zeros([nodes,nodes])

        ranks = rn.randint(minrank, maxrank)
        
        # mean of node binning
        # mu_rank =  0     # Convergent graphs 
        # mu_rank =  ranks # Divergent graphs
        mu_rank =  ranks/2 # Both 

        # variance of node binning
        sig_rank =  ranks/4
        
        nranks = np.zeros(nodes)
        while True:
            for i in range(1, nodes):
                nranks[i] =  int(limitednorm(mu_rank,sig_rank, 0, ranks)+0.5)
    #		    nranks(i) =  randi([1 ranks], 1);
            
            nranks=  np.sort(nranks)	

            if (len(np.unique(nranks)) >= min_actual_ranks): 
                break

        d_probability = 2 # determines how fast the probability drops with rank distancing
        for s in range(nodes):
            for d in range(s + 1, nodes):
                conn[s,d] =  (nranks[d] - nranks[s])/d_probability 
        # print(conn)
      
        # set connection probabilities 
        for s in range(nodes -1):
            for d in range(s + 1, nodes):
                if (conn[s,d] > 0):
                    conn[s,d] = pow(prob,conn[s,d])
                    conn[s,d] = rn.random() < conn[s,d] 
        
        conn = np.triu(conn)
                
        # set connected to be one with probability of conn
        # print(conn)
        
        # nonempty ranks
        aranks=np.unique(nranks);

        ideg = np.sum(conn, axis = 0); 
        odeg = np.sum(conn, axis = 1);
        
        ####################################################    
        # Repair unconnected nodes
        ####################################################    
        
        for i in range(1, nodes):
            rnode = nranks[i]
            Fi = np.where(rnode == aranks)[0]
            # for nodes that have in degree of zero and not in the first rank 
            if (not ideg[i]) and  (nranks[i] > aranks[0]):	
                Fi = np.where(aranks==rnode)[0]
                prev_rank = aranks[Fi - 1]
                prev_nodes = np.where(nranks == prev_rank)[0]
                Ai = np.argmin(odeg[prev_nodes])	
                conn[prev_nodes[Ai],i] = 1
            
            # for nodes that have out degree of zero and not in the last rank 
            if (not odeg[i] and  rnode < aranks[-1]):			
                Fi = np.where(aranks == rnode)[0]
                next_rank = aranks[Fi + 1]
                next_nodes = np.where(nranks == next_rank)[0]
                Ai = np.argmin(ideg[next_nodes])	
                conn[i, next_nodes[Ai]] = 1
                
        ####################################################    
        # check if the graph has any unconnected subgraphs 
        ####################################################   
        
        n=[]
        V=[]
        V0=[]
        V0.append(0)
        # Walk through all nodes starting from node 0
      
        while True:
            for i in range(len(V0)):
                n1 = np.where(conn[V0[i],:] > 0)[0]
                n.extend(n1)
                n2 = np.where(conn[:, V0[i]] > 0)[0]
                n.extend(n2)
            
            Vn = list(set(n) -  set(V))

            V.extend(Vn)
            if len(Vn) == 0: 
                break
            V0=Vn	
        

        if (len(V) == nodes):
            # print (V)
            valid = True
        else:
            print("invalid !");
   
    ####################################################    
    #                   Write GV and TG FILES
    ####################################################   
    # TG file is non standard text file used by us
    
    fname = "graphs/dag%04i" %g
    fgv = open(fname + '.gv', 'w')
    ftg = open(fname + '.tg.txt', 'w')
   
    f_str = "digraph {\n"
    f_str += "  splines=true;\n\r"
    f_str += "node [margin=0 fontname=arial fontcolor=black fontsize=20 shape=circle " 
    f_str += "width=0.9 fixedsize=true style=filled fillcolor=powderblue]\n\r" #%(label_fsize)

    # Nodes
    for i in range(1, nodes):
        # f_str += "  %i [label=\"V%i\"]\n"%(i,i)
        f_str += "  %i [label=\"P%i\"]\n"%(i,i)
    

    f_str += "rankdir=LR\n\r"
    f_str += "edge [margin=0 fontname=arial fontcolor=black fontsize=20]\n\r" #%(label_fsize)

    s, d = np.where(conn==1);

    no_edges = len(s); 

    # Edges 
    # loads = np.random.normal(MIN_BW, MAX_BW, no_edges)
    
    # interger or float loads ?
    # loads = [ limitednorm(MEAN_BW,VAR_BW,MIN_BW,MAX_BW)  for i in range(no_edges)] 
    loads = [ int(limitednorm(MEAN_BW,VAR_BW,MIN_BW,MAX_BW) + 0.5) for i in range(no_edges)] 
    
    
    # print(s)
    tg_str = ''
    for i in range(len(s)):
        # f_str += "	%i -> %i [label=\"%0.2f\"]\n" % (s[i],d[i], loads[i])
        f_str += "	%i -> %i [label=\"%i\"]\n" % (s[i],d[i], loads[i])
        tg_str+= "%i  %i  %i \n" % (s[i],d[i], loads[i])

   
    # Ranks
    for i in aranks:
        same = np.where(nranks==i)[0]

        f_str = f_str + "	{rank=same "
        for n in same:	
            f_str += " %i," % (n)
        
       
        f_str = f_str[:-1] + "}\n"
    
      
    f_str += "} \n\r"
    
    fgv.write(f_str)
    ftg.write(tg_str)
    
    fgv.close()
    ftg.close()
    
    
    
    ####################################################    
    # Write XML FILE
    ####################################################   
    
    # Nodes
    
    fid=open(fname+ '.xml', "w")
    f_str = ''
    f_str += "<mc:Graph>\n"
    f_str += " <mc:NodeList>\n"
    for i in range(1, nodes):
        f_str += "	<mc:Node>\n"
        f_str += "		<mc:id>%i</mc:id>\n"%i
        f_str += "		<mc:name>P%i</mc:name>\n"%i
        f_str += "		<mc:type>p</mc:type>\n"
        f_str += "		<mc:rank>%i</mc:rank>\n"%nranks[i]
    

    f_str += "	</mc:Node> \n\n"
    f_str += " </mc:NodeList>\n"

    # Edges 
    f_str += " <mc:EdgeList>\n"

    for i in range (len(s)): 
        f_str += "	<mc:Edge>\n"
        f_str += "      	<mc:sourceId>%i</mc:sourceId>\n"%s[i]
        f_str += "       	<mc:targetId>%i</mc:targetId>\n"%d[i]
        # f_str += "      	<mc:networkLoad>%0.2f</mc:networkLoad>\n"%loads[i]
        f_str += "      	<mc:networkLoad>%i</mc:networkLoad>\n"%loads[i]
        f_str += "	</mc:Edge>\n"
    
        f_str += " </mc:EdgeList>\n\n"
    # Ranks
        f_str += " <mc:RankList>\n"

    for r in aranks:
        f_str += "	<mc:RankGroup>\n"
        rnodes = np.where(nranks == r)[0]
        for n in rnodes:
            f_str += "		<mc:Rank>%i</mc:Rank>\n"% n
        

        f_str += "	</mc:RankGroup>\n"
    
    f_str += " </mc:RankList>\n\n"
    f_str += "</mc:Graph>\n"

    fid.write(f_str)
    fid.close()

    # Create pdf and view 
    
    s = Source.from_file(fname+'.gv')
    s.view()



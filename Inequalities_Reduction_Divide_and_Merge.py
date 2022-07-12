"""
Parameters 
1 - name of cipher (PRESENT/GIFT/WARP/TWINE/ASCON/FIDES-5/SC2000-5)
2 - sbox or prob   (This research uses sbox)
3 - DDT points (2/4/6), '-' for all  (This research uses '-')
4 - If do not want to introduce new inequalities then use '-' otherwise specify batch size(beta) in which inequalities are to be added e.g.'10' or '50', if batch size is '0' then every inequality will added with every other inequality (specifying '0' is same as len(inequalities), 0 must be used to get results equivalent to Boura and Coggia's approach(https://hal.inria.fr/hal-03046211/document) of inequalities reduction)
5 - k(=2/3/4) inequalities to be added (it should not be 1), if parameter 4 is '-' then use the parameter 5 as '-'

Examples

python Inequalities_Reduction_Divide_and_Merge.py WARP sbox - 10 2
python Inequalities_Reduction_Divide_and_Merge.py WARP sbox - 20 3
    
"""

import sys
import math
import numpy as np
#import gurobipy
#from gurobipy import GRB 
from numpy import array, hstack, ones, vstack
import cdd
from docplex.mp.model import Model

file_name = "_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4] +"_"+sys.argv[5] 

def compute_polytope_halfspaces(vertices):

    V = vstack(vertices)
    t = ones((V.shape[0], 1))  # first column is 1 for vertices
    tV = hstack([t, V])
    mat = cdd.Matrix(tV, number_type='fraction')
    mat.rep_type = cdd.RepType.GENERATOR
    P = cdd.Polyhedron(mat)
    bA = array(P.get_inequalities())
    #return bA
    if bA.shape == (0,):  # bA == []
        return bA
    # the polyhedron is given by b + A x >= 0 where bA = [b|A]
    b, A = array(bA[:, 0]), array(bA[:, 1:])
    
    return (A, b)

def MILP_Solve(ineq_list,impossible_diff_arr):
        
        ########################################GUROBI MODEL##################################################
#        if (IN_COLAB==True):
#          print("environment set")
#          m = gurobipy.Model(env=e)
#        else:
#          m = gurobipy.Model()
#        for i in range(0,len(ineq_list)):
#            m.addVar(vtype=GRB.BINARY,name="z%s" % str(i))
#
#        m.update()
#        variables = m.getVars()
#        
#        for i in range(0,len(impossible_diff_arr)):
#            ineq_solve_count = (np.multiply(np.array(impossible_diff_arr[i]),np.array(ineq_list))).sum(1)
#            less_than_zero = np.where(ineq_solve_count<0)[0]
#            m.addConstr(sum([variables[x] for x in less_than_zero]) >=1)
#        m.setObjective(sum(m.getVars()),GRB.MINIMIZE)
#        m.update()
#        end_3 = time.process_time()
#        print(">>Time>> TTime Taken to write lp model: " + str(end_3 - start_3))
#        m.write("GUROBI_Problem"+ file_name +".lp")
#        m.optimize()
#        m.write("GUROBI_Solution"+ file_name +".sol")
#        m.setParam(GRB.Param.SolutionNumber, 1)
#        sol_arr = np.array(m.X)
#        final_sol = np.where(sol_arr == 1)[0].tolist()
        
        #########################################################################################################
        ###################################MODEL CPLEX########################################################
        m = Model()
        m.binary_var
        for i in range(0,len(ineq_list)):
            m.binary_var(name="z%s" % str(i))

        variables = [i for i in m.iter_binary_vars()]
        for i in range(0,len(impossible_diff_arr)):
            ineq_solve_count = (np.multiply(np.array(impossible_diff_arr[i]),np.array(ineq_list))).sum(1)
            less_than_zero = np.where(ineq_solve_count<0)[0]
            m.add_constraint(sum([variables[x] for x in less_than_zero]) >=1)
        m.set_objective("min",sum(variables))
        m.export("CPLEX_Problem"+file_name +".lp")

        sol = m.solve(log_output=True)
        sol.export("CPLEX_Solution" +file_name+".sol")
        print(sol.get_objective_value())
        sol_arr = np.array(sol.get_values(variables)) #[ sol for var in variables ]
        final_sol = np.where(sol_arr == 1)[0].tolist()
        #########################################################################################################
        
        f = open("ineq"+ file_name +".txt","w")
        ineq_list = ineq_list.astype(int)
        #count = 1
        for ineqlitity in final_sol:
            #print(ineq_list[ineqlitity])
            inequality_rotated = ineq_list[ineqlitity].tolist()[1:] + [ineq_list[ineqlitity].tolist()[0]] 
            print(str(inequality_rotated)[1:-1]+",")
            f.write(str(inequality_rotated).replace("[","").replace("]",",\n"))    
            
            ############# To print in LaTex#######################################################################
#            print_str = str(count)+"&$\\;\\;\\;"
#            count += 1
#            for i in range(0,len(inequality_rotated)-1):
#                if (inequality_rotated[i]>=0):
#                    print_str += "+"
#                if( i <= 3):
#                    print_str += str(inequality_rotated[i]) + "*x_" + str(3-i) + "  "
#                if( i >= 4 and  i <= 7):
#                    print_str += str(inequality_rotated[i]) + "*y_" + str(3-(i%4)) + "  "
#                if( i >= 8):
#                    print_str += str(inequality_rotated[i]) + "*p_" + str(1-(i%2)) + "  "
#            print_str +=  " \\geq -" + str(str(inequality_rotated[-1])) + "$ \\\\"
#            print(print_str)
            ###########################################################################################################
            
        f.close()
        #########################################To print index of reduced set of inequalities###########################     
#        ineq_list_mod = []
#        print(final_sol)
#        for x in final_sol:
#            ineq_list_mod.append(ineq_list[x])
        #########################################################################################################
        return final_sol
        
    
def print_DDT(table):
    for row in range(len(table)):
        for col in range(len(table[row])):
            print(table[row][col],end='');
            if col == len(table[row])-1:
                print("\n");

#############################################List of S-boxes###################################################
if (sys.argv[1] == "PRESENT"):
  s_box = ((0xC , 0x5 , 0x6 , 0xB , 0x9 , 0x0 , 0xA , 0xD , 0x3 , 0xE , 0xF , 0x8 , 0x4 , 0x7 , 0x1, 0x2),);     # PRESENT
if (sys.argv[1] == "WARP"):
  s_box = ((0xc,0xa, 0xd,0x3, 0xe,0xb, 0xf, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6),);     # WARP
elif (sys.argv[1] == "GIFT"):
  s_box = ((0x1,0xa, 0x4,0xc, 0x6,0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe),);       # GIFT       
elif (sys.argv[1] == "TWINE"):
  s_box = ((0xc,0x0, 0xf,0xa, 0x2,0xb, 0x9, 0x5, 0x8, 0x3, 0xd, 0x7, 0x1, 0xe, 0x6, 0x4),);       # TWINE         
elif (sys.argv[1] == "ASCON"):
  s_box = ((0x4,0xb,0x1f,0x14,0x1a,0x15,0x9,0x2,0x1b,0x5,0x8,0x12,0x1d,0x3,0x6,0x1c,), (0x1e,0x13,0x7,0xe,0x0,0xd,0x11,0x18,0x10,0xc,0x1,0x19,0x16,0xa,0xf,0x17));       # ASCON  
elif (sys.argv[1] == "FIDES-5"):
  s_box = ((1,0,25,26,17,29,21,27,20,5,4,23,14,18,2,28),(15,8,6,3,13,7,24,16,30,9,31,10,22,12,11,19)); #FIDES-5
elif (sys.argv[1] == "SC2000-5"):
  s_box = ((20,26,7,31,19,12,10,15,22,30,13,14, 4,24, 9,18),(27,11, 1,21, 6,16, 2,28,23, 5, 8, 3, 0,17,29,25));

##############################################################################################################
  
######################################## DDT Construction#####################################################   
DDT_SIZE = (len(s_box)*len(s_box[0]))
input_size = int(math.log(DDT_SIZE,2))
DDT = np.zeros( (DDT_SIZE,DDT_SIZE) )
DDT = DDT.astype(int)
sbox_val = []

for p2 in range(DDT_SIZE):
    row = p2 >> 4
    col = p2 & 15
    sbox_val.append(s_box[row][col]);

for p1 in range(DDT_SIZE):
	for p2 in range(DDT_SIZE):
		XOR_IN = np.bitwise_xor(p1,p2);
		XOR_OUT = np.bitwise_xor(sbox_val[p1],sbox_val[p2]);
		DDT[XOR_IN][XOR_OUT] += 1

print(">> DDT of " + str(sys.argv[1]) +" is:" )
print(DDT)
unique_entries = np.unique(DDT)
unique_entries_count = len(np.unique(DDT))
print(">> DDT Entries are: " + str(unique_entries))
print(">> No of Possible Transitions in the DDT: " + str((DDT!=0).sum()))
print(">> No of Impossible Transitions in the DDT: " + str((DDT==0).sum()))
for entry in unique_entries:
    if(entry!=0):
        print(">> No of Possible Transitions in the DDT with Entry "+ str(entry) +": " + str((DDT==entry).sum()))

##############################################################################################################

##############################Generating Possilbe and Impossible Points From DDT################################
        
diff_arr = []
diff_arr_with_1 = []
impossible_diff_arr=[]

if (sys.argv[2] == "sbox"):
    check = False
    prob_point = 0
    if (sys.argv[3] != "-"):
       prob_point = int(sys.argv[3])
       check = True
    for row in range(len(DDT)):
            row_hex = bin(row)[2:].zfill(input_size);
            row_arr = [int(i) for i in row_hex];
            for col in range(len(DDT[row])):
                col_hex = bin(col)[2:].zfill(input_size);
                col_arr = [int(i) for i in col_hex];
                #if(DDT[row][col]!=0):
                    
                if ((DDT[row][col]==prob_point)==check):
                    diff_arr += [row_arr+col_arr];
                    diff_arr_with_1 += [[1]+row_arr+col_arr];
                else:
                    impossible_diff_arr += [[1]+row_arr+col_arr];


if (sys.argv[2] == "prob"):               
    diff_arr = []
    diff_arr_with_1 = []
    impossible_diff_arr=[]
    for row in range(len(DDT)):
            row_bin = bin(row)[2:].zfill(input_size);
            row_arr = [int(i) for i in row_bin];
            for col in range(len(DDT[row])):
                col_bin = bin(col)[2:].zfill(input_size);
                col_arr = [int(i) for i in col_bin];
                if(DDT[row][col]!=0):
                    DDT_bin = "".join(['0']*(unique_entries_count-2))
                    for num in range(1,len(unique_entries)-1):
                        if (DDT[row][col] == unique_entries[num]):
                            DDT_bin = DDT_bin[0:len(DDT_bin) - np.where(unique_entries==unique_entries[num])[0][0]] +  '1' + DDT_bin[len(DDT_bin) - np.where(unique_entries==unique_entries[num])[0][0] +1:] 

                    int_DDT_bin = int(DDT_bin,2);
                    DDT_arr = [int(i) for i in DDT_bin];
                    diff_arr += [row_arr+col_arr+DDT_arr];
                    diff_arr_with_1 += [[1]+row_arr+col_arr+DDT_arr];
                    for k in range(0,pow(2,unique_entries_count-2)):
                        if (k!=int_DDT_bin):
                            im_bin = bin(k)[2:].zfill(unique_entries_count-2);
                            im_arr = [int(i) for i in im_bin];
                            impossible_diff_arr += [[1] + row_arr+col_arr+im_arr];
                else:
                    for k in range(0,pow(2,unique_entries_count-2)):
                            im_bin = bin(k)[2:].zfill(unique_entries_count-2);
                            im_arr = [int(i) for i in im_bin];
                            impossible_diff_arr += [[1]+row_arr+col_arr+im_arr];

#############################################################################################################

###################################Convex Hull############################################################### 
  
ineq_list = []  
vertices = map(array, diff_arr)
A, b = compute_polytope_halfspaces(vertices)
ineq_list = ineq_list + np.column_stack((b,A)).astype(int).tolist()
print(">>> Number of Initial Linear Inequalities: "+ str(len(ineq_list)))

###########################To reduce ch inequalities which do not remove any impossible points################
########################## To get the results equivalent to Boura and Coggia Appraoch this section must be commented#########
ineq_to_remove = []
for i in range(0,len(ineq_list)):
    ineq_solve_count = (np.multiply(np.array(impossible_diff_arr),np.array(ineq_list[i]))).sum(1)
    if (len(np.where(ineq_solve_count<0)[0])==0):
         ineq_to_remove.append(ineq_list[i])

for ineq in ineq_to_remove:
    ineq_list.remove(ineq)
print(">> Number of Inequalities After Prelimnary Reduction: " + str(len(ineq_list)))
################################################################################################################    


if (sys.argv[4]!="-"):
    ineq_impoints_remove = []
    for i in range(0,len(ineq_list)):
        if(i%100==0):
                print(">> Step 1: Inequality Covered: " + str(i))
        ineq_solve_count = (np.multiply(np.array(impossible_diff_arr),np.array(ineq_list[i]))).sum(1)
        ineq_impoints_remove.append([i,len(np.where(ineq_solve_count<0)[0])])
    ineq_impoints_remove.sort(key=lambda x: x[1],reverse=True)
    #print(ineq_impoints_remove)
    ineq_points_add= []
    for i in range(0,len(ineq_list)):
        ineq_solve_count = (np.multiply(np.array(diff_arr_with_1),np.array(ineq_list[i]))).sum(1)
        ineq_points_add.append([i,len(np.where(ineq_solve_count>=0)[0])])
    #print(ineq_points_add)
    final_ineq_list = np.empty((0,len(ineq_list[0])))
    offset = int(sys.argv[4])
    if (offset == 0):
        offset = len(ineq_list)
    if (int(sys.argv[5]) != 1):
        for ineq in range(0,len(ineq_list),offset):
        #for ineq in [0,100,200]:
            if(ineq%100==0):
                print(">> Step 2: Inequality Covered 2: " + str(ineq))
            sub_ineq_list = [ineq_list[x[0]] for x in ineq_impoints_remove[ineq:min(ineq+offset,len(ineq_list))]]
            
            if(int(sys.argv[5])==2):
            
                for i in range(0,len(sub_ineq_list)):
                    vector1 = np.array(sub_ineq_list)
                    vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                    final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2)))
            
            if(int(sys.argv[5])==3):
                for i in range(0,len(sub_ineq_list)):
                    for j in range(0,len(sub_ineq_list)):
                        vector1 = np.array(sub_ineq_list)
                        vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                        vector3 = np.array(sub_ineq_list[j:] + sub_ineq_list[0:j] )
                        final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2 + vector3)))
            
            if(int(sys.argv[5])==4):
                for i in range(0,len(sub_ineq_list)):
                    for j in range(0,len(sub_ineq_list)):
                        for k in range(0,len(sub_ineq_list)):    
                            vector1 = np.array(sub_ineq_list)
                            vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                            vector3 = np.array(sub_ineq_list[j:] + sub_ineq_list[0:j] )
                            vector4 = np.array(sub_ineq_list[k:] + sub_ineq_list[0:k] )
                            final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2 + vector3 + vector4)))
            if(int(sys.argv[5])==5):
                for i in range(0,len(sub_ineq_list)):
                    for j in range(0,len(sub_ineq_list)):
                        for k in range(0,len(sub_ineq_list)):    
                            for l in range(0,len(sub_ineq_list)):  
                                vector1 = np.array(sub_ineq_list)
                                vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                                vector3 = np.array(sub_ineq_list[j:] + sub_ineq_list[0:j] )
                                vector4 = np.array(sub_ineq_list[k:] + sub_ineq_list[0:k] )
                                vector5 = np.array(sub_ineq_list[l:] + sub_ineq_list[0:l] )
                                final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2 + vector3 + vector4 + vector5)))
        
        ineq_list = final_ineq_list
    print(">> Number of Linear Inequalities After Introducing New Inequalities: "+ str(len(ineq_list)))

reduced_ineq_list = MILP_Solve(np.array(ineq_list),np.array(impossible_diff_arr))

print(">> Number of Linear Inequalities After Introducing New Inequalities: "+ str(len(ineq_list)))
print(">> Number of Minimzed Linear Inequalities: " + str(len(reduced_ineq_list)))

import os,sys
import numpy as np

class phonon_1bo_params:
    def __init__(self,alpha,beta):
        self.alpha=alpha
        self.beta=beta
        self.sigma_0=np.array([[1.0,0.0],[0.0,1.0]])
        self.sigma_x=np.array([[0.0,1.0],[1.0,0.0]])
        self.sigma_y=np.array([[0.0,-1j],[1j,0.0]])
        self.sigma_z=np.array([[1.0,0.0],[0.0,-1.0]])

        self.tau_0=np.array([[1.0,0.0],[0.0,1.0]])
        self.tau_x=np.array([[0.0,1.0],[1.0,0.0]])
        self.tau_y=np.array([[0.0,-1j],[1j,0.0]])
        self.tau_z=np.array([[1.0,0.0],[0.0,-1.0]])

        self.reprMat = np.array([[self.alpha,0.0],[0.0,self.beta]])

        self.VfX = np.kron(self.tau_x,self.sigma_x)
        self.VfY = np.kron(self.tau_y,self.sigma_x)
        self.VcX = np.kron(self.tau_x, np.kron(self.reprMat, self.sigma_x))
        self.VcY = np.kron(self.tau_y, np.kron(self.reprMat, self.sigma_x))

class State:                   
    def __init__(self,dagger,elec_type,elec_valley,elec_orb):
        self.dagger = dagger
        self.elec_valley = elec_valley + 1
        self.elec_orb = elec_orb + 1
        self.elec_type = elec_type

class onebody_op:
    def __init__(self,elec_type,coeff_matrix):
        self.elec_type = elec_type
        self.op_array = []
        self.op_dagger_array = []
        self.coeff_matrix = coeff_matrix
        self.op_terms = []
        if self.elec_type == "f":
            for ival in range(2):
                for iorb in range(2):
                    self.op_array.append(State(False,"f",ival,iorb))
                    self.op_dagger_array.append(State(True,"f",ival,iorb))
        else:
            for ival in range(2):
                for iorb in range(4):
                    self.op_array.append(State(False,"c",ival,iorb))
                    self.op_dagger_array.append(State(True,"c",ival,iorb))

        
    def create_1bo(self):
        op_terms=[]
        for ii in range(len(self.op_dagger_array)):
            for jj in range(len(self.op_array)):
                if self.coeff_matrix[ii][jj] != 0:
                    op_terms.append([self.coeff_matrix[ii][jj],self.op_dagger_array[ii],self.op_array[jj]])
        self.op_terms = op_terms
        
    def prettyprint(self):
        for iterm in range(len(self.op_terms)):
            term=self.op_terms[iterm]
            coeff=term[0]
            op_1=term[1]
            op_2=term[2]
            if op_1.dagger==True:
                suffix1="dg"
            else:
                suffix1=""
            if op_2.dagger==True:
                suffix2="dg"
            else:
                suffix2=""
            print(str(coeff) +"  "+ op_1.elec_type+suffix1+"_"+str(op_1.elec_valley)+"_"+str(op_1.elec_orb)+" "\
            +op_2.elec_type+suffix2+"_"+str(op_2.elec_valley)+"_"+str(op_2.elec_orb))

class twobody_op:
    def __init__(self,op_1,op_2,coeff):
        self.coeff = coeff
        self.op_1 = op_1
        self.op_2 = op_2
        self.op_terms = []
        self.mf_decoupled = []
        self.NO_FLAG = False
    
    def create_2bo(self):
        op_terms=[]
        for ii in range(len(self.op_1.op_terms)):
            firstop = self.op_1.op_terms[ii]
            for jj in range(len(self.op_2.op_terms)):
                secondop = self.op_2.op_terms[jj]
                overall_coeff = self.coeff * firstop[0]*secondop[0]
                if overall_coeff != 0:
                    op_terms.append([overall_coeff,firstop[1],firstop[2],secondop[1],secondop[2]])
        self.op_terms = op_terms
        
    def normalorder(self):
        dummylist=[]
        for iterm in range(len(self.op_terms)):
            term=self.op_terms[iterm]
           
            #save components of the term
            coeff=term[0]
            a=term[1]
            b=term[2]
            c=term[3]
            d=term[4]
            if (b.dagger == False) and (c.dagger == True):
                pass
            else:
                print("Either already NO or something wrong")
                return
            
            #is there the mean-field term?
            if (b.elec_valley == c.elec_valley) and\
            (b.elec_orb == c.elec_orb) and\
            (b.elec_type == c.elec_type):
                dummylist.append([coeff,a,d])    
            
            #anticommutation adds -1 to the coeff, swap dagger in internal terms
            coeff = -1.0*coeff
            b.dagger=False
            c.dagger=True

            #does the 4-operator term survive
            if (b.dagger == False) and (d.dagger == False) and\
            (b.elec_valley == d.elec_valley) and\
            (b.elec_orb == d.elec_orb) and\
            (b.elec_type == d.elec_type):
                pass
            elif (a.dagger == True) and (c.dagger == True) and\
            (a.elec_valley == c.elec_valley) and\
            (a.elec_orb == c.elec_orb) and\
            (a.elec_type == c.elec_type):
                pass
            else:
                dummylist.append([coeff,a,c,b,d])
        
    #copyback
        self.op_terms=dummylist[:]
        self.NO_FLAG=True
 
    def decouple(self):
        if self.NO_FLAG==False:
            print("Operator is not normal ordered")
            return
        else:
            for iterm in range(len(self.op_terms)):
                term=self.op_terms[iterm]
               
                #save components of the term
                coeff=term[0]
                a=term[1]
                b=term[2]
                try:
                    c=term[3]
                    d=term[4]
                    
                    # -b+d<a+c>
                    self.mf_decoupled.append([-1.0*coeff,b,d,a,c])
                    # -a+c<b+d>
                    self.mf_decoupled.append([-1.0*coeff,a,c,b,d])
                    # b+c<a+d>
                    self.mf_decoupled.append([coeff,b,c,a,d])
                    # a+d<b+c>
                    self.mf_decoupled.append([coeff,a,d,b,c])
                #except it's already mean-field
                except:
                    # -a+b
                    self.mf_decoupled.append([coeff,a,b])

        
    
    def prettyprint(self):
        for iterm in range(len(self.op_terms)):
            term=self.op_terms[iterm]
            coeff=term[0]
            op_1=term[1]
            op_2=term[2]
            try:
                op_3=term[3]
                op_4=term[4]
                if op_1.dagger==True:
                    suffix1="dg"
                else:
                    suffix1=""
                if op_2.dagger==True:
                    suffix2="dg"
                else:
                    suffix2=""
                if op_3.dagger==True:
                    suffix3="dg"
                else:
                    suffix3=""
                if op_4.dagger==True:
                    suffix4="dg"
                else:
                    suffix4=""
                print(str(coeff) +"  "\
                + op_1.elec_type+suffix1+"_"+str(op_1.elec_valley)+"_"+str(op_1.elec_orb)+" "\
                +op_2.elec_type+suffix2+"_"+str(op_2.elec_valley)+"_"+str(op_2.elec_orb)+" "\
                +op_3.elec_type+suffix3+"_"+str(op_3.elec_valley)+"_"+str(op_3.elec_orb)+" "\
                +op_4.elec_type+suffix4+"_"+str(op_4.elec_valley)+"_"+str(op_4.elec_orb))
            except:
                if op_1.dagger==True:
                    suffix1="dg"
                else:
                    suffix1=""
                if op_2.dagger==True:
                    suffix2="dg"
                else:
                    suffix2=""
                print(str(coeff) +"  "\
                + op_1.elec_type+suffix1+"_"+str(op_1.elec_valley)+"_"+str(op_1.elec_orb)+" "\
                +op_2.elec_type+suffix2+"_"+str(op_2.elec_valley)+"_"+str(op_2.elec_orb))
                
    def prettyprint_decoupled(self):
        for iterm in range(len(self.mf_decoupled)):
            term=self.mf_decoupled[iterm]
            coeff=term[0]
            op_1=term[1]
            op_2=term[2]
            #print which is dagger
            if op_1.dagger==True:
                suffix1="dg"
            else:
                suffix1=""
            if op_2.dagger==True:
                suffix2="dg"
            else:
                suffix2=""
            #if there is a density term
            try:
                dens_1=term[3]
                dens_2=term[4]
                densStr="dens_"+dens_1.elec_type+"_"+str(dens_1.elec_valley)+"_"+str(dens_1.elec_orb)+"_"+dens_2.elec_type+"_"+str(dens_2.elec_valley)+"_"+str(dens_2.elec_orb)
            except:
                densStr=""
            print(str(coeff) +"  "\
            +densStr+"  "
            +op_1.elec_type+suffix1+"_"+str(op_1.elec_valley)+"_"+str(op_1.elec_orb)+" "\
            +op_2.elec_type+suffix2+"_"+str(op_2.elec_valley)+"_"+str(op_2.elec_orb))

    def print_to_dc_decoupled(self):
        self.operator_array_for_dc=[]
        self.density_array_for_dc=[]
        self.coeff_for_dc=[]
        for iterm in range(len(self.mf_decoupled)):
            term=self.mf_decoupled[iterm]
            
            coeff=term[0]
            
            op_1_type=term[1].elec_type
            op_1_valley=term[1].elec_valley -1
            op_1_orb=term[1].elec_orb -1
            
            op_2_type=term[2].elec_type
            op_2_valley=term[2].elec_valley -1
            op_2_orb=term[2].elec_orb -1
            try:
                dens_1_type=term[3].elec_type
                dens_1_valley=term[3].elec_valley -1
                dens_1_orb=term[3].elec_orb -1
                
                dens_2_type=term[4].elec_type
                dens_2_valley=term[4].elec_valley -1
                dens_2_orb=term[4].elec_orb -1
                
                self.coeff_for_dc.append(coeff)
                self.operator_array_for_dc.append([op_1_orb,op_1_valley,op_1_type,op_2_orb,op_2_valley,op_2_type])
                self.density_array_for_dc.append([dens_1_orb,dens_1_valley,dens_1_type,dens_2_orb,dens_2_valley,dens_2_type])
            except:
                self.coeff_for_dc.append(coeff)
                self.operator_array_for_dc.append([op_1_orb,op_1_valley,op_1_type,op_2_orb,op_2_valley,op_2_type])
                self.density_array_for_dc.append(["NODENS"])
        
        return self.coeff_for_dc,self.operator_array_for_dc,self.density_array_for_dc





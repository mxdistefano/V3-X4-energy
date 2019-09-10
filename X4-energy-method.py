#!/usr/bin/python

import sys
import os

import yaml

from modeller import *
from modeller.automodel import *

from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.PDBParser import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1


#############################################################################
####   CONFIG   #############################################################
#############################################################################

Nmodels = 9
cutoff = 30.0
ruleresidues = { 11: 1.0, 22: 1.0, 25: 1.0 }

#############################################################################
####  GLOBALS   #############################################################
#############################################################################

with open("Data/distances_not_strict.yml", 'r') as f:
    distancias_no_estrictas = yaml.load(f)
with open("Data/distances_strict.yml", 'r') as g:
    distancias_estrictas = yaml.load(g)
with open("Data/energy.yml", 'r') as e:
    energia = yaml.load(e)

PIRtemplateX4 = """>P1;%s
sequence:%s:0::0::::-1.00: -1.00
MEGISIYTSDNYTEEMGSGDYDSMKEPCFREENANFNKIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAVANWYFGNFLCKAVHVIYTVNLYSSVWILAFISLDRYLAIVHATNSQRPRKLLAEKVVYVGVWIPALLLTIPDFIFANVSEADDRYICDRFYPNDLWVVVFQFQHIMVGLILPGIVILSCYCIIISKLSHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLEIIKQGCEFENTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKGKRGGHSSVSTESESSSFHSS/%s*
>P1;4CRX-%d
structure:4CRX-%d:1:A:35:B:::-1.00: -1.00
MEGISI.TSDN.TEEMGSGD.DSMKEPCFREENANFNKIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAVANWYFGNFLCKAVHVIYTVNLYSSVWILAFISLDRYLAIVHATNSQRPRKLLAEKVVYVGVWIPALLLTIPDFIFANVSEADDRYICDRFYPNDLWVVVFQFQHIMVGLILPGIVILSCYCIIISKLSHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLEIIKQGCEFENTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKGKRGGHSSVSTESESSSFHSS/CTRPNNNTRKRVSLGPGRVWYTTGQIVGDIRKAHC*"""

#############################################################################
####  METHOD ################################################################
#############################################################################

class EnergyMedia:
    def __init__(self,energia,cantContactos,aa):
        self.energia=energia
        self.cantContactos=cantContactos
        self.aa=aa

class Conctact:
    def __init__(self, atom1, atom2, energy):
        self.atom1=atom1
        self.atom2=atom2
        self.energy=energy
    
    def __eq__(self, listOther):
        for other in listOther:
            if (self.atom1 == other.atom1 and self.atom2 == other.atom2):
                return False
        return True

class ContactV3:
    def __init__(self, atom1, atom2, distance,energy):
        self.atomV3=atom2
        self.atomCR=atom1
        self.distance=distance
        self.energy = energy

def calcCOM(res):
    com = [0 for i in range(3)]
    atomCount = 0
    for a in res:
        com[0] = float(a.coord[0]) + com[0] 
        com[1] = float(a.coord[1]) + com[1]
        com[2] = float(a.coord[2]) + com[2]
        atomCount += 1
    return [i / atomCount for i in com] 
        
def distanceBetweenCOM(residue1,residue2):
    dist = 0
    com1 = calcCOM(residue1)
    com2 = calcCOM(residue2)
    for a_i, b_i in zip(com1, com2):
        dist += ((a_i - b_i)**2)
    return dist ** 0.5

def avgenergy(seq, N):
    p = PDBParser(QUIET=True)
    a="%s_%d.BL000%d0001.pdb"
    avg = 0.0
    positions = {}
    for t in range(1,11):
        for i in range(1,N+1):
            aa = a % ( seq, t, i )
            try:
                s = p.get_structure(aa, aa)
            except:
                continue
 
            ppb=PPBuilder()
 
            chains = s[0].get_list()
 
            ccr5 = chains[0]
            gp120 = chains[1]
 
            total_energy = 0
            for r1 in gp120:
                r1Code= str(protein_letters_3to1.get(r1.resname.title()))
                if (r1.get_id()[1]-352) not in ruleresidues: continue
                for r2 in ccr5:
                    k = distanceBetweenCOM(r1,r2)
                    r2Code= str(protein_letters_3to1.get(r2.resname.title()))
                    try:
                        cutoff_estricto = distancias_estrictas["distances"][r1Code+"_"+r2Code]
                        cutoff_no_estricto = distancias_no_estrictas["distances"][r1Code+"_"+r2Code]
                    except:
                        cutoff_estricto = distancias_estrictas["distances"][r2Code+"_"+r1Code]
                        cutoff_no_estricto = distancias_no_estrictas["distances"][r2Code+"_"+r1Code]

                    if k < cutoff_no_estricto:
                        try:
                            en = energia["energy"][r1Code+"_"+r2Code]
                        except:
                            en = energia["energy"][r2Code+"_"+r1Code]

			if (r1.get_id()[1]-352) not in positions: positions[r1.get_id()[1]-352] = []
                        positions[r1.get_id()[1]-352].append(float(en))
 
                        total_energy += float(en) * ruleresidues[r1.get_id()[1]-352]
 
            avg+= total_energy

    return -avg/(N*10)

class MyModel(loopmodel):
    def special_restraints(self, aln):
        for at in selection(self.chains['A']):
            if at.name=="CA":
                self.restraints.add(forms.gaussian(group=physical.xy_distance,feature=features.x_coordinate(at),mean=at._Atom__get_x(),stdev=0.1))
                self.restraints.add(forms.gaussian(group=physical.xy_distance,feature=features.y_coordinate(at),mean=at._Atom__get_y(),stdev=0.1))
                self.restraints.add(forms.gaussian(group=physical.xy_distance,feature=features.z_coordinate(at),mean=at._Atom__get_z(),stdev=0.1))
        for at in selection(self.chains['B']):
            if at.name=="CA" and 353 <= int(at.residue.num) <= 387:
                self.restraints.add(forms.gaussian(group=physical.xy_distance,feature=features.x_coordinate(at),mean=at._Atom__get_x(),stdev=0.1))
                self.restraints.add(forms.gaussian(group=physical.xy_distance,feature=features.y_coordinate(at),mean=at._Atom__get_y(),stdev=0.1))
                self.restraints.add(forms.gaussian(group=physical.xy_distance,feature=features.z_coordinate(at),mean=at._Atom__get_z(),stdev=0.1))

    def select_loop_atoms(self):
        sel = [ at for at in self.chains["B"].atoms if 353 < int(at.residue.num) < 387 ]
        return selection(sel)

def MakeModels(pirfile, template, seqID, n):
    env = environ()
    env.io.atom_files_directory = ['.', 'Templates/']
    a = MyModel(env, alnfile=pirfile, knowns=template, sequence=seqID, assess_methods=assess.DOPE)
    a.starting_model = 1
    a.ending_model = 1
    a.loop.starting_model = n
    a.loop.ending_model = n
    a.md_level = refine.fast
    a.make()


#############################################################################
####  MAIN   ################################################################
#############################################################################

if len(sys.argv)<2:
    print "Usage:"
    print "\tpython "+ sys.argv[0] + "<V3 loop sequence>"
    print "Example: python "+ sys.argv[0] + "CTRPNNNTRKRVSLGPGRVWYTTGQIVGDIRKAHC"
    exit()

seq = sys.argv[1]

for i in range(1, 11):
	with open("%s_%s.pir" % (seq, str(i)), "w") as f:
		print >> f, PIRtemplateX4 % ("%s_%d" % (seq, i), "%s_%d" % (seq, i), seq, i, i)

	for j in range(1,Nmodels+1):
		if not os.path.exists( "%s_%d.BL000%d0001.pdb" % (seq, i, j)  ):
 			MakeModels("%s_%d.pir" % (seq, i), "4CRX-%d" % i, "%s_%d" % (seq, i), j)

predictionScore = avgenergy(seq, Nmodels)

print "V3 loop sequence: ", seq
print "Predicted Tropism: ", "X4" if predictionScore > cutoff else "R5"
print "X4-Energy Score: ", predictionScore


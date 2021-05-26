import re
import tkinter as tk
import os.path
from matplotlib import pyplot as plt
from matplotlib import style
import numpy as np
from tkinter import messagebox
from pathlib import Path

#GlobalVariable to store names of all the Different Genome in the raw text file
Filenamelist=[]

#The genome class 
class genome(object):  
  
  #Object constructor
  def __init__(self,inp):

  	#if else block for constructor overloading
    
    if len(inp)==2:  #If reading from the raw file (First read)
      self.Name=inp[0]                                                                                        #stores the Name
      self.DNAseq=inp[1]																					  #stores the DNA Sequence
      self.whaterror={0:'No Error',1:"'Stop' In Middle of sequence",2:'DNA Sequence Not A Multiple of Three'} #An error dict used for storing DNA Sequence Errors
      self.CodonCountDict={                   #Initializing Dictionary map to store Codon counts to calculate required values  
        'ATA':0, 'ATC':0, 'ATT':0, 'ATG':0, 
        'ACA':0, 'ACC':0, 'ACG':0, 'ACT':0, 
        'AAC':0, 'AAT':0, 'AAA':0, 'AAG':0, 
        'AGC':0, 'AGT':0, 'AGA':0, 'AGG':0,                  
        'CTA':0, 'CTC':0, 'CTG':0, 'CTT':0, 
        'CCA':0, 'CCC':0, 'CCG':0, 'CCT':0, 
        'CAC':0, 'CAT':0, 'CAA':0, 'CAG':0, 
        'CGA':0, 'CGC':0, 'CGG':0, 'CGT':0, 
        'GTA':0, 'GTC':0, 'GTG':0, 'GTT':0, 
        'GCA':0, 'GCC':0, 'GCG':0, 'GCT':0, 
        'GAC':0, 'GAT':0, 'GAA':0, 'GAG':0, 
        'GGA':0, 'GGC':0, 'GGG':0, 'GGT':0, 
        'TCA':0, 'TCC':0, 'TCG':0, 'TCT':0, 
        'TTC':0, 'TTT':0, 'TTA':0, 'TTG':0, 
        'TAC':0, 'TAT':0, 'TAA':0, 'TAG':0, 
        'TGC':0, 'TGT':0, 'TGA':0, 'TGG':0,}

      self.ProtienSeq=self.translate()  #calls the member function transate() to get the protien seq and store it in a instance variable
      self.ProtienCodonCount=self.GetCodonDistribution() #get the codon count distribution
      if self.errorno==0:    #check if any errors were set during protien sequence generation if yes set Nc to 0 else
      	self.Nc=self.calcNcValue()    #Compute Nc Value
      else:
      	self.Nc=0
    else:             #if reading from a processed file/db table get individual info and discard unnecessary info
      self.Name=inp[0]
      self.DNAseq=inp[1]
      self.ProtienSeq=inp[2]
      self.Nc=float(inp[4])
      self.errorno=int(inp[3])
      self.whaterror={0:'No Error',1:"'Stop' In Middle of sequence",2:'DNA Sequence Not A Multiple of Three'}
  
  #function to generate the protien sequence from DNA sequence

  def translate(self):       
    table = {                                         #A map to translate individual codons into corresponding protien
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'Stop', 'TAG':'Stop', 
        'TGC':'C', 'TGT':'C', 'TGA':'Stop', 'TGG':'W', 
    } 
    protein =""                                       #temp variable to store the protien seq
    length=len(self.DNAseq)                           #check DNA seq Length for errors
    if length%3 == 0: 
      for i in range(0, length, 3):                   #take set of three genes from the start and translate it to corresponding protien
        codon = self.DNAseq[i:i + 3] 
        if  i<(length-3) and table[codon] == 'Stop' : #check of "stop in middle" error
          self.errorno=1
          return 'N/A' 
        protein+= table[codon] 
        self.storeCodonCount(codon)
    else:
      self.errorno=2
      return 'N/A'
    self.errorno=0
    return protein
  
  #similar to tostring() in java, a string representation of an object
  def __str__(self):                  
    return ('\nGenome Name : ' + self.Name + '\nGenome DNA sequence :\n' + self.DNAseq + '\nGenome Protien Sequence :\n'+ self.ProtienSeq + '\nError if any : ' + self.whaterror[self.errorno] + '\nNc Value of the Codon : ' + str(self.Nc))
  
  #increment the codon counts of each gene as an when the occur
  def storeCodonCount(self,codname):      
    if (codname in  self.CodonCountDict):
        self.CodonCountDict[codname]+=1

  #Function to get the codon counts (number only) grouped by the protien they generate		
  def GetCodonDistribution(self):       
    self.CodonDist=[
    [self.CodonCountDict['ATG']],
    [self.CodonCountDict['TGG']],
    [self.CodonCountDict['TTT'],self.CodonCountDict['TTC']],
    [self.CodonCountDict['TAT'],self.CodonCountDict['TAC']],
    [self.CodonCountDict['CAT'],self.CodonCountDict['CAC']],
    [self.CodonCountDict['CAA'],self.CodonCountDict['CAG']],
    [self.CodonCountDict['AAT'],self.CodonCountDict['AAC']],
    [self.CodonCountDict['AAA'],self.CodonCountDict['AAG']],
    [self.CodonCountDict['GAT'],self.CodonCountDict['GAC']],
    [self.CodonCountDict['GAA'],self.CodonCountDict['GAG']],
    [self.CodonCountDict['TGT'],self.CodonCountDict['TGC']],
    [self.CodonCountDict['ATT'],self.CodonCountDict['ATC'],self.CodonCountDict['ATA']],
    [self.CodonCountDict['GTT'],self.CodonCountDict['GTC'],self.CodonCountDict['GTA'],self.CodonCountDict['GTG']],
    [self.CodonCountDict['CCT'],self.CodonCountDict['CCC'],self.CodonCountDict['CCA'],self.CodonCountDict['CCG']],
    [self.CodonCountDict['ACT'],self.CodonCountDict['ACC'],self.CodonCountDict['ACA'],self.CodonCountDict['ACG']],
    [self.CodonCountDict['GCT'],self.CodonCountDict['GCC'],self.CodonCountDict['GCA'],self.CodonCountDict['GCG']],
    [self.CodonCountDict['GGT'],self.CodonCountDict['GGC'],self.CodonCountDict['GGA'],self.CodonCountDict['GGG']],
    [self.CodonCountDict['TTA'],self.CodonCountDict['TTG'],self.CodonCountDict['CTT'],self.CodonCountDict['CTC'],self.CodonCountDict['CTA'],self.CodonCountDict['CTG']],
    [self.CodonCountDict['TCT'],self.CodonCountDict['TCC'],self.CodonCountDict['TCA'],self.CodonCountDict['TCG'],self.CodonCountDict['AGT'],self.CodonCountDict['AGC']],
    [self.CodonCountDict['CGT'],self.CodonCountDict['CGC'],self.CodonCountDict['CGA'],self.CodonCountDict['CGG'],self.CodonCountDict['AGA'],self.CodonCountDict['AGG']]
    ]
    #map to store the codon counts grouped by the protiens they generate
    ProtienDist={
    'Met':self.CodonDist[0],   #eg this is a 'list' of counts of all individual codons that makes up the protien Met
    'Trp':self.CodonDist[1],
    'Phe':self.CodonDist[2],
    'Tyr':self.CodonDist[3],
    'His':self.CodonDist[4],
    'Gln':self.CodonDist[5],
    'Asn':self.CodonDist[6],
    'Lys':self.CodonDist[7],
    'Asp':self.CodonDist[8],
    'Glu':self.CodonDist[9],
    'Cys':self.CodonDist[10],
    'Lle':self.CodonDist[11],
    'Val':self.CodonDist[12],
    'Pro':self.CodonDist[13],
    'Thr':self.CodonDist[14],
    'Ala':self.CodonDist[15],
    'Gly':self.CodonDist[16],
    'Leu':self.CodonDist[17],
    'Ser':self.CodonDist[18],
    'Arg':self.CodonDist[19]}
    return ProtienDist

  #get the Pi squared value of the different individual 'types'of same protien,
  def getpisqvec(self):
    subtotal=0
    temppi2vec=[]
    for i in range(len(self.CodonDist)):
      temp=[]
      for j in range(len(self.CodonDist[i])):
        subtotal+=self.CodonDist[i][j]
      for j in range(len(self.CodonDist[i])):
        if subtotal == 0 or self.CodonDist[i][j]== 0:
          temp.append(0)
        else:
          temp.append((self.CodonDist[i][j]/subtotal)**2)
      temppi2vec.append(temp)
      subtotal=0
    return temppi2vec
  
  #generate a list of Fk values of each protien
  def getfk(self):
    subtotal=0
    pi2sumvec=[]
    pi2vec=self.getpisqvec()
    for i in range(len(pi2vec)):
      for j in range(len(pi2vec[i])):
        subtotal+=pi2vec[i][j]
      if subtotal == 0:
        pi2sumvec.append(0)
      else:
        pi2sumvec.append(subtotal)
      subtotal=0
    for i in range(len(pi2sumvec)):
      if pi2sumvec[i]!=0:
        pi2sumvec[i]=1/pi2sumvec[i]
    return pi2sumvec
  
  #compute Nc value of this genome instance (sum of fk vlaues)
  def  calcNcValue(self):
    subtotal=0
    pi2sumvec=self.getfk()
    fk=0
    for i in range(len(pi2sumvec)):
      fk+=pi2sumvec[i]
    return fk
  #function to write the compute heavy processed data to a file (later to be upgraded to work with a database backend)
  def writetofile(self):
    if (not os.path.isdir('GenomeObjFiles')):     #checks if output folder exists and makes one if it doesn't
      os.mkdir('GenomeObjFiles')
    filetoopen=os.path.join('GenomeObjFiles' , self.Name[-11:]) 
    filepath=Path(filetoopen+'.txt')
    #if filepath.is_file():
      #return
    fileopen=open(filetoopen + ".txt", 'w')
    fileopen.write(self.Name + '\n')
    fileopen.write(self.DNAseq + '\n')
    fileopen.write(self.ProtienSeq + '\n')
    fileopen.write(str(self.errorno) +'\n')
    fileopen.write(str(self.Nc) + '\n')
    fileopen.close


class fullGenome(object):

  def __init__(self,fname):
    dna_str=[]
    self.check = 0
    self.No_A=0
    self.No_T=0
    self.No_G=0
    self.No_C=0
    self.ATarr=[]
    self.GCarr=[]
    self.ATcumu=[]
    self.GCcumu=[]
    self.ATarr.append(0)
    self.GCarr.append(0)
    self.ATcumu.append(0)
    self.GCcumu.append(0)
    dnaArr=[]
    seqstr=''
    eofcheck=False
    my_file = Path(fname+".txt")
    if not my_file.is_file():
      messagebox.showinfo("ERROR!!!", "Wrong Input File Name Given, Please Give correct Filename")
    with open(fname + ".txt",'r') as dna:
      if self.check==0:
        dna_str=dna.readline().replace("\n","")
        dna1st=''.join(str(e) for e in dna_str)
        dna1stS=dna1st.lower()
        self.check=1
        if not dna1stS[-18:] == "complete sequence.":
          messagebox.showinfo("ERROR!!!", "Wrong type of file uploaded, Please Give correct File")
          return
      while eofcheck==False:
        for var in range(18):
          dna_str=dna.readline().replace("\n","")
          if dna_str=="":
            eofcheck=True
          dna_str = dna_str.replace("\r", "")
          dna_str = dna_str.replace("u", "t")
          dnaArr.append(dna_str)
        seqstr = "".join(str(e) for e in dnaArr)
        seqstr2 = seqstr.upper()
        self.calc_skew(seqstr2)
        dnaArr=[]
    self.get_skewPlots()

  def calc_skew(self,dnastr):
    A=self.No_A
    T=self.No_T
    C=self.No_C
    G=self.No_G
    for i in dnastr:
      if i == 'A':
        A=A+1
      elif i == 'T':
        T=T+1
      elif i == 'G':
        G=G+1
      else:
        C=C+1
    self.No_A=A
    self.No_T=T
    self.No_G=G
    self.No_C=C

    self.ATarr.append((self.No_A-self.No_T))
    self.GCarr.append((self.No_G-self.No_C))
    #self.AT.append(self.No_A-self.No_T)
    #self.GC.append(self.No_G-self.No_C)

  def get_skewPlots(self):
    style.use('ggplot')
    x1 = [int(i*1080) for i in range(len(self.ATarr))]
    x2 = [int(i*1080) for i in range(len(self.GCarr))]
    y1=[self.ATarr[int(i/1080)] for i in x1]
    y2=[self.GCarr[int(i/1080)] for i in x2]
    #y3=[self.AT[i] for i in x]
    plt.subplot(2,1,1)
    plt.plot(x1,y1,'g',label='AT-Skew', linewidth=1)
    plt.subplot(2,1,2)
    plt.plot(x2,y2,'c',label='GC-Skew',linewidth=1)
    #plt.plot(x,y3,'r',label='A-T cumulative', linewidth=1)
    plt.title('Cumulative DNA Skew Graph')
    plt.ylabel('Skewness')
    plt.xlabel('No of Nucleotides In the DNA')
    plt.legend()
    plt.grid(True,color='k')
    #plt.savefig('OutputGraphs\\AT and GC Skew_' + self.Name[-11:-1] + '.png')
    plt.show()
  

def readfromfile(Name):
  if (not os.path.isdir('GenomeObjFiles')):
      messagebox.showinfo('ERROR', 'NO GENERATED FILES FOUND, PLEASE UPLOAD FILE FIRST')
      return
  filetoopen=os.path.join('GenomeObjFiles' ,Name[-11:])
  inpfile=open(filetoopen + '.txt' , 'r')
  Genomename=inpfile.readline()[:-1]
  GenomeDnaseq=inpfile.readline()[:-1]
  GenomeProtienseq=inpfile.readline()[:-1]
  Genomeerrno=inpfile.readline()[:-1]
  GenomeNcVal=inpfile.readline()[:-1]
  inpvec=[Genomename,GenomeDnaseq,GenomeProtienseq,Genomeerrno,GenomeNcVal]
  TempGenome=genome(inpvec)
  inpfile.close()
  return TempGenome
  
  

def getDnaList(fname):
  dna_str=[]
  matchObj=re.compile(r'>\w{4}_\w{2}\d{4}:\|\d{1,4}\|\w{6}_\d{4}')
  my_file = Path(fname+".txt")
  if not my_file.is_file():
    messagebox.showinfo("ERROR!!!", "Wrong Input File Name Given, Please Give correct Filename")
  with open(fname + ".txt",'r') as dna:
    dna_str=dna.read().replace("\n","")
  dna_str = dna_str.replace("\r", "")
  dna_str = dna_str.replace("u", "t")
  seqstr = ''.join(str(e) for e in dna_str)
  seqstr2 = seqstr.upper()
  MatchList=matchObj.finditer(seqstr2)
  span_list=[]
  for match in MatchList:
    span_list.append(list(match.span()))
  if span_list==[]:
  	messagebox.showinfo("ERROR!!!", "Wrong Type of file, Please use correct Format for the file")
  	return 
  dna_list=[]
  last=0
  for var in range(len(span_list)-1):
    dna_list.append(seqstr2[span_list[var][0]:span_list[var][1]])
    dna_list.append(seqstr2[span_list[var][1]:span_list[var+1][0]])
    last=var+1
  dna_list.append(seqstr2[span_list[last][0]:span_list[last][1]])
  dna_list.append(seqstr2[span_list[last][1]:])
  return dna_list
	
def preprocessfile():
  fname=filenameVar.get()
  dna_list=getDnaList(fname)
  for var in range(0,len(dna_list)-1,2):
    inplist=[dna_list[var],dna_list[var+1]]
    GenomeObj=genome(inplist)
    GenomeObj.writetofile()
    if var%2==0:
      Filenamelist.append(dna_list[var])
  #insert code to remove previous files entries
  dropDownMenu.delete(0, 'end')
  for var in range(len(Filenamelist)):
    dropDownMenu.insert(var,Filenamelist[var])
  del Filenamelist[:]

def printtable():
  data=gendata.get()
  headers=['Genome Name','DNA Sequence','Generated Protien Sequence', 'Input error if any', 'Nc Value']
      

def fetchShow():
  if dropDownMenu.get(tk.ACTIVE)=="":
    messagebox.showinfo("ERROR!!!", "No Gene Sequence Selected. Please select gene sequence or upload the genome file first.")
    return
  tempvar=dropDownMenu.get(tk.ACTIVE)
  gen=readfromfile(tempvar)
  Outputbox.config(state="normal")
  Outputbox.delete('1.0', tk.END)
  Outputbox.insert(tk.INSERT,"Name:\n")
  Outputbox.insert(tk.INSERT,gen.Name)
  Outputbox.insert(tk.INSERT,"\nDNA sequence:\n")
  Outputbox.insert(tk.INSERT,gen.DNAseq)
  Outputbox.insert(tk.INSERT,"\nProtien Sequence:\n")
  Outputbox.insert(tk.INSERT,gen.ProtienSeq)
  Outputbox.insert(tk.INSERT,"\nInput Error if any:\n")
  Outputbox.insert(tk.INSERT,gen.whaterror[gen.errorno])
  Outputbox.insert(tk.INSERT,"\nNc Value:\n")
  Outputbox.insert(tk.INSERT,gen.Nc)
  Outputbox.config(state="disabled")

def returnShow():
  tempvar=dropDownMenu.get(tk.ACTIVE)
  gen=readfromfile(tempvar)
  return gen
  
def genGraph():
  if filenamePlotVar.get()=="":
    messagebox.showinfo("ERROR!!!", "Filename Doesn't match with any file in the local directory. Please Check if you entered the correct Filename and it is present in the current Directory.")
    return
  Com_gene_filename=filenamePlotVar.get()
  fullGen=fullGenome(Com_gene_filename)


if __name__=="__main__":
  root=tk.Tk()
  root.geometry("1200x1000")

  root.title('DNA sequence Processor')
  
  frame1=tk.Frame(root)
  l=tk.Label(frame1,text="Upload the DNA sequence File Name :")
  l.grid(row=0,column=0,padx=10,pady=10)
  filenameVar=tk.StringVar(frame1)
  entryb=tk.Entry(frame1,textvariable=filenameVar).grid(row=0,column=1,padx=10,pady=10)
  fileb=tk.Button(frame1,text= "Process",command=preprocessfile)
  fileb.grid(row=0,column=2,padx=10,pady=10)
  frame1.grid(row=0,column=0,padx=10,pady=10)
  
  frame2=tk.Frame(root)  
  l2=tk.Label(frame2,text="Choose which file info is needed")
  l2.grid(row=0,column=0,padx=10,pady=10)
  fr=tk.Frame(frame2)
  dropDownMenu = tk.Listbox(fr,font=("Verdana",10),width=100)
  dropDownMenu.pack(side='left',fill="both",expand=True)
  sbr=tk.Scrollbar(fr)
  sbr.pack(side='right',fill='y')
  sbr.config(command=dropDownMenu.yview)
  dropDownMenu.config(yscrollcommand=sbr.set)
  fr.grid(row=1,column=0,padx=10,pady=10)
  frame2.grid(row=1,column=0,padx=10,pady=10)

  frame3=tk.Frame(root)
  dispbtn=tk.Button(frame3,text="Get Details",command=fetchShow)
  dispbtn.grid(row=0,column=1,padx=10,pady=10)
  btnl=tk.Label(frame3,text="Get The Details of Selected Genome :")
  btnl.grid(row=0,column=0,padx=10,pady=10)
  frame3.grid(row=2,column=0)

  frame4=tk.Frame(root)
  l3=tk.Label(frame4,text="Details of the selected Geneome Data :")
  l3.grid(row=0,column=0,padx=10,pady=10)
  fr2=tk.Frame(frame4)
  Outputbox = tk.Text(fr2,font=("Verdana",10),height=15,width=100)
  Outputbox.pack(side='left',fill="both",expand=True)
  sbr2=tk.Scrollbar(fr2)
  sbr2.pack(side='right',fill='y')
  sbr2.config(command=Outputbox.yview)
  Outputbox.config(yscrollcommand=sbr2.set)
  fr2.grid(row=1,column=0,padx=10,pady=10)
  frame4.grid(row=3,column=0,padx=10,pady=10)

  frame5=tk.Frame(root)
  l2=tk.Label(frame5,text="Give the Name of the File with the Complete genome Sequence: ")
  l2.grid(row=0,column=0,padx=10,pady=10)
  filenamePlotVar=tk.StringVar(frame5)
  entryb=tk.Entry(frame5,textvariable=filenamePlotVar).grid(row=0,column=1,padx=10,pady=10)
  graphbtn=tk.Button(frame5,text = "Generate Graph",command=genGraph)
  graphbtn.grid(row=0,column=2,padx=10,pady=10)
  frame5.grid(row=4,column=0,padx=10,pady=10)
  root.mainloop()

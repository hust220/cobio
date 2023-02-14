import tarfile
import sys
import cobio as jp

proteins_tar = 'proteins.tar'
def get_protein_sequences(outfile):                                                                                                                                              
  t = tarfile.open(f'{proteins_tar}', mode='r')                                                                                                                                  
  f = open(outfile, 'w+')                                                                                                                                                        
  for m in t.getmembers():                                                                                                                                                       
    print(m.name)                                                                                                                                                                
    sys.stdout.flush()                                                                                                                                                           
    pdbfile = t.extractfile(m.name)                                                                                                                                              
    pdb = jp.Pdb(pdbfile)                                                                                                                                                        
    seq = pdb.seq()                                                                                                                                                              
    f.write(f'>{m.name}\n{seq}\n')                                                                                                                                                 
#    pdbfile.close()                                                                                                                                                             
  f.close()                                                                                                                                                                      
  t.close()

if __name__ == '__main__':
  get_protein_sequences(sys.argv[1])                                                                                                                                                                     
        
